// *********************************** //
// ******* Process Definitions ******* //
// *********************************** //

// Write process directive (https://www.nextflow.io/docs/latest/dsl2.html#process)
process guppy_basecaller {
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/summary/*' // "**" means include subdirecteries in glob pattern"
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/bams/*'
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/fastq/*'
    container  "dhspence/docker-guppy"
    cpus 6
    memory '12 GB'
    
    input:
        path fast5_path
    output:
        path "basecall/summary/*", emit: summary
        path "basecall/bams", emit: basecall_bams_path
        path "basecall/fastq/*", emit: fastqs

    script:
        """
        guppy_basecaller -i $fast5_path --bam_out -s unaligned_bam -c /opt/ont/guppy/data/${params.basecall_config} --num_callers ${task.cpus}
        
        cat unaligned_bam/pass/*.fastq > unaligned_bam/pass/${params.sample_id}.fastq
        gzip unaligned_bam/pass/${params.sample_id}.fastq

        mkdir basecall
        mkdir basecall/summary && mv unaligned_bam/sequencing_* basecall/summary
        mkdir basecall/bams && mv unaligned_bam/pass/*.bam basecall/bams
        mkdir basecall/fastq && mv unaligned_bam/pass/${params.sample_id}.fastq.gz basecall/fastq/${params.sample_id}.fastq.gz
        rm -rf unaligned_bam

        """
}

process guppy_basecaller_gpu {
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/summary/*' // "**" means include subdirecteries in glob pattern"
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/fastq/*'
    container  "dhspence/docker-gguppy"
    cpus 2
    memory '32 GB'
    
    input:
        path fast5_path
    output:
        path "basecall/summary/*", emit: summary
        path "basecall/bams", emit: basecall_bams_path
        path "basecall/fastq/*", emit: fastqs

    script:
        """
        guppy_basecaller -i $fast5_path --bam_out -s unaligned_bam -c /opt/ont/guppy/data/${params.basecall_config} --num_callers ${task.cpus} -x cuda:all:100%
        
        cat unaligned_bam/pass/*.fastq > unaligned_bam/pass/${params.sample_id}.fastq
        gzip unaligned_bam/pass/${params.sample_id}.fastq

        mkdir basecall
        mkdir basecall/summary && mv unaligned_bam/sequencing_* basecall/summary
        mkdir basecall/bams && mv unaligned_bam/pass/*.bam basecall/bams
        mkdir basecall/fastq && mv unaligned_bam/pass/${params.sample_id}.fastq.gz basecall/fastq/${params.sample_id}.fastq.gz
        rm -rf unaligned_bam

        """

}


process guppy_aligner {
    publishDir "${params.outdir}", mode:'copy', pattern: 'alignment/**', enabled: false // disabled by default. These bam files are too big, and better to save the merged sorted bam from the next process (merge_bams)
    container  "dhspence/docker-guppy"
    cpus 6
    memory '12 GB'
    
    input:
        path unaligned_bams
        path reference_fasta
    output:
        path "alignment/*.bam", emit: bams
        path "alignment/*.bai", emit: bais
    
    script:
        """
        guppy_aligner -i $unaligned_bams -t ${task.cpus} --bam_out --index -a $reference_fasta -s alignment       
        
        """

}

process merge_bams {
    publishDir "${params.outdir}/aligned_bams", mode:'copy', pattern: "${params.sample_id}.ba*"
    container  "dhspence/docker-baseimage:latest"
    cpus 6
    memory '12 GB'
    
    input:
        path aligned_bams
        path aligned_bams_index
    output:
        path "${params.sample_id}.bam", emit: bam
        path "${params.sample_id}.bam.bai", emit: bai
    
    script:
        """
        samtools merge -@ ${task.cpus} -o ${params.sample_id}_unsorted.bam $aligned_bams &&
        samtools sort -o ${params.sample_id}.bam ${params.sample_id}_unsorted.bam &&
        samtools index ${params.sample_id}.bam
                
        """

}

process pepper {
    publishDir "${params.outdir}/pepper/haplotagged_bam", mode:'copy', pattern: "pepper_out/${params.sample_id}.haplotagged.bam"
    container  "kishwars/pepper_deepvariant:r0.8"
    cpus 1
    memory '12 GB'
    
    input:
        path aligned_merged_bam
        path aligned_merged_bam_index
        path reference_fasta
    output:
        path "pepper_out/${params.sample_id}.haplotagged.bam", emit: bam
    
    script:
        """
        run_pepper_margin_deepvariant call_variant \\
        -b $aligned_merged_bam \\
        -f $reference_fasta \\
        -o pepper_out \\
        -p ${params.sample_id} \\
        -t ${task.cpus} \\
        --${params.nanopore_reads_type} \\
        --phased_output \\
        --keep_intermediate_bam_files
                        
        """

}




// *********************************** //
// ******* Workflow Definitions ******* //
// *********************************** //

workflow {
    // Define inputs
    fast5_channel = Channel.fromPath(params.input) // Make input channels via "Channel Factory" constructors (https://www.nextflow.io/docs/latest/channel.html#frompath)
    reference_fasta_channel = Channel.fromPath(params.genome_fasta)

    // Workflow logic
    if( params.use_gpu ) {
        guppy_basecaller_gpu(fast5_channel)
        guppy_aligner(guppy_basecaller_gpu.out.basecall_bams_path, reference_fasta_channel)
    }
    else {
        guppy_basecaller(fast5_channel)
        guppy_aligner(guppy_basecaller.out.basecall_bams_path, reference_fasta_channel)
    }
    merge_bams (guppy_aligner.out.bams, guppy_aligner.out.bais)
    pepper (merge_bams.out.bam, merge_bams.out.bai, reference_fasta_channel)
}