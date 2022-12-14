// *********************************** //
// ******* Process Definitions ******* //
// *********************************** //

// Write process directive (https://www.nextflow.io/docs/latest/dsl2.html#process)
process guppy_basecaller {
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/summary/*' // "**" means include subdirectories in glob pattern"
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/bams/*'
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/fastq/*'
    container  "dhspence/docker-guppy"
    cpus 16
    memory '48 GB'
    
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
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/fastq/*'
    container  "dhspence/docker-gguppy"
    cpus 2
    memory '16 GB'
    
    input:
        path fast5_path
    output:
        path "basecall/summary/sequencing_summary.txt", emit: summary
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

process pycoqc {
    container  "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    cpus 4
    memory '8 GB'
    
    input:
        path sequencing_summary

    output:
        path('*.json') , emit: pycoqc_json
    
    script:
        """
        pycoQC -f $sequencing_summary --json_outfile ${params.sample_id}.json
                
        """

}



process guppy_aligner {
    publishDir "${params.outdir}", mode:'copy', pattern: 'alignment/**', enabled: false // disabled in this example (but will be saved if enabled: true). These bam files are too big, and better to save the merged sorted bam from the next process (merge_bams)
    container  "dhspence/docker-guppy"
    cpus 16
    memory '72 GB'
    
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
    memory '24 GB'
    
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
    publishDir "${params.outdir}/pepper/haplotagged_bam", mode:'copy', pattern: "pepper_out/${params.sample_id}.haplotagged.bam", saveAs: { "${params.sample_id}.haplotagged.bam" }
    publishDir "${params.outdir}/pepper/vcf", mode:'copy', pattern: "pepper_out/${params.sample_id}*.vcf*", saveAs: { filename -> file(filename).getName() }
    container  "kishwars/pepper_deepvariant:r0.8"
    cpus 16
    memory '48 GB'
    
    input:
        path aligned_merged_bam
        path aligned_merged_bam_index
        path reference_fasta
    output:
        path "pepper_out/${params.sample_id}.haplotagged.bam", emit: bam
        path "pepper_out/${params.sample_id}*.vcf*", emit: vcf

    script:
        """
        export PATH=/opt/margin_dir/build/:$PATH
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

process index_pepper {
    publishDir "${params.outdir}/pepper/haplotagged_bam", mode:'copy', pattern: "*.bai"
    container  "dhspence/docker-baseimage:latest"
    cpus 12
    memory '72 GB'
    
    input:
        path pepper_bam

    output:
        path "*.bam.bai", emit: bai
    
    script:
        """
        samtools index $pepper_bam
                
        """

}



process mosdepth {
    publishDir "${params.outdir}/mosdepth", mode:'copy', pattern: "*", enabled: false
    container  "ghcr.io/dhslab/docker-mosdepth:latest"
    cpus 6
    memory '12 GB'
    
    input:
        path pepper_bam
        path pepper_bam_index
    output:
        path('*.global.dist.txt')      , emit: global_txt
        path('*.summary.txt')          , emit: summary_txt
        path('*.region.dist.txt')      , emit: regions_txt
        path('*.regions.bed.gz')       , emit: regions_bed
        path('*.regions.bed.gz.csi')   , emit: regions_csi
        path('*.quantized.bed.gz')     , emit: quantized_bed
        path('*.quantized.bed.gz.csi') , emit: quantized_csi
    
    script:
        """
        export PATH=/opt/conda/bin/:$PATH

        export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
        export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
        export MOSDEPTH_Q2=CALLABLE      # 5..149
        export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

        mosdepth -t ${task.cpus} -n -x -Q 1 --by 500 --quantize 0:1:5:150: ${params.sample_id} $pepper_bam
                
        """

}

process multiqc {
    publishDir "${params.outdir}/multiqc", mode:'copy', pattern: "multiqc_*"
    container  "quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0"
    cpus 4
    memory '12 GB'
    
    input:
        path multiqc_input_files

    output:
        path('multiqc_*')     , emit: multiqc_report

    
    script:
        """

        multiqc .
                
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
        basecall_summary_ch = guppy_basecaller_gpu.out.summary
    }
    else {
        guppy_basecaller(fast5_channel)
        guppy_aligner(guppy_basecaller.out.basecall_bams_path, reference_fasta_channel)
        basecall_summary_ch = guppy_basecaller.out.summary
    }

    pycoqc (basecall_summary_ch)
    merge_bams (guppy_aligner.out.bams, guppy_aligner.out.bais)
    pepper (merge_bams.out.bam, merge_bams.out.bai, reference_fasta_channel)
    index_pepper (pepper.out.bam)
    mosdepth (pepper.out.bam, index_pepper.out.bai)
    

    // MultiQc collection: the logic in the follwoing channel operation is to "dump" all files required for multiqc report in one channel (ch_multiqc_files). this channel will be input for multiqc process
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix( pycoqc.out.pycoqc_json)
    ch_multiqc_files = ch_multiqc_files.mix( mosdepth.out)
    ch_multiqc_files = ch_multiqc_files.collect()
    multiqc (ch_multiqc_files)
}
