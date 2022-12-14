// *********************************** //
// ******* Process Definitions ******* //
// *********************************** //

// Write process directive (https://www.nextflow.io/docs/latest/dsl2.html#process)
process guppy_basecaller {
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/**' // "**" means include subdirecteries in glob pattern"
    container  "dhspence/docker-guppy"
    cpus 12
    memory '48 GB'
    
    input:
        path fast5_path
    output:
        path "basecall/summary/*", emit: summary
        path "basecall/bams/*", emit: bams
        path "basecall/fastq/*", emit: fastqs
    
    script:
        """
        guppy_basecaller -i $fast5_path --bam_out -s unaligned_bam -c /opt/ont/guppy/data/dna_r10.4.1_e8.2_260bps_modbases_5mc_cg_sup.cfg --num_callers ${task.cpus}
        
        mkdir basecall
        mkdir basecall/summary && mv unaligned_bam/sequencing_* basecall/summary
        mkdir basecall/bams && mv unaligned_bam/pass/*.bam basecall/bams
        mkdir basecall/fastq && mv unaligned_bam/pass/*.fastq basecall/fastq

        """
}

process guppy_basecaller_gpu {
    publishDir "${params.outdir}", mode:'copy', pattern: 'basecall/**'
    container  "dhspence/docker-gguppy"
    cpus 2
    memory '32 GB'
    
    input:
        path fast5_path
    output:
        path "basecall/summary/*", emit: summary
        path "basecall/bams/*", emit: bams
        path "basecall/fastq/*", emit: fastqs
    
    script:
        """
        guppy_basecaller -i $fast5_path --bam_out -s unaligned_bam -c /opt/ont/guppy/data/dna_r10.4.1_e8.2_260bps_modbases_5mc_cg_sup.cfg --num_callers ${task.cpus} -x cuda:all:100%
        
        mkdir basecall
        mkdir basecall/summary && mv unaligned_bam/sequencing_* basecall/summary
        mkdir basecall/bams && mv unaligned_bam/pass/*.bam basecall/bams
        mkdir basecall/fastq && mv unaligned_bam/pass/*.fastq basecall/fastq
        """

}


// *********************************** //
// ******* Workflow Definitions ******* //
// *********************************** //

workflow {
    // Define inputs
    fast5_channel = Channel.fromPath(params.input) // Make input channels via "Channel Factory" constructors (https://www.nextflow.io/docs/latest/channel.html#frompath)

    // Workflow logic
    if( params.use_gpu )
        guppy_basecaller_gpu(fast5_channel)
    else
        guppy_basecaller(fast5_channel)
}