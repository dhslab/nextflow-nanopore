<br>

# Nanopore WGS analysis pipeline

**nextflow-nanoproe** is a Nextflow pipeline for analysis of Nanopore Whole Genome Sequencing.


## Pipeline summary
1. Basecalling ([`Guppy`](https://nanoporetech.com/nanopore-sequencing-data-analysis)) - with GPU run option
1. Basecalling QC ([`PycoQC`](https://a-slide.github.io/pycoQC/))
1. Alignment ([`Guppy`](https://nanoporetech.com/nanopore-sequencing-data-analysis) with [`minimap2`](https://github.com/lh3/minimap2))
1. Merge all aligned bam files into asingle file ([`samtools`](http://www.htslib.org/doc/samtools.html))
1. Haplotyping and phased variants calling ([`PEPPER-Margin-DeepVariant`](https://github.com/kishwarshafin/pepper))
1. Depth calculation ([`mosdepth`](https://github.com/brentp/mosdepth))
1. MultiQC ([`MultiQC`](https://multiqc.info/)) for Basecalling (PycoQC) and Depth (mosdepth)





## Usage
### Input files:
1. fast5 raw reads provided as a full path for the directory containing all fast5 files, either in a configuration file or as (--input path/to/fast5) command line parameter.
2. Path for reference genome fasta file, either in a configuration file or as (--genome_fasta path/to/genome.fasta) command line parameter.

### Running pipeline in LSF cluster (configured to WashU RIS cluster environment)
#### Example: Test run in RIS cluster
- This test run takes input of:
  - few ".fast5" files (6 files: ~3.5 GB)
  - chr22.fasta as reference genome

#### 1) Running directly from GitHub:
```bash
LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active $HOME:$HOME" bsub -g /dspencer/nextflow -G compute-dspencer -q dspencer -e nextflow_launcher.err -o nextflow_launcher.log -We 2:00 -n 2 -M 12GB -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -a "docker(mdivr/centos:v0.1)" "NXF_HOME=${PWD}/.nextflow ; nextflow run dhslab/nextflow-nanopore -r main -profile ris,dhslab_test"
```

#### 2) Alternatively, clone the repository and run the pipeline from local directory:
```bash
git clone https://github.com/dhslab/nextflow-nanopore.git
cd nextflow-nanopore/
LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active $HOME:$HOME" bsub -g /dspencer/nextflow -G compute-dspencer -q dspencer -e nextflow_launcher.err -o nextflow_launcher.log -We 2:00 -n 2 -M 12GB -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -a "docker(mdivr/centos:v0.1)" "NXF_HOME=${PWD}/.nextflow ; nextflow run main.nf -profile ris,dhslab_test"
```

### Note:
If the pipeline is intended to be run from local code (after being cloned), instead of running: 
```
nextflow run main.nf -profile ris,dhslab_test
```
you can run:
```
nextflow run main.nf -profile ris -c conf/dhslab_test.config
```
The above two examples are interchangeable. As **dhslab profile** (defined in nextflow.config file) is basically just importing (or including in nextflow language) **conf/dhslab_test.config** file to the pipeline scope and append it to the configurations.

However "```-profile ris```" is still required in both cases as it is important to define the LSF runtime commands.
- Output:
  - "results/" is the desired output from the test run
  - "work/" is the working directory for all tasks, can be removed if the pipeline ran successfully

- Example for results output for sample "aml476081" in the test workflow
```
results/
├── aligned_bams
│   ├── aml476081.bam
│   └── aml476081.bam.bai
├── basecall
│   └── fastq
│       └── aml476081.fastq.gz
├── multiqc
│   ├── multiqc_data
│   │   ├── mosdepth_cov_dist.txt
│   │   ├── mosdepth_cumcov_dist.txt
│   │   ├── mosdepth_perchrom.txt
│   │   ├── multiqc.log
│   │   ├── multiqc_citations.txt
│   │   ├── multiqc_data.json
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc_sources.txt
│   │   └── pycoqc.txt
│   └── multiqc_report.html
├── pepper
│   ├── haplotagged_bam
│   │   ├── aml476081.haplotagged.bam
│   │   └── aml476081.haplotagged.bam.bai
│   └── vcf
│       ├── aml476081.phased.vcf.gz
│       ├── aml476081.phased.vcf.gz.tbi
│       ├── aml476081.vcf.gz
│       └── aml476081.vcf.gz.tbi
└── pipeline_info
    ├── pipeline_report.html
    └── pipeline_timeline.html
```
