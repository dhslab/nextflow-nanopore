<br>

# Dev branch
# Nanopore WGS analysis pipeline  

nextflow-nanoproe is a Nextflow pipeline for analysis of Nanopore Whole Genome Sequencing.


## Pipeline summary
1. Basecalling and Alignment ([`Guppy`](https://nanoporetech.com/nanopore-sequencing-data-analysis))




## Usage
### Input files:
1. fast5 raw reads provided as a full path for the directory containing all fast5 files, either in a configuration file or as (--input path/to/fast5) command line parameter.
2. Path for reference genome fasta file, either in a configuration file or as (--genome path/to/genome.fasta) command line parameter.

## Running pipeline in LSF cluster (configured to WashU RIS cluster environment)

To run directly from GitHub, only a config file is needed (for testing: download and save test_ris.config file to launch directory). Use the following command:
```bash
NXF_HOME=${PWD}/.nextflow LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active $HOME:$HOME" bsub -g /dspencer/nextflow -G compute-dspencer -q dspencer -e nextflow_launcher.err -o nextflow_launcher.log -We 2:00 -n 2 -M 12GB -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -a "docker(mdivr/centos:v0.1)" nextflow run dhslab/nextflow-nanopore -r main -profile ris -c test_ris.config
```

Alternatively,  clone the repository and run the pipeline from local directory:
```bash
git clone https://github.com/dhslab/nextflow-nanopore.git
cd nextflow-nanopore/
LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active $HOME:$HOME" bsub < lsf_launch.sh
```
- Output:
  - "results/" is the desired output from the test run
  - "work/" is the working directory for all tasks, can be removed if the pipeline ran successfully
