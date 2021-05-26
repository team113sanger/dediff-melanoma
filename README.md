# The Clinicopathologic Spectrum and Genomic Landscape of De- and Trans-differentiated Melanoma

This repository contains the code used for the sample processing and analysis for the paper "The clinicopathologic spectrum and genomic landscape of de-/trans-differentiated melanoma".

## Citation

If you use any data or code from this repository, please cite the paper.

## Funding

The work presented here was supported by the UKRI MRC grant MR/V000292/1 – "The Genomic Atlas of Dermatological Tumours (DERMATLAS)".

## QSubsec Scripts

Many of the scripts in this project use the [`qsubsec`](https://github.com/alastair-droop/qsubsec) template system. This allows the scripts to run either as bash scripts on a local machine, or to be submitted via a queueing engine on shared compute resources (such as [`LSF`](https://www.ibm.com/products/hpc-workload-management) or [`SGE`](https://arc.liv.ac.uk/trac/SGE)). Please see the [`documentation`](https://github.com/alastair-droop/qsubsec/tree/master/docs) for more information.

Most of the qsubsec scripts require a token file (`.tff`) to specify the variable used for in the script. Many of the parameters in these TFF files will need to be modified for your specific computation environment. Please see the documentation within the TFF files for guidance.

## GATK Image

The GATK processing scripts are performed using the GATK Docker image running via singularity 3.2.0.

The GATK image was obtained as:

~~~bash
singularity pull docker://broadinstitute/gatk:latest
singularity run gatk_latest.sif gatk --version
~~~

The version data ias as follows:

~~~plain
Using GATK jar /gatk/gatk-package-4.1.6.0-local.jar
Running:
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.1.6.0-local.jar --version
The Genome Analysis Toolkit (GATK) v4.1.6.0
HTSJDK Version: 2.21.2
Picard Version: 2.21.9
~~~

## Input Data

Raw sequence data are available from the [ENA](https://www.ebi.ac.uk/ena/browser/home):

* DNA: [EGAS00001003471](https://ega-archive.org/studies/EGAS00001003471)
* RNA: [EGAS00001003601](https://ega-archive.org/studies/EGAS00001003601)
