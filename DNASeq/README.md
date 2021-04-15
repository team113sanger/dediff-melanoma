# Project 5625 "Dedifferentiated Melanoma SNP Analysis"

This project has two main parts:

1. Automated "pipeline" step to generate SNP data via Mutect2;
2. Downstream (manual) analysis of the Mutect2 data.

NB: All the code snippets below assume that `$BASE` is the project base directory on lustre.

## Sequenza CNV Analysis

The CNV profiles of the samples are assessed using [sequenza](http://www.cbs.dtu.dk/biotools/sequenza/).

All the code and scripts for this analysis are in the `$BASE/sequenza` directory.

## Mutect2 Data Generation

### 1: Sample Metadata Collection

BAM statistics were collected from the canapps pipeline. These data are stored in `$BASE/metadata/2103-bamstats.json`.

### 2: Calculate BAM QC Data

The BAM idxstats and flagstats data were created using:

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-bamqc.qsubsec 5625.tff SAMPLE='FILE(../metadata/sample-ids.txt)' MUTECT2_TO_MAX_MEM='2048'
~~~

### 3: Estimate Contamination

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-estimate-contamination.qsubsec 5625.tff SAMPLE='FILE(../metadata/patients.txt)' STYPE='tumour,dediff'
~~~

### 4: Collect Artifacts

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-collect-artifacts.qsubsec 5625.tff SAMPLE='FILE(../metadata/patients.txt)' STYPE='tumour,dediff'
~~~

### 5: Run Variant Detection

**Notes**:

1. Due to an odd bug in Mutect2 for sample `PD42798c`, the single region `20 2673457 2673777` was removed from the original baitset file. The modified baitset file used for the analysis was `$BASE/reference/baitset/sureselect_5_9-padded-trimmed.bed`.
2. As patient `PD42799` does not have a tumour sample, it can not be run in the same submission command (hence the two commands below)

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-mutect2.qsubsec 5625.tff SAMPLE='PD42795,PD42796,PD42797,PD42798,PD42800,PD45781,PD45890' STYPE='tumour,dediff' MUTECT2_QUEUE='normal' MUTECT2_TO_MAX_MEM='8192'
qsubsec -ps 5625-mutect2.qsubsec 5625.tff SAMPLE='PD42799' STYPE='dediff' MUTECT2_QUEUE='normal' MUTECT2_TO_MAX_MEM='8192'
~~~

### 6: Merge Variant Files

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-merge.qsubsec 5625.tff MUTECT2_QUEUE='normal' MUTECT2_TO_MAX_MEM='8192'
~~~

### 7: Identify MNPs in all Samples

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-identify-mnps.qsubsec 5625.tff STATE="PRE,POST" MUTECT2_QUEUE='normal' MUTECT2_TO_MAX_MEM='8192'
~~~

### 8: Run VEP on merged data:

~~~bash
cd $BASE/scripts
qsubsec -ps 5625-vep.qsubsec 5625.tff
~~~

## Downstream Data Analysis

### 9: Run downstream Analyses

See the `$BASE/downstream` directory for the different downstream analyses run of the merged data.
