# Project 5625 "Dedifferentiated Melanoma SNP Analysis" Sequenza Analysis

This folder contains the scripts used to run [sequenza](http://www.cbs.dtu.dk/biotools/sequenza/) on the project 5625 DNA samples.

**NB**: When running these scripts, make sure that you update the variables contained in the `5625-sequenza.tff` file.

## Prepare the Reference GC Track

~~~bash
cd $BASE/sequenza/scripts
qsubsec -ps 5625-sequenza-gcprep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff
~~~

## Pull down the BAMs

**NB**: This step is specific to Sanger.  Simply get your BAM files and place them in the directory specified by `{SEQUENZA_BAM_DIR}` in the TFF file.

~~~bash
cd $BASE/sequenza/scripts
qsubsec -ps 5625-sequenza-getbam.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE="FILE(../../metadata/sample-ids.txt)"
~~~

## Generate Pileups from BAMS

~~~bash
cd $BASE/sequenza/scripts
qsubsec -ps 5625-sequenza-pileup.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE="FILE(../../metadata/sample-ids.txt)"
~~~

## Prepare Sequenza .seqz Files

~~~bash
cd $BASE/sequenza/scripts
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42795 STYPE=tumour NORMAL=PD42795b TUMOUR=PD42795a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42795 STYPE=dediff NORMAL=PD42795b TUMOUR=PD42795c
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42796 STYPE=tumour NORMAL=PD42796b TUMOUR=PD42796a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42796 STYPE=dediff NORMAL=PD42796b TUMOUR=PD42796c
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42797 STYPE=tumour NORMAL=PD42797b TUMOUR=PD42797a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42797 STYPE=dediff NORMAL=PD42797b TUMOUR=PD42797c
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42798 STYPE=tumour NORMAL=PD42798b TUMOUR=PD42798a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42798 STYPE=dediff NORMAL=PD42798b TUMOUR=PD42798c
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42799 STYPE=dediff NORMAL=PD42799b TUMOUR=PD42799a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42800 STYPE=tumour NORMAL=PD42800b TUMOUR=PD42800a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42800 STYPE=dediff NORMAL=PD42800c TUMOUR=PD42800d
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD45890 STYPE=tumour NORMAL=PD45890b TUMOUR=PD45890a
qsubsec -ps 5625-sequenza-prep.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD45890 STYPE=dediff NORMAL=PD45890b TUMOUR=PD45890c
~~~

## Run Sequenza on .seqz Files

~~~bash
cd $BASE/sequenza/scripts
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42795 STYPE=tumour NORMAL=PD42795b TUMOUR=PD42795a GENDER=male SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42795 STYPE=dediff NORMAL=PD42795b TUMOUR=PD42795c GENDER=male SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42796 STYPE=tumour NORMAL=PD42796b TUMOUR=PD42796a GENDER=male SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42796 STYPE=dediff NORMAL=PD42796b TUMOUR=PD42796c GENDER=male SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42797 STYPE=tumour NORMAL=PD42797b TUMOUR=PD42797a GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42797 STYPE=dediff NORMAL=PD42797b TUMOUR=PD42797c GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42798 STYPE=tumour NORMAL=PD42798b TUMOUR=PD42798a GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42798 STYPE=dediff NORMAL=PD42798b TUMOUR=PD42798c GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42799 STYPE=dediff NORMAL=PD42799b TUMOUR=PD42799a GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42800 STYPE=tumour NORMAL=PD42800b TUMOUR=PD42800a GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD42800 STYPE=dediff NORMAL=PD42800c TUMOUR=PD42800d GENDER=female SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD45890 STYPE=tumour NORMAL=PD45890b TUMOUR=PD45890a GENDER=male SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
qsubsec -ps 5625-sequenza-analyse.qsubsec ../../scripts/5625.tff 5625-sequenza.tff SAMPLE=PD45890 STYPE=dediff NORMAL=PD45890b TUMOUR=PD45890c GENDER=male SQ_CPU_MIN=20 SQ_CPU_MAX=25 MAX_MEM=SQ_32768
~~~
