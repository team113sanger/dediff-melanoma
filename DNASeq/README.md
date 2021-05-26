# Project 5625 "Dedifferentiated Melanoma SNP Analysis"

This folder contains the scripts to process the dedifferentiated melanoma SNV analysis.

## Processing steps

**NB**:

* The bash variable `$SAMPLE` should contain the (comma-separated) list of samples to process;
* The base variable `$DTYPE` should contain the (comma-separated) list of disease types (`tumour`, or `dediff`)
* The bash variable `$BASE` should point to the DNASeq directory.
* These scripts use the [`qsubsec`](https://github.com/alastair-droop/qsubsec) template system. Make sure that the values specified in the token file `dediff-melanoma.tff` match your computational setup.

~~~bash
cd $BASE/scripts
qsubsec -ps dediff-mutect2.qsubsec dediff-melanoma.tff MEM_MAX=32768 CPU_MAX=20 SAMPLE=$SAMPLE DTYPE=$DTYPE
qsubsec -ps dediff-orientation-model.qsubsec dediff-melanoma.tff SAMPLE=$SAMPLE DTYPE=$DTYPE
qsubsec -ps dediff-pileup-summary.qsubsec dediff-melanoma.tff SAMPLE=$SAMPLE
qsubsec -ps dediff-calculate-contamination.qsubsec dediff-melanoma.tff SAMPLE=$SAMPLE DTYPE=$DTYPE
qsubsec -ps dediff-filter-mutect2.qsubsec dediff-melanoma.tff SAMPLE=$SAMPLE DTYPE=$DTYPE
qsubsec -ps dediff-annotate-mutect2.qsubsec dediff-melanoma.tff SAMPLE=$SAMPLE DTYPE=$DTYPE
~~~
