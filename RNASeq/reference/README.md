# RNASeq Reference Files

This folder contains the reference data necessary for the RNASeq analysis.

## 1000genomes Reference

We need to contact canapps to get details of exactly how their reference was generated.

## STAR-Fusion Reference

The `star-fusion` directory contains the singularity image and reference data for running [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) on Ingrid's RNAseq data.

To run the STAR-Fusion singularity image, we need two things:

1. The singularity image; and
2. The reference dataset

### 1: Obtain the Singularity Image

~~~bash
cd $BASE/reference/star-fusion
curl --output star-fusion.v1.8.1.simg https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/star-fusion.v1.8.1.simg
~~~

### 2: Obtain the STAR-Fusion Reference Dataset

~~~bash
curl --output ${BASE}/reference/star-fusion/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz
tar -zxvf ${BASE}/reference/star-fusion/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz --directory ${BASE}/reference/star-fusion
mv ${BASE}/reference/star-fusion/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir ${BASE}/reference/star-fusion/sflib-GRCh37-19
rmdir ${BASE}/reference/star-fusion/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play
rm ${BASE}/reference/star-fusion/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz
~~~
