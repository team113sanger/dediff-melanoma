# 5627-dediff-rnaseq

Analysis scripts for the Dedifferentiated Melanoma RNAseq project

Alastair Droop, 2019-07-09

## Input Data

The input dataset consists of samples from 7 patients:

| patient ID   | Normal                  | Tumour     | Dedifferentiated Tumour |
|:-------------|:------------------------|:-----------|:------------------------|
| `PD42795`    | `PD42795b`              | `PD42795a` | `PD42795c`              |
| `PD42796`    | `PD42796b`              | `PD42796a` | `PD42796c`              |
| `PD42797`    | `PD42797b`              | `PD42797a` | `PD42797c`              |
| `PD42798`    | `PD42798b`              | `PD42798a` |                         |
| `PD42799`    | `PD42799b`              |            | `PD42799a`              |
| `PD42800`    | `PD42800b` & `PD42800c` | `PD42800a` | `PD42800d`              |
| `PD45781`    | `PD45781b` & `PD45781d` | `PD45781c` | `PD45781a`              |

**NB**:
* Patient `PD42798` has no dediff sample (as it failed QC)
* Patient `PD42799` has no tumour sample (no sample from patient)
* Patient `PD42800` has 2 normal samples
* Patient `PD45781` has 2 normal samples

The `$BASE/metadata/sample-map.json` file contains these data:

~~~json
{
	"PR42795": {
		"normal": ["PR42795b"],
		"tumour": ["PR42795a"],
		"dediff": ["PR42795c"]
	},
	"PR42796": {
		"normal": ["PR42796b"],
		"tumour": ["PR42796a"],
		"dediff": ["PR42796c"]
	},
	"PR42797": {
		"normal": ["PR42797b"],
		"tumour": ["PR42797a"],
		"dediff": ["PR42797c"]
	},
	"PR42798": {
		"normal": ["PR42798b"],
		"tumour": ["PR42798a"],
		"dediff": []
	},
	"PR42799": {
		"normal": ["PR42799b"],
		"tumour": [],
		"dediff": ["PR42799a"]
	},
	"PR42800": {
		"normal": ["PR42800b", "PR42800c"],
		"tumour": ["PR42800a"],
		"dediff": ["PR42800d"]
	},
	"PR45781": {
		"normal": ["PR45781b", "PR45781d"],
		"tumour": ["PR45781c"],
		"dediff": ["PR45781a"]
	}
}
~~~

### Sample List

The list of (valid) sample IDs is located at `$BASE/metadata/samples.txt`:

## Input BAM file statistics

Input BAM file statistics are retrieved from the core processing (canapps) metadata and stored in `$BASE/metadata/canapps-stats.json`.

## Processing Steps

Several analyses are performed on the data. For all of the below steps, the environment variable `$BASE` is assumed to be the project base directory on lustre.

### Pull down BAM files from canapps

It is far more efficient to work with local BAM files, so we download them to the project directory ($BASE/bam).


### Generate transcript counts with featureCounts

We use [`featureCounts`](http://subread.sourceforge.net) to count hits to transcripts. This populates the `counts` directory:

~~~bash
cd $BASE/scripts
qsubsec -ps 5627-featurecounts.qsubsec 5627.tff 
~~~

### Generate Tissue Profiles from Counts Data

~~~bash
cd $BASE/scripts
qsubsec -ps 5627-profiles.qsubsec 5627.tff
~~~

### Generate Fusion Profiles with STAR-Fusion

We use [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) to identify fusions. Before starting, make sure that the STAR-Fusion data is downloaded. The instructions for this are in `$BASE/reference/README.md`.

#### 1: Convert BAMs to FASTQ

STAR-Fusion requires input FASTQ, not BAM. To generate paired-end FASTQ, run:

~~~bash
cd $BASE/scripts
qsubsec -ps 5627-bam2fastq.qsubsec 5627.tff SAMPLE="FILE(../metadata/5627-samples.txt)"
~~~

#### 2: Run STAR-Fusion

~~~bash
cd $BASE/scripts
qsubsec -ps 5627-fusion.qsubsec 5627.tff SAMPLE="FILE(../metadata/5627-samples.txt)" MAX_MEM=65536 CPU_MIN=10
~~~

After running STAR-Fusion on all the samples, generate the collated fusion dataset with:

~~~bash
cd $BASE/scripts
# module load ISG/R # This is necessary for the internal canapps R version.
./5627-collate-fusions ../fusions ../metadata/5627-sample-map.json ../fusions/5627-fusions.csv
~~~

**NB**: This step is quick, so it has not been given its own qsubsec template.
