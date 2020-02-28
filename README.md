# **(V)ariant (A)nalysis (P)ipeline**

Thank you for your interest in using the Variant Analysis Pipeline.
VAP is a comprehensive workflow for reference mapping and variant detection of genomic and transcriptomic reads using a suite of bioinformatics tools.

_**Article Source:**_ 

Variant analysis pipeline for accurate detection of genomic variants from transcriptome sequencing data
Adetunji MO, Lamont SJ, Abasht B, Schmidt CJ (2019) Variant analysis pipeline for accurate detection of genomic variants from transcriptome sequencing data. PLOS ONE 14(9): e0216838. https://doi.org/10.1371/journal.pone.0216838


## Bioinformatic tools 
Bioinformatic tools are grouped based on sequencing reads

### Genomic Sequencing 
* BOWTIE2
* BWA

### Transcriptomic Sequencing
* TOPHAT2
* STAR (2-PASS)
* HISAT2

### Variant Calling (for both Genomic/Transcriptomic Sequencing)
* PICARD + GATK HaplotypeCaller
  * sort, addreadgroups, markduplicates using **Picard Tools**.
  * split cigar reads using **GATK** from Transcriptomic Sequencing reads.
  * variant detection using **GATK**.

**N.B.** : parameters of all tools are set to default.

### Software used to design the VAP workflow are:

Software | Version
-------- | -------------
TopHat2 | 2.1.1
HiSAT2 | 2.1.0
STAR | 2.5.2b
SAMtools | 1.4.1
Picard tools | 2.13.2
GATK | 3.8
BWA-mem | 0.7.17
BOWTIE2 | 2.3.5.1

*Current pipeline is not compatible with GATK v4*

Contact maintainer to make custom changes to the different tools


## Things to be aware of 

### Job File
* change **_config_job.file_** file with settings or renamed as required.
* If parameters are not needed, they must be either removed or changed to *false*
  * Needed workflows (prefix: run) must be change to *true* (case-sensitive) 
	e.g : SAM = false (this means there is no sam file) else input the file directory ```/path/to/samfiles/*sam```
	e.g : runTopHAT = true (this means the pipeline should run TopHAT2)

### Indexes for Assembly tools and Variant Calling tools
Before running the pipeline. Create indexes for the different assemblers specified
**REFERENCE GENOME INDEX SYNTAXS:**
- GATK : 
```
	java -jar <picard directory>/picard.jar CreateSequenceDictionary R=<reference.fa> O=<reference.dict>
```
- HISAT :
```
	hisat2-build <reference.fa> <path to reference.fa>/<index_name>
```
- BOWTIE/TOPHAT :
```	
	bowtie2-build <reference.fa> <path to reference.fa>/<index_name>
```
- BWA :
```
	bwa index <reference.fa> 
```
_**N.B.** For easy use make sure all ```<index_name>``` should be the same and stored in the ```<reference.fa>``` directory_


### Downstream Merge and Filter Step (```runMergeFilter```)
The downstream step performs the following:
1. Merge SNPs from all variant calling tools initially specified to execute (_TopHAT2/HiSAT2/STAR_ or _BOWTIE/BWA_). 
1. Pre-set filtering criteria using GATK-VariantFiltration tool.
   1.   ReadRankPosSum (RRPS) < -8
   1.   Quality by depth (QD) < 5
   1.   Read depth (DP) < 10
   1.   Fisher’s exact test p-value (FS) > 60
   1.   Mapping Quality (MQ) < 40
   1.   SnpCluster (3 SNPs in 35bp)
   1.   Mann-Whitney Rank-Sum (MQRankSum) < -12.5
1. Exploratory statistics of all variant files.


## To run workflow
perl VariantAnalysisPipeline.pl -c config_job.file
