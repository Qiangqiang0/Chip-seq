# CHIP-seq


ref: https://github.com/macs3-project/MACS

ref: http://www.bio-info-trainee.com/2257.html

ref: https://www.sciencedirect.com/science/article/pii/S1046202320300591

ref: https://drompaplus.readthedocs.io/en/latest/

ref: https://github.com/hbctraining/Intro-to-ChIPseq
## 1. overview

process: QC --> BWA/bowtie -->  MACS --> bedtools merge --> deeptools(ploting)

chromosome annotation

motif analysis

prediction of gene expression level

chromatin loops


QC: 



## 2. Procudure

__pass QC, bwa/bowtie__ starting with MACS.

detail options: https://github.com/macs3-project/MACS

### 1. filter and map

```bash
# sambamba filter
# sambamba: http://www.360doc.com/content/19/0905/08/62751463_859209982.shtml
sambamba view -h  -f bam -F "[XS] == null and not unmapped and not duplicate" $align_sorted > $align_filtered

#peak calling
macs2 callpeak \
	-c control.bam \
	-t target.bam \
	-g hs \
	-f BAM \
	-n $target > ${name}.log


```

### 2. ChIPQC(R)

```bash


```

```bash
#bw ref: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#usage-example-for-chip-seq
bamCoverage -b $target.bam -o $target.bw

#
computeMatrix

```
output files:

1. name_peaks.xls

2. NAME_peaks.narrowPeak

3. NAME_summits.bed

4. NAME_peaks.broadPeak

5. NAME_peaks.gappedPeak

6. NAME_model.r

7. NAME_treat_pileup.bd



P.S.

Bed format: genome interval file, 

	for example, _chr, start_pos, end_pos_

bedGraph format: 1.must have track line; 2.the score is placed in col 4; 

	for example, _chr, start_pos, end_pos, score_

bigWig: generated from bedGraph