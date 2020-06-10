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

# using Rscript to get peak model and cross correlation
# from the peak model, we could determin the shift size
```

### 2. ChIPQC(R)

```r
## Load libraries
library(ChIPQC)

## Load sample data
samples <- read.csv('meta/samplesheet_chr12.csv')

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="hg19") 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: Nanog and Pou5f1", reportFolder="ChIPQCreport")


```

### 2.1 merge replicates

```bash
# merge replicates with IDR
idr --samples Rep1_sorted_peaks.narrowPeak Rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file idr \
--plot \
--log-output-file idr.log

```

### 3. differential peaks

```r
library(DiffBind)
library(tidyverse)

# load data
samples <- read.csv('samplesheet_chr12.csv')
dbObj <- dba(sampleSheet=samples)

# pca plot
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)

#heatmap plot
plot(dbObj)

# establishing a contrast
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)

# differential analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

# different method of results
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS) 

# MA plots are a useful way to visualize the effect of normalization on data, as well as seeing which of the data points are being identified as differentially bound. 
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)

pvals <- dba.plotBox(dbObj)

# Extracting results
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

```

### 3 visulization deeptools

```bash
# bamCoverage
bamCoverage -b aln.bam \
-o visualization/bigWig/aln.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 

# if you have input bam file
# bamCompare against the input bam

# compute matrix
# reference point or scale
$ computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R ~/chipseq/results/visualization/refGenes.{bed|gtf} \
-S /n/groups/hbctraining/chip-seq/full-dataset/bigWig/Encode_Nanog*.bw \
--skipZeros \
-o ~/chipseq/results/visualization/matrixNanog_TSS_chr12.gz \
-p 6 \
--outFileSortedRegions ~/chipseq/results/visualization/regions_TSS_chr12.bed

# profile
plotProfile -m visualization/matrixNanog_TSS_chr12.gz \
-out visualization/figures/TSS_Nanog_profile.png \
--perGroup \
--colors green purple \
--plotTitle "" --samplesLabel "Rep1" "Rep2" \
--refPointLabel "TSS" \
-T "Nanog read density" \
-z ""

# heatmap
plotHeatmap -m visualization/matrixNanog_TSS_chr12.gz \
-out visualization/figures/TSS_Nanog_heatmap.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -4 --zMax 4  

# profile and heatmap
plotHeatmap -m visualization/matrixPou5f1_TSS_chr12.gz \
-out visualization/figures/TSS_Pou5f1_profile-heatmap.png \
--colorMap RdBu \
--zMin -2 --zMax 2  
```

### 4. functional analysis

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)


# Load data
samplefiles <- list.files("data/idr-bed", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("control", "Treat")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# annotatePeak
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

nanog_annot <- data.frame(peakAnnoList[["Nanog"]]@anno)

entrez <- nanog_annot$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                         keys = entrez,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)


# Run GO enrichment analysis 
ego <- enrichGO(gene = entrez, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

# Dotplot visualization
dotplot(ego, showCategory=50)

ekegg <- enrichKEGG(gene = entrez,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

dotplot(ekegg)

# Create a list with genes from each sample
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

# Run KEGG analysis
compKEGG <- compareCluster(geneCluster = genes, 
                         fun = "enrichKEGG",
                         organism = "human",
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

```

### 5. Motif discovery

[DREME website](http://meme-suite.org/tools/dreme) 

Tomtom: 
To determine if the identified motifs resemble the binding motifs of known transcription factors, we can submit the motifs to Tomtom, which searches a database of known motifs to find potential matches and provides a statistical measure of motif-motif similarity. We can run the analysis individually for each motif prediction by performing the following steps:

1. Click on the `Submit / Download` button for motif `ATGYWAAT` in the DREME output
2. A dialog box will appear asking you to Select what you want to do or Select a program. Select `Tomtom` and click `Submit`. This takes you to the input page. 
3. Tomtom allows you to select the database you wish to search against. Keep the default parameters selected, but keep in mind that there are other options when performing your own analyses.
4. Enter your email address and job description and start the search.

MEME-ChIP

MEME-ChIP is a tool that is part of the MEME Suite that is specifically designed for ChIP-seq analyses. MEME-ChIP performs DREME and Tomtom analysis in addition to using tools to assess which motifs are most centrally enriched (motifs should be centered in the peaks) and to combine related motifs into similarity clusters. It is able to identify longer motifs < 30bp, but takes much longer to run. The report generated by MEME-ChIP is available [here](https://github.com/hbctraining/Intro-to-ChIPseq/raw/master/chipseq_MEME%20ChIP_report.pdf).


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