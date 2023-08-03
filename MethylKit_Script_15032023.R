library(methylKit)
########## Reading the methylation call files.

#### Define the list containing the bismark coverage files.
file.list1 <- list("SRR2102345_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz","SRR2102346_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz","SRR2102347_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz","SRR2102348_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz","SRR2102349_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz","SRR2102350_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz")

#### read the listed files into a methylRawList object making sure the other parameters are filled in correctly.
myobj <- methRead(file.list1, sample.id=list("vascular_1","vascular_2","vascular_3","nonvascular_1","nonvascular_2","nonvascular_3"), pipeline = "bismarkCoverage", assembly="Sb", treatment=c(1,1,1,0,0,0), mincov = 10)

#### check number of samples
myobj

#### What type of data is stored here?
head(myobj[[1]])

#### Get a histogram of the methylation percentage per sample Here for sample 1
pdf(file="historgram.pdf",width=12,height=12)
getMethylationStats(myobj[[1]], plot=TRUE, both.strands=FALSE)
dev.off()

#### Get a histogram of the read coverage per sample
pdf("coverage_vascular1.pdf",width=10,height=10)
getCoverageStats(myobj[[1]], plot=TRUE, both.strands=FALSE)
dev.off()

#### Get percentile data by setting plot=FALSE  (Note down the results or plot for each sample)
getCoverageStats(myobj[[1]], plot=FALSE, both.strands=FALSE)

#### Filtering samples based on read coverage.
myobj.filt=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

#### Again get percentile data by setting plot=FALSE after the filering samples (Note down the results or plot for each sample)
getCoverageStats(filtered.myobj[[1]], plot=FALSE, both.strands=FALSE)

#### Normalization
myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")


######## Comparative analysis.

#### Merging samples
meth <- unite(myobj.filt.norm, destrand=FALSE)

######## Further Filtering

#### get percent methylation matrix
pm=percMethylation(meth)

#### calculate standard deviation of CpGs
pdf("sds_plot.pdf",width=10,height=10)
sds=matrixStats::rowSds(pm)
dev.off()

######## Visualize the distribution of the per-CpG standard deviation

#### to determine a suitable cutoff
pdf("sds_after_fragment.pdf",width=10,height=10)
hist(sds, breaks = 100)
dev.off()

#### keep only CpG with standard deviations larger than 2%
meth <- meth[sds > 2]

#### This leaves us with this number of CpG sites
nrow(meth)

#### Sample Correlation
pdf("Correlation_plot.pdf",width=10,height=10)
getCorrelation(meth,plot=TRUE)
dev.off()

#### Clustering samples
pdf("Cluster_sample.pdf",width=10,height=10)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
pdf("screeplot.pdf",width=10,height=10)
PCASamples(meth, screeplot=TRUE)
dev.off()
PCASamples(meth)

######## Batch effects  (Not do for our analysis)
#sampleAnnotation=data.frame(batch_id=c("a","a","b","b"), age=c(1,2,3,4,5,6))
#as=assocComp(mBase=meth,sampleAnnotation)
# as

######## Tiling windows analysis (Not necessary for our analysis)
#myobj_lowCov = methRead(file.list,sample.id=list("Root_1","Root_2","Root_3","Shoot_1","Shoot_2","Shoot_3"),assembly="Sb",treatment=c(1,1,1,0,0,0),mincov = 3)

######## Finding differentially methylated bases or regions


#### Test for differential methylation... This might take a few minutes using multiple-cores.  (BH= Benjamini-Hochberg, MN=McCullagh and Nelder # (https://rdrr.io/github/al2na/methylKit/man/calculateDiffMeth-methods.html))
myDiff <- calculateDiffMeth(meth, overdispersion = "MN", adjust="BH", test="Chisq", mc.cores=20)    #myDiff=calculateDiffMeth(meth)
myDiff

#### Simple volcano plot to get an overview of differential methylation
plot(myDiff$meth.diff, -log10(myDiff$qvalue))
#abline(v=0) 

#### Overview of percentage hyper and hypo CpGs per chromosome (plot).
diffMethPerChr(myDiff)

#### Get hyper methylated bases and order by qvalue
pdf("hyper.pdf",width=10,height=10)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hyper <- myDiff25p.hyper[order(myDiff25p.hyper$qvalue),]
dev.off()

#### Get hypo methylated bases and order by qvalue
pdf("hypo.pdf",width=10,height=10)
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
myDiff25p.hypo <- myDiff25p.hypo[order(myDiff25p.hypo$qvalue),]
dev.off()

#### Get all differentially methylated bases and order by qvalue
myDiff25p <- getMethylDiff(myDiff,difference=25,qvalue=0.01)
myDiff25p <- myDiff25p[order(myDiff25p$qvalue),]

#### Plot the percentage hyper and hypo CpGs per chromosome.
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

######## Annotating differentially methylated bases or regions
library(genomation)

#### First load the annotation data; i.e the coordinates of promoters, TSS, intron and exons
refseq_anot <- readTranscriptFeatures("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.56.chr.bed.txt")

#### Annotate hypermethylated CpGs ("target") with promoter/exon/intron information ("feature"). This function operates on GRanges objects, so we first coerce the methylKit object to GRanges.
myDiff25p.hyper.anot <- annotateWithGeneParts(target = as(myDiff25p.hyper,"GRanges"),feature = refseq_anot)

#### Summary of hypre target set annotation
myDiff25p.hyper.anot
write.table(myDiff25p.hyper.anot, "myDiff25p.hyper.anot_txt", sep =" ")

#### Annotate hypomethylated CpGs ("target") with promoter/exon/intron information ("feature"). This function operates on GRanges objects, so we first coerce the methylKit object to GRanges.
myDiff25p.hypo.anot <- annotateWithGeneParts(target = as(myDiff25p.hypo,"GRanges"),feature = refseq_anot)

#### Summary of hypo target set annotation
myDiff25p.hypo.anot
write.table(myDiff25p.hypo.anot, "myDiff25p.hypo.anot_txt", sep =" ")

#### Total (hypo and hypre) annotate
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),refseq_anot)

#### target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))

#### To get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
diffmethregion = getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
write.table(diffmethregion, "diffmethregion_txt", sep=" ")


#### Plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
pdf("diffmethbase",width=10,height=10)
diffmethbase = plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")
dev.off()

#### convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")

#### Plot the CpG island annotation
plotTargetAnnotation(diffCpGann,col=c("blue","orange","pink"), main="differential methylation annotation")

#### To get percentage of intron/exon/promoters that overlap with differentially methylated bases.
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)


######## Reading the methylation call files and store them as flat file database

myobjDB=methRead(file.list1,sample.id=list("Root_1","Root_2","Root_3","Shoot_1","Shoot_2","Shoot_3"),assembly="Sb",treatment=c(1,1,1,0,0,0),context="CpG",dbtype = "tabix",dbdir = "methylDB")

#### Read the shores and flanking regions and name the flanks as shores and CpG islands as CpGi
cpg.obj=readFeatureFlank("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.56.chr.bed.txt",feature.flank.name=c("CpGi","shores"))

#### convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")

# Load the CpG info  (Same as above but not if you want to extract flanking region wise the we can use given below command (from line 131-138))
cpg_anot <- readFeatureFlank("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.56.chr.bed.txt", feature.flank.name = c("CpGi", "shores"), flank=2000)
diffCpGann2 <- annotateWithFeatureFlank(as(myDiff25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

# See whether the CpG in myDiff25p belong to a CpG Island or Shore
head(getMembers(diffCpGann2))
write.table(diffCpGann2, "diffCpGann2_txt", sep=" ")

# This can also be summarized for all differentially methylated CpGs
pdf("diffCpGann2",width=10,height=10)
plotTargetAnnotation(diffCpGann2, main = "Differential Methylation Annotation")
dev.off()

#### Regional analysis
promoters=regionCounts(myobj,refseq_anot$promoters)
head(promoters[[1]]) or promoters
write.table(promoters, "promoters_txt", sep=" ")

finalreport = getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
write.table(finalreport, "finalreport", sep=" ")

