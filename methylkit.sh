
#call the library
library(methylKit)
library(IRanges)
library(GenomeInfoDb)
#input the data
file.list=list("SRR2102342_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz", "SRR2102343_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz", "SRR2102344_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz","SRR2102348_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz", "SRR2102349_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz", "SRR2102350_output_forward_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz")
#define the parameters
myobj <- methRead(file.list,
           sample.id=list("shoot_1","shoot_2","shoot_3","non_vascular_1","non_vascular_2","non_vascular_3"),
           pipeline = "bismarkCoverage",
           assembly="mm10",
           treatment=c(1,1,1,0,0,0),
           mincov = 10
           )
#call the data
myobj
#plot the  Descriptive statistics on samples
getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
# plots the histogram for percent methylation distribution.
pdf(file="shoot_vs_non_vascular_MethylationStats.Sample01.pdf",width=12,height=12) ;
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE); 
dev.off()
#plot the read coverage per base
pdf("shoot_vs_non_vascular_coverage_non_vascular2.pdf",width=10,height=10)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
dev.off()
#Filtering samples based on read coverage
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
  #Comparative analysis                                  
 #Merging samples                                    
 meth=unite(myobj, destrand=FALSE)
 head(meth)
 #Sample Correlation
 pdf("shoot_vs_non_vascular_correlation.pdf",width=10,height=10)
 getCorrelation(meth,plot=TRUE)
 dev.off()
 #Clustering samples
 pdf("shoot_vs_non_vascular_cluster_sample.pdf",width=10,height=10)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
# PCA analysis 
pdf("shoot_vs_non_vascular_PCA_sample.pdf",width=10,height=10)
PCASamples(meth, screeplot=TRUE)
dev.off()
#plot PC1 and PC2 axis and a scatter plot 
pdf("shoot_vs_non_vascular_PCA_sample2.pdf",width=10,height=10)
PCASamples(meth)
dev.off()
mat=percMethylation(meth)
mat[mat==100]=80
newobj=reconstruct(mat,meth)
#Tiling windows analysis
myobj_lowCov = methRead(file.list,
           sample.id=list("shoot_1","shoot_2","shoot_3","non_vascular_1","non_vascular_2","non_vascular_3"),
           pipeline = "bismarkCoverage",
           assembly="mm10",
           treatment=c(1,1,1,0,0,0),
           context="CpG",
           mincov = 3
           )
tiles = tileMethylCounts(myobj_lowCov,win.size=1000,step.size=1000,cov.bases = 10)
head(tiles[[1]],3)
#Finding differentially methylated bases or regions
myDiff=calculateDiffMeth(meth)
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
write.table(myDiff25p.hyper, "shoot_vs_non_vascular_myDiff25p.hyper.txt", sep =" ")
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
write.table(myDiff25p.hypo, "shoot_vs_non_vascular_myDiff25p.hypo.txt", sep =" ")
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
write.table(myDiff25p,"shoot_vs_non_vascular_myDiff25p.txt", sep ="\t")
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
#Correcting for overdispersion
sim.methylBase1<-dataSim(replicates=6,sites=1000,
                         treatment=c(rep(1,3),rep(0,3)),
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
                        )
my.diffMeth<-calculateDiffMeth(sim.methylBase1,
                                overdispersion="MN",test="Chisq",mc.cores=1)
 #Finding differentially methylated bases using multiple-cores
 myDiff=calculateDiffMeth(meth,mc.cores=2)
#Annotating differentially methylated bases or regions
library(genomation)
# read the gene BED file
gene.obj=readTranscriptFeatures("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.56.chr.bed.txt")
annotateWithGeneParts(as(myDiff25p,"GRanges"), gene.obj)
# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank("GCA_000003195.3_Sorghum_bicolor_NCBIv3_genomic.bed",feature.flank.name=c("CpGi","shores"))
# convert methylDiff object to GRanges and annotate
diffCpGann <- annotateWithFeatureFlank(as(myDiff25p,"GRanges"), feature = cpg.obj$CpGi, flank = cpg.obj$shores, feature.name = "CpGi", flank.name = "shores")
# Regional analysis
promoters=regionCounts(myobj,gene.obj$promoters)
head(promoters[[1]])
#Convenience functions for annotation objects
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
pdf("shoot_vs_non_vascular_diffAnn.pdf",width=10,height=10)
plotTargetAnnotation(diffAnn,precedence=TRUE,
main="differential methylation annotation")
dev.off()
pdf("shoot_vs_non_vascular_diffCpAnn.pdf",width=10,height=10)
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
main="differential methylation annotation")
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)
 

 

 
