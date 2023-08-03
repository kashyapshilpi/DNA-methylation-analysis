
#input data
inputmethfile <- Methfile_read(paths = "/home/scbb/GeneDMRs_result", suffix = ".gz")
inputrefseqfile <- Bedfile_read(paths = "/home/scbb/GeneDMRs_result", bedfile = "refseq", suffix = ".txt", feature = FALSE)

#If get all differentially methylated genes (DMGs) quickly
controls <- c("/home/scbb/shilpi_trainee/GeneDMRs_result1_1.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result1_2.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result1_3.gz")
cases <- c("/home/scbb/shilpi_trainee/GeneDMRs_result2_1.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result2_2.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result2_3.gz")
allDMGs <- Quick_GeneDMRs(paths = "/home/scbb/GeneDMRs_result", control_paths = controls, case_paths = cases)

#If get all differentially methylated cytosine sites (DMCs) quickly
controls <- c("/home/scbb/shilpi_trainee/GeneDMRs_result1_1.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result1_2.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result1_3.gz")
cases <- c("/home/scbb/shilpi_trainee/GeneDMRs_result2_1.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result2_2.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result2_3.gz")
allDMCs <- Quick_DMCs(control_paths = controls, case_paths = cases)


# read the methylation file #
inputmethfile <- Methfile_read(paths = paste(system.file(package = "GeneDMRs_result"), sep="", suffix = ".gz"))

# or if it is a case-control design #
controls <- c("/home/scbb/shilpi_trainee/GeneDMRs_result/1_1.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result/1_2.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result/1_3.gz")
cases <- c("/home/scbb/shilpi_trainee/GeneDMRs_result/2_1.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result/2_2.gz", "/home/scbb/shilpi_trainee/GeneDMRs_result/2_3.gz")
inputmethfile <- Methfile_read(control_paths = controls, case_paths = cases)

# quality control #
inputmethfile_QC <- Methfile_QC(inputmethfile)

# read the bedfile #
inputrefseqfile <- Bedfile_read(paths = "/home/scbb/GeneDMRs_result", bedfile = "refseq", suffix = ".txt", feature = FALSE)

 #methylation mean 
regiongeneall <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "all")

# sifnificant filter #
regiongeneall_significant <- Significant_filter(regiongeneall_Qvalue)

getCorrelation(meth,plot=FLASE)





