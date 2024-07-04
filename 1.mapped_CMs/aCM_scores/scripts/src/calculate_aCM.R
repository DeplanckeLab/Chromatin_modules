## Packages
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(ggplot2))

## QTLtools paths (sv-en)
qtltools_path <- "/software/qtltools-1.3/bin/QTLtools"
bgzip_path <- "/software/htslib-1.12/bgzip"
tabix_path <- "/software/htslib-1.12/tabix"

# Parameters
set.seed(1234)

## Options
options(echo=F, error=traceback)
args <- commandArgs(trailingOnly = TRUE)

## Get args
args <- commandArgs(trailingOnly = T)
if(length(args) != 9) stop("You should give NINE arguments to this R script")
param.crd_content <- args[1]
param.crd_position <- args[2]
param.peak_data <- args[3]
param.genotypes_vcf <- args[4]
param.output_merge_marks <- args[5]
param.output_proportion_variance_plot <- args[6]
param.output_aCRD_matrix <- args[7]
param.output_crdQTL <- args[8]
sample <- args[9]

# required input and output files
message("Running crdQTL R script:")
message("1: [INPUT] CRD content file: ",param.crd_content)
message("2: [INPUT] CRD track BED file: ",param.crd_position)
message("3: [INPUT] List of normalized count matrices used for generating the CRDs (clomics/qtltools *.bed.gz files, separated by commas):",param.peak_data)
message("4: [INPUT] genotypes VCF file: ",param.genotypes_vcf)
message("5: [OUTPUT] file name for merged count matrix with all marks: ", param.output_merge_marks)
message("6: [OUTPUT] proportion_variance_plot (pdf extension) output directory: ", param.output_proportion_variance_plot)
message("7: [OUTPUT] aCRD_matrix output directory:", param.output_aCRD_matrix)
message("8: [OUTPUT] crdQTLs output file: ", param.output_crdQTL)
message("9: [GENERAL] sample name is set to: ",sample)

if(length(args) != 9){
  message("Running normalization and regression R script:")
  message("ERROR: Requires 8 parameters!")
  quit(save="no")
}

# in advance, copy the input files for the specific sample to path/QTLtools_output/crdQTLs/ and rename to the required
# mark extension as present in the .vcm.content.txt file

# ##param.crd_content <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/CRDs_Clomics_regressed/LCL/vcm.content.txt")
# ##param.crd_position <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/CRDs_Clomics_regressed/LCL/vcm.tracks.bed")
# ##param.peak_data <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/normalized_matrices/PCs_regressed/LCL/H3K27ac_LCL_RPKMnorm_reg##ressed_qqnorm_forCRD.bed.gz",
#                     ##"~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/normalized_matrices/PCs_regressed/LCL/H3K4me1_LCL_RPKMnorm_regre##ssed_qqnorm_forCRD.bed.gz")
# ##param.genotypes_vcf <- c("~/SVRAW1/pushkare/PHM/genotypes/delaneau_genotypes/gencord_kgpgen_kgpseq.vcf.gz")
# ##param.output_merge_marks <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/aVCM_matrices_regressed/LCL/LCL_ALL_marks_matrix.txt")
# ##param.output_proportion_variance_plot <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/aVCM_matrices_regressed/LCL/LCL_aVCM.proportion.variance.PC1.p##df")
# ##param.output_aCRD_matrix <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/aVCM_matrices_regressed/LCL/LCL_aVCM.txt")
# ##param.output_crdQTL <- ##c("~/NAS1/vanmierl/CRDs_commonpeaksets_matrices/aVCM_matrices_regressed/LCL/LCL_crdQTLs.txt")
# ##sample <- c("LCL") # set sample name here 

## Functions
# For each CRD, generate a matrix with all peaks of the CRD, and perform PCA to compute aCRD
compute_aCRD = function(data.crd, data.all_marks){
  data.crd[,colnames(data.all_marks)] = NA # Prepare the aCRD matrix
  print(paste0("Computing the aCRD matrix for ", nrow(data.crd), " CRDs..."))
  prop.variance.pc1 = c() # Compute the proportion of variance explained by PC1
  progress = 1
  for(i in 1:nrow(data.crd)) # For each CRD
  {
    if(progress %% 500 == 0) message(paste0(progress, " CRDs processed (", round(100*(progress / nrow(data.crd)), 2), "%)"))
    peak.names = limma::strsplit2(data.crd[i,"Peaks"], ",")
    # Perform a PCA to summarize all peaks from this CRD
    data.pca = prcomp(na.omit(t(data.all_marks[peak.names,])))
    # Compute the percent variance explain by PC1
    prop.variance.pc1 = c(prop.variance.pc1, summary(data.pca)$importance[2,1])
    # Fill aCRD matrix
    data.aCRD = data.pca$x[,1]
    data.crd[i,names(data.aCRD)] = data.aCRD
    progress = progress + 1
  }
  message(paste0(progress-1, " CRDs processed (100%)"))
  data.crd$PropVariancePC1 = prop.variance.pc1
  data.crd
}

## Load CRD content from files
data.crd <- fread(param.crd_content, header=F, sep="\t", data.table=F)
colnames(data.crd) = c("CRD", "NbPeaks", "Peaks")

# Load Peak/Count matrix data
data.count_matrix <- list()
for(file in limma::strsplit2(param.peak_data, split = ",")){
    mark_file <- limma::strsplit2(file, ":")
    mark <- mark_file[[1]]  # Just get the name of the mark/TF
    file <- mark_file[[2]]
    data.count_matrix[[mark]] <- fread(file, header=T, sep="\t", data.table=F)
    data.count_matrix[[mark]]$pid <- paste0(mark, ":", data.count_matrix[[mark]]$pid)
}

## Create summary matrix
# Only focus on common patients across peaks
if(length(data.count_matrix) == 0) stop("There should be at least one peak file")
all_patients = colnames(data.count_matrix[[1]])[7:ncol(data.count_matrix[[1]])]

if(length(data.count_matrix) > 1) for(i in 2:length(data.count_matrix)) all_patients = intersect(all_patients, colnames(data.count_matrix[[i]])[7:ncol(data.count_matrix[[i]])])
all_patients = sort(unique(all_patients))
message(paste0(length(all_patients), " samples/patients overlap across ", length(data.count_matrix)," files"))

# Preparing the cols of the output matrix (all samples)
all_patients = c()
for(i in 1:length(data.count_matrix)) all_patients = c(all_patients, colnames(data.count_matrix[[i]])[7:ncol(data.count_matrix[[i]])])
all_patients = sort(unique(all_patients))
message(paste0(length(all_patients), " different samples/patients found in total"))

# Preparing the rows of the output matrix (all peaks)
all_peaks = c()
for(i in 1:length(data.count_matrix)) all_peaks = c(all_peaks, data.count_matrix[[i]][,"pid"])
n = length(all_peaks)
all_peaks = sort(unique(all_peaks))
message(paste0(n, " peaks found in total [", length(all_peaks), " different]"))
if(n != length(all_peaks)) stop("Peaks should have unique names")

# Creating the output matrix
data.all_marks = data.frame(row.names=all_peaks)
data.all_marks[,all_patients] = NA # Prepare the Peak matrix with NA for non existing values

# Fill the matrix
for(i in 1:length(data.count_matrix)) {
  data.all_marks[data.count_matrix[[i]][,"pid"],colnames(data.count_matrix[[i]])[7:ncol(data.count_matrix[[i]])]] = data.count_matrix[[i]][, 7:ncol(data.count_matrix[[i]])]
}
message("Writing merged count matrix of all peaks")
write.table(
    x = data.all_marks,
    file = param.output_merge_marks,
    sep="\t",
    row.names = T,
    col.names = T,
    na = "NA",
    quote = F)

# adapt rownames
rownames(data.all_marks) <- gsub("\\_.*gz","",rownames(data.all_marks))
if(!rownames(data.all_marks)[1] %like% "chr" & data.crd$Peaks[1] %like% "chr"){
  data.crd$Peaks <- gsub("chr","",data.crd$Peaks)
}

## Create aCRD matrix
message("Now generate aCRD scores...")
data.aCRD <- compute_aCRD(data.crd, data.all_marks)

message("Plotting distribution of PC1 proportion of explained variance (should be high, since peaks are correlated)")
p <- ggplot(
    data.aCRD,
    aes(x = PropVariancePC1 * 100)) + geom_density() + geom_vline(aes(xintercept=median(PropVariancePC1 * 100)),
                                                                  color="blue",
                                                                  linetype="dashed",
                                                                  size=1)
ggsave(
    p,
    filename = paste0(param.output_proportion_variance_plot, '/PC1_proportion_of_explained_variance.pdf'),
    width = 7,
    height = 7)

message("Writing aCRD matrix")
rownames(data.aCRD) = data.aCRD$CRD
data.aCRD <- data.aCRD[,colnames(data.all_marks)]
write.table(
    x = data.aCRD,
    file = paste0(param.output_aCRD_matrix, "/aCM_matrix.bed"),
    sep="\t",
    row.names = T,
    col.names = T,
    na = "NA",
    quote = F)

## Load CRD position files
crdspos <- fread(param.crd_position, data.table=F, sep="\t", header = F)
crdspos <- crdspos[, c("V4", "V1", "V2", "V3")]
colnames(crdspos) <- c("geneid", "chr", "start", "stop")
rownames(crdspos) <- crdspos$geneid
crdspos$chr <- gsub(x = crdspos$chr, pattern = "chr", replacement = "")

# Prepare the BED pheno file
data.crd.qtltools <- crdspos[, c("chr", "start", "stop", "geneid", "geneid")]
data.crd.qtltools$chr <- paste0("chr", data.crd.qtltools$chr)
colnames(data.crd.qtltools) <- c("#Chr", "start", "end", "pid", "gid")
data.crd.qtltools$strand <- "+"
data.crd.qtltools <- cbind(data.crd.qtltools, data.aCRD[rownames(data.crd.qtltools), ])
data.crd.qtltools2 <- data.crd.qtltools


# check if 'chr' should be added
message("Checking consistency of input file for mapping crdQTLs and .vcf file")

# Monocytes, Neutrophils, Tcells: no 'chr' in .vcf file
# iPSCs: 'chr' in .vcf file
# Grubert: 'chr' in .vcf file
# delaneau: 'chr' in .vcf file

if(sample == "Monocytes" | sample == "Neutrophils" | sample == "Tcells" | sample == "Gaffney") {
    message("In the related .vcf files there is no 'chr' so checking for presence here")
    if(data.crd.qtltools2$`#Chr`[1] %like% "chr"){
        message("chr found so removing it from the `#Chr` column")
        data.crd.qtltools2$`#Chr` <- gsub("chr","",data.crd.qtltools2$`#Chr`)
    }
   else {
       message("chr not present so no changes needed")
   }
} else if(sample == "Grubert" | sample == "LCL" | sample == "FIB" | sample == "iPSCs") {
    message("In the related .vcf files there is 'chr' so checking for presence here")
    if(data.crd.qtltools2$`#Chr`[1] %like% "chr") {
        message("chr found so no action required") 
    }
    else {
        message("adding 'chr' into the `#Chr` column")
        data.crd.qtltools2$`#Chr` <- paste0("chr", data.crd.qtltools2$`#Chr`)}
  }

# write matrices
message("Writing matrices to ", param.output_aCRD_matrix)
write.table(
    data.crd.qtltools2,
    file=paste0(param.output_aCRD_matrix, "/data_for_qtltools_on_CMs.bed"),
    sep="\t",
    quote=F,
    row.names=F,
    na = "NA"
)

if(file.exists(paste0(param.output_aCRD_matrix, "/data_for_qtltools_on_CMs.bed.gz"))){
    system(paste0("rm ", param.output_aCRD_matrix, "/data_for_qtltools_on_CMs.bed.gz"))
}

system(paste0(
    "bedtools sort -i ",
    param.output_aCRD_matrix,"/data_for_qtltools_on_CMs.bed -header > ",
    param.output_aCRD_matrix, "/sorted_data_for_qtltools_on_CMs.bed"))
system(paste0(bgzip_path, " ", param.output_aCRD_matrix, "/sorted_data_for_qtltools_on_CMs.bed")) # bgzip it
system(paste0(tabix_path, " -p bed ", param.output_aCRD_matrix, "/sorted_data_for_qtltools_on_CMs.bed.gz")) # tabix it (index)


# calculate crdQTLs
message("Computing crdQTLs using QTLtools")
system(
    paste0(
        qtltools_path,
        " cis --vcf ",
        param.genotypes_vcf,
        " --bed ",
        param.output_aCRD_matrix,
        "/sorted_data_for_qtltools_on_CMs.bed.gz --out ",
        param.output_crdQTL,
        " --permute 1000 --seed 42 --normal"
    )
)

