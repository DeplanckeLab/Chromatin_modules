# #################################################
# # Project: CM mapping
# # Script purpose: running the CRD pipeline for building CMs with Clomics
# # Date: 2023 Oct 12
# # Author: Vincent Gardeux (vincent.gardeux@epfl.ch), latest edits: Olga Pushkarev
# #################################################

### Options
options(echo=F, error=traceback)
args <- commandArgs(trailingOnly = TRUE)

install.packages("data.table")
# Packages
suppressPackageStartupMessages(require(data.table))

parse_mixed_symbols <- function(input_string){
    if (grepl(",", input_string)) {
        # Comma-separated values: Split by comma and keep them as strings
        elements <- unlist(strsplit(input_string, ","))
        chromosomes_to_iter <- elements
    }
    else if (grepl("-", input_string)){
        # Dash-separated values: Split by dash and convert to integers
        elements <- unlist(strsplit(input_string, "-"))
        chromosomes_to_iter <- as.integer(elements)[[1]]:as.integer(elements)[[2]]
    }
    else{
        chromosomes_to_iter <- c(input_string)
    }
    return(chromosomes_to_iter)
}

# Parameters
if(length(args) == 4){
    message("Running CRD (Clomics) R script:")
    message("1. Path to Clomics tool")
    message("2: Normalized matrices in bed.gz format from qtltools (separated by commas)")
    message("3. Chromosomes to use")
    message("4: Output folder")
}
else{
    message("ERROR: The script requires 4 parameters! Check the input parameters")
    quit(save="no")
}
param.path_to_clomics <- args[1]
param.matrix_files <- args[2]
param.chromosomes <- args[3]
param.output_folder <- args[4]

# In case of multiple files
param.matrix_files <- strsplit(param.matrix_files, split = ",")[[1]]
message(length(param.matrix_files), " matrix file(s) are found.")

# +
# Check if CLomics is installed
param.path_to_clomics <- gsub(
    x = param.path_to_clomics,
    pattern="\\\\",
    replacement="/"
)
if(!dir.exists(param.path_to_clomics)){
    message(
        "ERROR: Clomics cannot be found! Check if the path to the tool is correct or install Clomics (https://github.com/odelaneau/clomics)")
    quit(save="no")
}

# Get chromosomes
chromosomes <- parse_mixed_symbols(param.chromosomes)

# Checking output folder
param.output_folder <- gsub(
    x = param.output_folder,
    pattern="\\\\",
    replacement="/"
)

if(!endsWith(param.output_folder, "/")){
    param.output_folder <- paste0(param.output_folder, "/")
}
if(!dir.exists(param.output_folder)){
    dir.create(param.output_folder, recursive = T)
}
# -

generated.files <- ""
sep <- ""
for(param.matrix_file in param.matrix_files){
    # Checking sample name from vcmtools file name
    param.matrix_file <- trimws(param.matrix_file)
    data.filename <- strsplit(param.matrix_file, split = "/")[[1]]
    data.filename <- data.filename[length(data.filename)]
    data.filename <- strsplit(data.filename, "\\.")[[1]][1]
    data.filename <- strsplit(data.filename, "_")[[1]][1]
    message("Processing sample ", data.filename, "...")

    # Reading the normalized matrix
    data.matrix <- fread(param.matrix_file, sep="\t", header=T, data.table=F)
    message(ncol(data.matrix) - 6, " samples are found")
    message(nrow(data.matrix), " peaks are found")
    data.matrix$pid <- paste0(data.filename, ":", data.matrix$pid)
    data.matrix$gid <- data.matrix$pid

    # Prepare the BED pheno file for Clomics
    new.file <- paste0(param.output_folder, data.filename, ".clomics.bed")
    generated.files <- paste0(generated.files, sep, new.file, ".gz")
    sep <- " "
    data.clomics <- data.matrix
    colnames(data.clomics) <- c(
      "#chr", "start", "end", "id", "dummy", "strand",
      colnames(data.matrix)[7:(ncol(data.matrix))]
    )
    data.clomics <- data.clomics[rowSums(data.clomics[7:ncol(data.clomics)]) != 0, ] # Remove genes/peaks full of 0 because this makes Clomics to fail
    fwrite(
        data.clomics,
        file = new.file,
        sep = "\t",
        quote = F,
        col.names = T,
        row.names = F,
        scipen = 50
    )
    system(paste0("bgzip -f ", new.file)) # bgzip it
    system(paste0("tabix ", new.file, ".gz")) # tabix it (index)
}

# Running Clomics by chromosome
for(chr in chromosomes){
    cmd <- paste0(
        param.path_to_clomics,
        " build --region ",
        chr,
        " --bed ",
        generated.files,
        " --out ",
        param.output_folder,
        "clomics.tree.chr",
        chr,
        ".out"
    )
    message("Running command: ", cmd)
    system(command = cmd)

    cmd <- paste0(
        path_to_clomics,
        " call --tree ",
        param.output_folder,
        "clomics.tree.chr",
        chr,
        ".out --out ",
        param.output_folder,
        "clomics.tree.chr",
        chr,
        ".call.out"
    )
    message("Running command: ", cmd)
    system(command = cmd)

    cmd <- paste0(
        param.path_to_clomics,
        " locate --bed ",
        generated.files,
        " --tree ",
        param.output_folder,
        "clomics.tree.chr",
        chr,
        ".call.out --out ",
        param.output_folder,
        "clomics.tree.chr",
        chr,
        ".final.out"
    )
    message("Running command: ", cmd)
    system(command = cmd)
}

# Merging output of clomics by chr
data.clomics <- NULL
for(chr in chromosomes){
    message("Merging chr", chr)
    data.tmp <- fread(
      paste0(
          param.output_folder,
          "clomics.tree.chr",
          chr,
          ".final.out"
      ),
      sep = " ",
      header = T,
      data.table = F,
      fill = T
    )
    data.tmp$CRD <- data.tmp$MOD == 1 & data.tmp$ANN1 > 1
    data.tmp$chr <- chr
    data.tmp <- data.tmp[, c("chr", "START", "STOP", "UUID", "IDX", "CHILD1", "CHILD2", "CRD")]
    data.tmp$IDX = paste0("chr", chr, "_", data.tmp$IDX)
    data.tmp$CHILD1[!is.na(data.tmp$CHILD1)] = paste0(
        "chr",
        chr,
        "_",
        data.tmp$CHILD1
    )[!is.na(data.tmp$CHILD1)]
    data.tmp$CHILD2[!is.na(data.tmp$CHILD2)] = paste0(
        "chr",
        chr,
        "_",
        data.tmp$CHILD2
    )[!is.na(data.tmp$CHILD2)]
    data.tmp <- data.tmp[with(data.tmp, order(START, STOP, UUID)),]
    data.clomics <- rbind(data.clomics, data.tmp)
}  
fwrite(
    data.clomics,
    paste0(
        param.output_folder,\
        "CRD.all.clomics.bed"
    ),
    sep = "\t",
    col.names = F,
    row.names = F,
    na = "NA",
    quote = F,
    scipen = 50
)
