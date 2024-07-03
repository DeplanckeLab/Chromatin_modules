# #################################################
# Project: CM mapping across cell types
# Script purpose: running the Clomics pipeline for building CMs
# Date: 2022 Oct 23
# Author: Vincent Gardeux (vincent.gardeux@epfl.ch), latest edits: Olga Pushkarev
# #################################################

# Options
options(echo = F, error = traceback)

# Packages
suppressPackageStartupMessages(require(data.table))

parse_arguments <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    # Set default values for parameters
    path_to_clomics <- NULL
    path_to_src <- NULL

    bg_threshold <- "3"
    n_peaks <- "200"
    n_cores <- "1"
    chromosomes <- "1-22"

    dataset_name <- NULL
    matrix_files <- NULL
    output_folder <- NULL

    # Parse arguments
    if ("-path_to_clomics" %in% args) {
        arg1_index <- which(args == "-path_to_clomics")
        path_to_clomics <- args[arg1_index + 1]
    } else {
        message("ERROR: -path_to_clomics parameter is missing!")
        quit(save = "no")
    }

    if ("-path_to_src" %in% args) {
        arg2_index <- which(args == "-path_to_src")
        path_to_src <- args[arg2_index + 1]
    } else {
        message("ERROR: -path_to_src parameter is missing!")
        quit(save = "no")
    }

    if ("-bg_threshold" %in% args) {
        arg3_index <- which(args == "-bg_threshold")
        bg_threshold <- args[arg3_index + 1]
    }
    if ("-n_peaks" %in% args) {
        arg4_index <- which(args == "-n_peaks")
        n_peaks <- args[arg4_index + 1]
    }

    if ("-n_cores" %in% args) {
        arg5_index <- which(args == "-n_cores")
        n_cores <- args[arg5_index + 1]
    }

    if ("-chromosomes" %in% args) {
        arg6_index <- which(args == "-chromosomes")
        chromosomes <- args[arg6_index + 1]
    }

    if ("-dataset_name" %in% args) {
        arg7_index <- which(args == "-dataset_name")
        dataset_name <- args[arg7_index + 1]
    } else {
        message("ERROR: -dataset_name parameter is missing!")
        quit(save = "no")
    }

    if ("-matrix_files" %in% args) {
        arg8_index <- which(args == "-matrix_files")
        matrix_files <- args[arg8_index + 1]
    } else {
        message("ERROR: -matrix_files parameter is missing!")
        quit(save = "no")
    }

    if ("-output_folder" %in% args) {
        arg9_index <- which(args == "-output_folder")
        output_folder <- args[arg9_index + 1]
    } else {
        message("ERROR: -output_folder parameter is missing!")
        quit(save = "no")
    }

    message("Running Clomics:")
    message("1. Path to Clomics executable: ", path_to_clomics)
    message("2. Path to CRD to CM convertion python script: ", path_to_src)
    message("3. Background threshold: ", bg_threshold)
    message("4. Number of peaks to the right of the focal peak: ", n_peaks)
    message("5. Number of cores to use for duilding CMs: ", n_cores)
    message("6. Running Clomics on chromosomes: ", chromosomes)
    message("7. Dataset: ", dataset_name)
    message("8. Used count matrices: ", matrix_files)
    message("9. Output folder: ", output_folder)

    list(
        path_to_clomics = path_to_clomics,
        path_to_src = path_to_src,
        bg_threshold = bg_threshold,
        n_peaks = n_peaks,
        n_cores = n_cores,
        chromosomes = chromosomes,
        dataset_name = dataset_name,
        matrix_files = matrix_files,
        output_folder = output_folder
    )
}

# Parameters
args <- parse_arguments()
param_path_to_clomics <- args$path_to_clomics
param_path_to_src <- args$path_to_src
param_bg_threshold <- args$bg_threshold
param_n_peaks <- args$n_peaks
param_n_cores <- args$n_cores
param_chromosomes <- args$chromosomes
param_dataset_name <- args$dataset_name
param_matrix_files <- args$matrix_files
param_output_folder <- args$output_folder

# If multiple count matrices are found, separate them
param_matrix_files <- strsplit(param_matrix_files, split = ",")[[1]]
message(length(param_matrix_files), " matrix file(s) are found.")

# If multiple thresholds, create a list of thresholds
bg_thresholds <- strsplit(param_bg_threshold, split = ";")[[1]]
message(length(bg_thresholds), " thresholds are found.")

# Defining chromosomes
if (grepl("-", param_chromosomes)) {
    chromosomes_ab <- strsplit(param_chromosomes, split = "-")[[1]]
    chromosomes_list <- chromosomes_ab[1]:chromosomes_ab[2]
} else if (grepl(",", param_chromosomes)) {
    chromosomes_list <- strsplit(param_chromosomes, split = ",")[[1]]
} else {
    chromosomes_list <- param_chromosomes
}

# Checking output folder
param_path_to_clomics <- gsub(x = param_path_to_clomics, pattern = "\\\\", replacement = "/")
if (endsWith(param_path_to_clomics, "/")) {
    param_path_to_clomics <- substr(param_path_to_clomics, 1, nchar(param_path_to_clomics) - 1)
}

param_output_folder <- gsub(x = param_output_folder, pattern = "\\\\", replacement = "/")
if (endsWith(param_output_folder, "/")) {
    param_output_folder <- substr(param_output_folder, 1, nchar(param_output_folder) - 1)
}

param_output_folder <- paste0(param_output_folder, "/", "n_peaks_", param_n_peaks)

if (!dir.exists(param_output_folder)) {
    dir.create(
        paste0(
            param_output_folder, "/",
            param_dataset_name
        ),
        recursive = T
    )
}


generated.files <- ""
sep <- ""
for (param_matrix_file in param_matrix_files) {
    param_matrix_file <- trimws(param_matrix_file)
    data.filename <- strsplit(param_matrix_file, ":")[[1]][1]
    param_matrix_file <- strsplit(param_matrix_file, ":")[[1]][2]
    message("Processing ", data.filename, "...")

    new.file <- paste0(
        param_output_folder, "/",
        param_dataset_name, "/",
        data.filename, ".clomics.bed"
    )
    generated.files <- paste0(generated.files, sep, new.file, ".gz")
    if (!file.exists(new.file)) {
        # Reading the normalized matrix
        data.matrix <- fread(param_matrix_file, sep = "\t", header = T, data.table = F)
        message(ncol(data.matrix) - 6, " samples are found")
        message(nrow(data.matrix), " peaks are found")
        data.matrix$pid <- paste0(data.filename, ":", data.matrix$pid)
        data.matrix$gid <- data.matrix$pid

        # Prepare the BED pheno file for Clomics
        sep <- " "
        data.clomics <- data.matrix
        colnames(data.clomics) <- c("#chr", "start", "end", "id", "dummy", "strand", colnames(data.matrix)[7:(ncol(data.matrix))])

        # Remove genes/peaks full of 0 because this makes Clomics to fail:
        data.clomics <- data.clomics[rowSums(data.clomics[7:ncol(data.clomics)]) != 0, ]
        fwrite(data.clomics, file = new.file, sep = "\t", quote = F, col.names = T, row.names = F, scipen = 50)
        system(paste0("bgzip -f ", new.file)) # bgzip it
        system(paste0("tabix ", new.file, ".gz")) # tabix it (index)
    }
}

# Running Clomics by chromosome
for (chr in chromosomes_list) {
    build_cmd <- paste0(
        param_path_to_clomics, " build --region ", chr,
        " --bed ", generated.files,
        " --out ", param_output_folder, "/", param_dataset_name, "/clomics.tree.chr", chr, ".out",
        " --npeaks ", param_n_peaks
    )
    message("Running command: ", build_cmd)
    system(command = build_cmd)

    for (threshold in bg_thresholds) {
        if (!dir.exists(
            paste0(
                param_output_folder, "/",
                "bg_threshold_", threshold, "/",
                param_dataset_name
            )
        )
        ) {
            dir.create(
                paste0(
                    param_output_folder, "/",
                    "bg_threshold_", threshold, "/",
                    param_dataset_name
                ),
                recursive = T
            )
        }

        call_cmd <- paste0(
            param_path_to_clomics, " call ",
            "--tree ", param_output_folder, "/", param_dataset_name, "/clomics.tree.chr", chr, ".out ",
            "--threshold ", threshold, " ",
            "--out ", param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name, "/clomics.tree.chr", chr, ".call.out"
        )
        message("Running command: ", call_cmd)
        system(command = call_cmd)

        locate_cmd <- paste0(
            param_path_to_clomics, " locate --bed ", generated.files, " ",
            "--tree ", param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name, "/clomics.tree.chr", chr, ".call.out ",
            "--out ", param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name, "/clomics.tree.chr", chr, ".final.out"
        )
        message("Running command: ", locate_cmd)
        system(command = locate_cmd)
    }
}


# Merging output of Clomics by chr
for (threshold in bg_thresholds) {
    data.clomics <- NULL
    for (chr in chromosomes_list) {
        message("Merging chr", chr)
        data.tmp <- fread(
            paste0(
                param_output_folder, "/",
                "bg_threshold_", threshold, "/",
                param_dataset_name,
                "/clomics.tree.chr", chr, ".final.out"
            ),
            sep = " ",
            header = T,
            data.table = F,
            fill = T
        )
        data.tmp$CRD <- data.tmp$MOD == 1 & data.tmp$ANN1 > 1
        data.tmp$chr <- chr
        data.tmp <- data.tmp[, c("chr", "START", "STOP", "UUID", "IDX", "CHILD1", "CHILD2", "CRD")]
        data.tmp$IDX <- paste0("chr", chr, "_", data.tmp$IDX)
        data.tmp$CHILD1[!is.na(data.tmp$CHILD1)] <- paste0("chr", chr, "_", data.tmp$CHILD1)[!is.na(data.tmp$CHILD1)]
        data.tmp$CHILD2[!is.na(data.tmp$CHILD2)] <- paste0("chr", chr, "_", data.tmp$CHILD2)[!is.na(data.tmp$CHILD2)]
        data.tmp <- data.tmp[with(data.tmp, order(START, STOP, UUID)), ]
        data.clomics <- rbind(data.clomics, data.tmp)
    }
    fwrite(
        data.clomics,
        paste0(param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name, "/CRD.all.clomics.bed"),
        sep = "\t",
        col.names = F,
        row.names = F,
        na = "NA",
        quote = F,
        scipen = 50
    )
}

for (threshold in bg_thresholds) {
    message("Converting CRDs into CM format...")
    crd_to_cm_cmd <- paste0(
        "python3 ", param_path_to_src, "/0.crds_to_cms.py",
        " -d ", param_dataset_name,
        " -i ", param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name, "/CRD.all.clomics.bed ",
        " -o ", param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name,
        " -p ", param_n_cores
    )
    message("Running command: ", crd_to_cm_cmd)
    system(command = crd_to_cm_cmd)

    message("Merging batches...")
    merge_batches_cmd <- paste0(
        "python3 ", param_path_to_src, "/1.merge_batches.py",
        " -d ", param_dataset_name,
        " -o ", param_output_folder, "/", "bg_threshold_", threshold, "/", param_dataset_name
    )
    message("Running command: ", merge_batches_cmd)
    system(command = merge_batches_cmd)
}
