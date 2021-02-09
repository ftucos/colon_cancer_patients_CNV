##############################################################
#
# This script will automatically generate the bash script
# to process and aggregates all the samples
#
##############################################################

setwd("/Users/tucos/Documents/bioinformatic_projects/colon_cancer_patients_CNV")
library(tidyverse)

# Load files metadata
metadata <- read.delim("~/Documents/bioinformatic_projects/colon_cancer_patients_CNV/data/E-MTAB-8410.sdrf.txt")

# Function to write line bye line the script
writescript <- function(x){write(x, file="scripts/01-align_and_BAFextract", append=TRUE)}


# Build bash script ------------------------------------------------------
# add boilerplate code from template to the script owervriting any previous file
write(readLines("data/script_header.txt"), file="scripts/01-align_and_BAFextract", append=FALSE)

# Loop over file group
file_groups <- unique(metadata$Source.Name)
for(i in 1:length(file_groups)) {
  
  samples_to_aggregate <- metadata %>% filter(Source.Name == file_groups[[i]])
  writescript(paste0("\n\n# Processing sample: ", samples_to_aggregate[1,"Source.Name"], " ----------------------------------------------"))
    # download all the fastq files
  writescript("# download all the required fastq files\n")
  writescript("cd $FASTQ_DIR\n")
  for(j in 1:nrow(samples_to_aggregate)) {
    ## Download Read 1
    writescript(paste0("wget ", samples_to_aggregate[j, "Comment.FASTQ_URI."]))
    ## Download Read 2
    writescript(paste0("wget ", samples_to_aggregate[j, "Comment.FASTQ_URI..1"]))
  }
  
  # STAR solo align fastq files
  writescript("\n# STAR solo alignment\n")
  writescript(
    paste0(
      "$STAR --runThreadN 96 --genomeDir $GENOME_INDEX \\
      --readFilesIn ",
      paste(samples_to_aggregate[1:nrow(samples_to_aggregate), "Comment.read2.file."], collapse=","),
      " ",
      # Barcode + UMI reads
      paste(samples_to_aggregate[1:nrow(samples_to_aggregate), "Comment.read1.file."], collapse=","),
      "\\
      --readFilesCommand gunzip -c \\
      --outSAMtype BAM SortedByCoordinate \\
      --outFileNamePrefix $BAM_DIR/", samples_to_aggregate[1, "Source.Name"], ". \\
      --soloType CB_UMI_Simple --soloCBwhitelist $WD/data/737K-august-2016.txt \\
      --soloFeatures Gene"
      ))

  writescript("\ncd $BAM_DIR\n")

  # Move alignment logs 
  writescript("\n # Move Log Files \n")
  writescript("mv *Log* $LOG_DIR/")
  
  # Generate index
  writescript("\n# Generate index\n")
  writescript("samtools index *.bam")
  
  # BAFextract
  writescript("\n# BAF extract\n")
  writescript(paste0("mkdir -p $BAF_DIR/", samples_to_aggregate[1, "Source.Name"], "\n",
  "samtools view -h $BAM_DIR/", samples_to_aggregate[1, "Source.Name"], ".Aligned.sortedByCoord.out.bam | \\
  $BAFExtract -generate_compressed_pileup_per_SAM stdin $CHROM_SIZE $BAF_DIR/", samples_to_aggregate[1, "Source.Name"], "/ 10 0 \n
  $BAFExtract -get_SNVs_per_pileup $CHROM_SIZE $BAF_DIR/", samples_to_aggregate[1, "Source.Name"], "/ $GENOME_PILEUP 20 4 0.1 $BAF_DIR/", samples_to_aggregate[1, "Source.Name"], ".baf \n
  rm -rf $BAF_DIR/", samples_to_aggregate[1, "Source.Name"]))
  
  # Remove .fastq and .bam files
  writescript("rm -rf $FASTQ_DIR/* $BAM_DIR/*")
}

# make it executable
Sys.chmod("scripts/01-align_and_BAFextract", mode = "755")

