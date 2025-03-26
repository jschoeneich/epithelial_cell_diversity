#!/usr/bin/Rscript
library(Rsubread)
library(Rsamtools)

#Name of the analyis
NAME <-"bulk_full"

# where to save them
fc <- "/data/analysis_bulk"

if (!dir.exists(fc)){
	dir.create(fc,recursive = T)
}

# GTF File
gtf_file <- "~/genomes/mouse_genome_m38/Mus_musculus.GRCm38.102.gtf"

# Bam files
bam_path <- "/data/bams"

bamFiles <- list.files(bam_path, pattern = ".*bam$",full.names = T)

# works only for indexed bam-files
for (bf in bamFiles){
  print(testPairedEndBam(bf))
}

# setwd(bam_path)
# specifying the gtf file containing the annotations and the format of the file
# counting reads that match overlap exons and grouping exons by gene_id
count_matrix <- featureCounts(files = bamFiles,
                              GTF.featureType = "exon",
                              GTF.attrType = "gene_id",
                              annot.ext = gtf_file, 
                              isGTFAnnotationFile = TRUE,
                              isPairedEnd = TRUE,
                              nthreads = 10) 

print(head(count_matrix))

print(summary(count_matrix))

sample_ID <- colnames(count_matrix$counts) 
print(sample_ID)

write.csv(
  count_matrix$counts,
  file.path(fc,paste0("counts_",NAME,".csv"))
)

write.csv(
  count_matrix$annotation,
  file.path(fc,paste0("annotation_",NAME,".csv")),
  row.names=FALSE
)

write.csv(
  count_matrix$targets,
  file.path(fc,paste0("samples_",NAME,".csv")),
  row.names=FALSE
)

write.csv(
  count_matrix$stat,
  file.path(fc,paste0("stats_counts_",NAME,".csv")),
  row.names=FALSE
)

# Finally, as a good practice, always save the R image to your working directory
save.image(paste0(fc, "/featureCounts_workspace_",NAME, ".RData"))

