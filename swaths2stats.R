#!/usr/bin/env Rscript

.libPaths("/opt/Rlibs/") 

suppressPackageStartupMessages(library(SWATH2stats))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("PECA"))

parser = ArgumentParser()

parser$add_argument("--input", 
                    help = "DIA-matrix from MPP.")

parser$add_argument("--design-file",
	dest="designfile",
	help="Design file for SWATH2stats")

args = parser$parse_args()

# Main Entry

if (is.null(args$input)) {
  print("--input is missing")
  stop()
}

# Read DIA data
data <- data.frame(fread(args$input, header=TRUE), stringsAsFactors=FALSE)
data$run_id <- basename(data$filename)
data <- reduce_OpenSWATH_output(data)
data = data[grep('iRT', data$ProteinName, invert=TRUE),]
data = data[grep('DECOY_', data$ProteinName, invert=TRUE),]

peptide_matrix = write_matrix_peptides(data,filename = "DIA-peptide-matrix")
write.table(peptide_matrix, sep="\t", file = "DIA-peptide-matrix.tsv", row.names=FALSE)

protein_matrix = write_matrix_proteins(data,filename = "DIA-protein-matrix")
write.table(protein_matrix, sep="\t", file = "DIA-protein-matrix.tsv", row.names=FALSE)

if (!is.null(args$designfile)) { 

   # Read design
   design <- read.table(args$designfile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
   data <- sample_annotation(data, design)

   # Convert for PECA
   data <- convert4PECA(data)

   # Run PECA for all possible combinations
   comb <- combn(unique(design$Condition),2)
   for(i in 1:ncol(comb)) {
   	 group1 <- paste(design$Condition,design$BioReplicate,sep="_")[design$Condition==comb[1,i]]
  	 group2 <- paste(design$Condition,design$BioReplicate,sep="_")[design$Condition==comb[2,i]]
  	 peca.out <- PECA_df(data, group1, group2, id="ProteinName", normalize="median", test="rots", progress=TRUE)
  	 write.table(peca.out, file=paste("PECA_",comb[1,i],"-",comb[2,i],".txt",sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
   }
}

