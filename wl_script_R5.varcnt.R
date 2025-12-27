#!/usr/bin/env RScript

##### Parse parameters #####

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
        cat("\n")
        cat("Usage: Rscript --vanilla gds.r <Out_GDS> ncpus varcnt <VCF> [<VCF>...]\n")
        q(save="no") 
}

library(SeqArray,quiet=T)
#library(SeqVarTools)

##### Inputs #####
error = FALSE
ncpus <- as.numeric(args[2])
varcnt <- as.numeric(args[3]) 

vcffile <- paste0(args[4])
if (length(args) >4) {
        for (vcf in 5:length(args)) {
                vcffile <- c(vcffile,paste0(args[vcf]))
        }
}

for (vcf in vcffile) {
        if (!file.exists(vcf)) {
                cat (paste0("ERROR: ", vcf, " does NOT exist.\n"))
                error = TRUE
        }
}

if (error) {
        cat("\n")
        cat("ERROR(s) are detected. Terminate thr process.\n")
        q(save="no")
}


##### Output #####
gdsfile <- args[1]


# Load VCF
#seqVCF2GDS(vcffile, gdsfile, info.import="ABHet_One_subgroup", fmt.import=c("GT","AD","GQ"), storage.option="LZMA_RA", parallel=ncpus, verbose=TRUE)
seqVCF2GDS(vcffile, gdsfile, info.import="ABHet_One_subgroup", fmt.import=c("GT","AD","GQ"), storage.option = "LZ4_RA" , optimize = TRUE, parallel = ncpus, verbose = FALSE, variant_count = varcnt)
