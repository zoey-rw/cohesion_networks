
# load NEON otu table
# otu <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/seqTab/CPER/CPER_otuTable.rds")
# # load NEON taxonomy and reformat
# tax <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/seqTab/CPER/CPER_taxTable.rds")
# for (i in 1:ncol(tax)) {  tax[, i] <- substring(tax[, i], 4) }
# tax <- as.matrix(tax)
# output.path <- "/projectnb/talbot-lab-data/zrwerbin/cohesion/CPER_16S_phyloseq.rds"
# 
# ps_cper_16S <- make_phyloseq(otu = otu, tax = tax)

make_phyloseq <- function(otu = otu, tax = tax, output.path = NULL){
  
  require(phyloseq)
  

# clean up taxonomy to match
tax <- tax[rownames(tax) %in% colnames(otu),]

# create metadata df
sample.names <- rownames(otu)
siteID <- substr(sample.names, 1, 4)
split <- strsplit(sample.names, c("\\-"))
dates <- lapply(split, "[[", 5)
month <- substr(dates, 1, 6)
metadata <- data.frame(sampleID = sample.names,  siteID = siteID, month = month)
metadata$horizon <- ifelse(grepl("\\-M\\-", sample.names), "M", "O")
rownames(metadata) <- metadata$sampleID
#metadata$names <- NA
#metadata[which(metadata$site=="HARV"),]$names <- "Harvard Forest, MA"

ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), #sample_names(metadata$sample),
                   sample_data(metadata), 
                   tax_table(tax))


# ps <- phyloseq(otu_table(t(otu), taxa_are_rows=TRUE), #sample_names(metadata$sample),
#                sample_data(metadata), 
#                tax_table(tax))

# remove long OTU names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)

# give the finest known taxonomic assignment
taxonomy <- tax_table(ps)
taxonomy[which(taxonomy[,7] != ""),7] <- paste(taxonomy[which(taxonomy[,7] != ""),][,6], taxonomy[which(taxonomy[,7] != ""),][,7], sep="_")
#tax_names <- t(apply(taxonomy, 1, function(x) x[max(which(x != "" & !is.na(x)))]))
tax_names <- t(apply(taxonomy, 1, function(x) x[max(which(x != ""))]))

#esv_num <- paste0("ESV", seq(ntaxa(ps)))
tax_names <- make.unique(as.character(tax_names), sep = "")
taxa_names(ps) <- tax_names
#esv_num(ps) <- esv_num # this doesn't work lol

if(!is.null(output.path)) saveRDS(tax, output.path, version=2)
else return(ps)
}
