# calculating cohesion using spiec-easi output


source("helperFunctions.r")
# 
# iter <- 30 # number of iterations to run for each taxon
# tax.shuffle <- T # whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
# use.custom.cors <- F # input your own correlation table
# #frequency <- .22 # fraction of samples that an OTU must be present in
# pers.cutoff <- 0.22
# 
# otu1 <- ps_harv_16S.filt@otu_table
# otu2 <- ps_harv_ITS.filt@otu_table

calcCohesionMulti <- function(otu1, otu2, pers.cutoff = 0.1, tax.shuffle = T, iter = 30, custom.cor = NULL){
  
  # Read in custom correlation matrix, if desired. Must set "use.custom.cors" to TRUE
  
    if(is.null(custom.cor)) {
      cat("No custom correlation path provided, so use.custom.cors must be set to FALSE.")
    } else {
      custom.cor.mat <- as.matrix(custom.cor)
      if((dim(otu1)[2] + dim(otu2)[2]) != dim(custom.cor.mat)[2]){
        cat("Dimensions of custom correlation matrix do not match those of abundance matrix.")
      }
    }
  
  # Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, but not if matrix c is count data

  # Create relative abundance matrices
  rel.1 <- otu1/rowSums(otu1)
  rel.2 <- otu2/rowSums(otu2)
  rel.d <- cbind(rel.1, rel.2)

  obs.exp.cors.mat <- custom.cor.mat
  
  diag(obs.exp.cors.mat) <- 0
  # Calculate connectedness by averaging positive and negative observed - expected correlations
  connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
  connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)
  # Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
  cohesion.pos <- rel.d %*% connectedness.pos
  cohesion.neg <- rel.d %*% connectedness.neg
  output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
  names(output) <- c("NegativeConnectedness", "PositiveConnectedness", "NegativeCohesion", "PositiveCohesion")
  return(output)
}




