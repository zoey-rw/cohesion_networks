# calculating cohesion, optionally using custom correlation table (i.e. spiec-easi output)

source("helperFunctions.r")

# iter <- 30 # number of iterations to run for each taxon
# tax.shuffle <- T # whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
# use.custom.cors <- F # input your own correlation table
# #frequency <- .22 # fraction of samples that an OTU must be present in
# pers.cutoff <- 0.22
# tax.shuffle = T

calcCohesion <- function(otu, pers.cutoff = 0.1, tax.shuffle = T, iter = 30, use.custom.cors = F, custom.cor = NULL){

if(!is.null(custom.cor)) {
  use.custom.cors <- T
}
    
# Read in custom correlation matrix, if desired. Must set "use.custom.cors" to TRUE
if(use.custom.cors == T) {
  if(is.null(custom.cor)) {
    cat("No custom correlation path provided, so use.custom.cors must be set to FALSE.")
  } else {
  custom.cor.mat <- as.matrix(custom.cor)
  if(dim(otu)[2] != dim(custom.cor.mat)[2]){
    cat("Dimensions of custom correlation matrix do not match those of abundance matrix.")
  } #Check that corr matrix and abun matrix have same dimensions
  }
}

c <- as.matrix(otu)
# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Remove taxa that are below the persistence cutoff, and any samples w/ no taxa
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1]) # number of zeroes allowed in a taxon's distribution
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
d <- d[rowSums(d) > 0, ]

# Create relative abundance and observed matrix.  
rel.d <- d / rowsums.orig
cor.mat.true <- cor(rel.d)
# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()

# Run this loop for the null model to get expected pairwise correlations
# Bypass null model if the option to input custom correlation matrix is TRUE
if(use.custom.cors == F) {
  if(tax.shuffle){
    for(which.taxon in 1:dim(rel.d)[2]){
      #create permuted correlations vector matrix to hold correlations from every permutation for each single otu
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
       # print(i)
        #Create empty matrix of same dimension as rel.d
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        #For each otu
        for(j in 1:dim(rel.d)[2]){ 
          # Replace the original taxon vector with a permuted taxon vector
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        # Do not randomize focal column 
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      if(which.taxon %% 10 == 0){print(paste(which.taxon,"/", dim(rel.d)[2], "taxa"))}
    }
  } else { # row shuffle approach
    for(which.taxon in 1:dim(rel.d)[2]){
      #create vector to hold correlations from every permutation for each single otu
      perm.cor.vec.mat <- vector()
      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d 
        #For each taxon
        for(j in 1:dim(rel.d)[1]){ 
          which.replace <- which(rel.d[j, ] > 0 ) 
          # if focal taxon is greater than zero, take it out of the replacement vector
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #Replace original taxon vector with a vector where the values greater than 0 have been randomly permuted 
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      if(which.taxon %% 20 == 0){print(paste(which.taxon,"/", dim(rel.d)[2], "taxa"))}
    }
  }
  obs.exp.cors.mat <- cor.mat.true - med.tax.cors # Save observed minus expected correlations
} else if(use.custom.cors == T){
  obs.exp.cors.mat <- custom.cor.mat
}


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




