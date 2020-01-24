# Online script to generate cohesion metrics for a set of samples 
# CMH 26Apr17; cherren@wisc.edu
# adapted into function (minimal changes) by Zoey Werbin Jan 2020; zrwerbin@bu.edu

# User instructions: read in a sample table (in absolute or relative abundance) as object "b".
# If using a custom correlation matrix, read in that matrix at the designated line.
# Run the entire script, and the 4 vectors (2 of connectedness and 2 of cohesion) are generated for each sample at the end.
# Parameters that can be adjusted include pers.cutoff (persistence cutoff for retaining taxa in analysis), iter (number of iterations for the null model), tax.shuffle (whether to use taxon shuffle or row shuffle randomization), and use.custom.cors (whether to use a pre-determined correlation matrix)

source("helperFunctions.r")
###################################################################
###################################################################
### Workflow options ####
###################################################################
###################################################################


# for testing
# otu <- readRDS("CPER_otuTable.rds")
# tax.shuffle <- T # whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
# iter <- 30 # Larger values of iter mean the script takes longer to run (>= 200 is recommended)
# pers.cutoff <- 0.4 # Choose a persistence cutoff (min. fraction of taxon presence) for retaining taxa in the analysis

calcCohesionOrig <- function(otu, pers.cutoff = 0.1, tax.shuffle = T, iter = 30){

# Read in dataset
## Data should be in a matrix where each row is a sample. 
b <- otu


# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no empty samples or blank taxon columns. 
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# Optionally re-order dataset to be in chronological order. Change date format for your data. 
#c <- c[order(as.Date(rownames(c), format = "%m/%d/%Y")), ]

# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])

# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

# Create relative abundance matrix.  
rel.d <- d / rowsums.orig
# Optionally, check to see what proportion of the community is retained after cutting out taxa
#hist(rowSums(rel.d))

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()

# Run this loop for the null model to get expected pairwise correlations
if(tax.shuffle) {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
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
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 10 == 0){print(paste(which.taxon,"/", dim(rel.d)[2], "taxa"))}
    }
  } else {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d 
        
        #For each taxon
        for(j in 1:dim(rel.d)[1]){ 
          which.replace <- which(rel.d[j, ] > 0 ) 
          # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          
          #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 10 == 0){print(paste(which.taxon,"/", dim(rel.d)[2], "taxa"))}
    }
  }
  

# Save observed minus expected correlations. 
obs.exp.cors.mat <- cor.mat.true - med.tax.cors
diag(obs.exp.cors.mat) <- 0

#### 
#### Produce desired vectors of connectedness and cohesion 

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

####
#### Combine vectors into one list and print 
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("NegativeConnectedness", "PositiveConnectedness", "NegativeCohesion", "PositiveCohesion")

return(output)
}
