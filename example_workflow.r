# Example cohesion routine
require(SpiecEasi)
require(phyloseq)
require(Matrix)

source("helperFunctions.r")
source("calcCohesion_original.r")
source("calcCohesion.r")
source("calcCohesion_multi.r")
source("calcBrayCurtis.r")
ps_16S <- readRDS("data/NEON_ORNL_16S_ps.rds") # "ps" stands for PhyloSeq, a type of sequencing data object
ps_ITS <- readRDS("data/NEON_ORNL_ITS_ps.rds")

otu_16S = as(otu_table(ps_16S), "matrix") # script runs slow if you create null model on phyloseq object
otu_ITS = as(otu_table(ps_ITS), "matrix") # script runs slow if you create null model on phyloseq object

# calculate cohesion using the Herren and McMahon 2017 code/approach (using a null model)
cohesion_orig <- calcCohesionOrig(otu = otu_16S, pers.cutoff = .2, iter = 50)

# calculate cohesion using the Herren and McMahon 2017 approach (using a null model), just checking that my code produces the same result
cohesion <- calcCohesion(otu = otu_16S, pers.cutoff = .2, iter = 50)

# create a network via SpiecEasi, and use that as the cohesion custom-correlation input
lambda.min.ratio = .1
nlambda = 30
network_16S <- spiec.easi(ps_16S, method='glasso', lambda.min.ratio=lambda.min.ratio, sel.criterion='bstars',
                      nlambda=nlambda, pulsar.params = list(ncores=3))
evalSpiecEasi(network_16S)
cor_16S <- cov2cor(as.matrix(getOptCov(network_16S))) # create correlation matrix from output
cohesion_network_16S <- calcCohesion(otu = otu_16S, 
                                     custom.cor = cor_16S) # can't use persistence cutoff unless you apply it to spiec-easi input first

# create a multi-domain network via  SPIEC-EASI, get correlation matrix, use for multi-domain cohesion
lambda.min.ratio = .1
nlambda = 30
network_16S_ITS <- spiec.easi(list(ps_16S, ps_ITS), method='glasso', lambda.min.ratio = lambda.min.ratio, sel.criterion='bstars',
                     nlambda=30, pulsar.params = list(ncores=3))
evalSpiecEasi(network_16S_ITS)
cor_16S_ITS <- cov2cor(as.matrix(getOptCov(network_16S_ITS))) # create correlation matrix from output
cohesion_network_16S_ITS <- calcCohesionMulti(otu1 = otu_16S, 
                                     otu2 = otu_ITS, 
                                     custom.cor = cor_16S_ITS) # can't use persistence cutoff unless you apply it to spiec-easi input first

# compare results
# original cohesion script and mine (both using null models for correlations) - some randomness at low iterations
plot(cohesion_orig$NegativeCohesion ~ cohesion$NegativeCohesion)
plot(cohesion_orig$NegativeConnectedness ~ cohesion$NegativeConnectedness)
plot(cohesion_orig$PositiveCohesion ~ cohesion$PositiveCohesion)
plot(cohesion_orig$PositiveConnectedness ~ cohesion$PositiveConnectedness)

# check null-model results vs network correlation results
plot(cohesion_network_16S$NegativeCohesion ~ cohesion$NegativeCohesion)
plot(cohesion_network_16S$PositiveCohesion ~ cohesion$PositiveCohesion)

# check single-domain results vs multi-domain results
plot(cohesion_network_16S$NegativeCohesion ~ cohesion_network_16S_ITS$NegativeCohesion)
plot(cohesion_network_16S$PositiveCohesion ~ cohesion_network_16S_ITS$PositiveCohesion) 
# nonlinear strong relationship w/ positive cohesion, not much for negative cohesion

coh.merged <- cbind.data.frame(cohesion_orig$`PositiveCohesion`, 
                               cohesion_orig$`NegativeCohesion`,
                               cohesion_network_16S$`PositiveCohesion`,
                               cohesion_network_16S$`NegativeCohesion`,
                               cohesion_network_16S_ITS$`PositiveCohesion`,
                               cohesion_network_16S_ITS$`NegativeCohesion`)
colnames(coh.merged) <- c("posCoh.orig", "negCoh.orig","posCoh.network", "negCoh.network","posCoh.multi.network", "negCoh.multi.network")
coh.merged$sample1 <- rownames(coh.merged)


# calculate Bray-Curtis similarity - function is tailored for NEON data, not generic
bc_16S <- calcBrayCurtis(otu_16S, rarefy = 2000)
# merge with cohesion
df.merged <- merge(coh.merged, bc_16S)

# check how each can predict Bray-Curtis similarity
summary(lm(df.merged$bcSim ~ df.merged$negCoh.orig * df.merged$posCoh.orig))
summary(lm(df.merged$bcSim ~ df.merged$negCoh.network * df.merged$posCoh.network))
summary(lm(df.merged$bcSim ~ df.merged$negCoh.multi.network * df.merged$posCoh.multi.network))
# here, 16S-network-based outperforms original or multi-domain




