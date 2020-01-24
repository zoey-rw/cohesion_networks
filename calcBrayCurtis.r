# calculate Bray-curtis and (optionally) spatial distance using NEON OTU table

# otu table should have samples as rows, and taxa as columns
# rownames should be sampleIDs
# 
# otu <- as.data.frame(ps_16s@otu_table)
# bc.df <- bray_curtis_calc(otu)
# bc.df <- bray_curtis_calc(otu, spatial.distance=T)
# bc.df <- bray_curtis_calc(otu, spatial.distance=T, rarefy=2000)

calcBrayCurtis <- function(otu = otu, 
                             rarefy = NULL, 
                             spatial.distance = FALSE, 
                             phys.path = "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/soil_phys_allsites.rds"){

  require(stringr)
  require(vegan)
  require(spatstat)
  require(fields)

# isolate sampling dates, put in order
sampleID <- rownames(otu)
siteID <- substr(sampleID, 1, 4)
YYYYMMDD <- str_extract(sampleID, "\\d\\d\\d\\d\\d\\d\\d\\d")
date <- as.Date(YYYYMMDD, "%Y%m%d")
metadata <- cbind.data.frame(siteID, sampleID, YYYYMMDD, date)
metadata <- metadata[order(metadata$date),]

# create matrix of date differences
days <- as.numeric(metadata$date)
i <- outer(days, days, '-') # get days since previous sample
i[upper.tri(i)] <- NA
date.df <- as.data.frame(i)
date.df$time2 <- metadata$date
colnames(i) <- metadata$date
rownames(i) <- metadata$date
melted <- reshape2::melt(i)
colnames(melted) <- c("time2","time1", "daysBetween")

# create bray-curtis similarity matrix
bray.dis <- vegan::vegdist(otu, method='bray', na.rm = T)
bray.sim <- 1 - as.matrix(bray.dis)
bray.sim[upper.tri(bray.sim)] <- NA
melted2 <- reshape2::melt(bray.sim)
colnames(melted2) <- c("sample2","sample1", "bcSim")


# combine BC with date matrix
out <- cbind(melted, melted2)
out <- out[!is.na(out$bcSim),]
out <- out[out$bcSim != 1,]
out$soilSample1 <- gsub("-GEN-DNA.", "", out$sample1)
out$soilSample2 <- gsub("-GEN-DNA.", "", out$sample2)


if(!is.null(rarefy)){
  cat(paste0("Rarefying to ", rarefy, " reads"))
  otu_rare <- otu[which(rowSums(otu) > rarefy),]
  otu_rare <- vegan::rrarefy(as.matrix(otu_rare), rarefy)
  bray.dis_rare <- vegan::vegdist(otu_rare, method='bray', na.rm = T)
  bray.sim_rare <- 1 - as.matrix(bray.dis_rare)
  bray.sim_rare[upper.tri(bray.sim_rare)] <- NA
  melted3 <- reshape2::melt(bray.sim_rare)
  colnames(melted3) <- c("sample2","sample1", "bcSim_rare")
  out <- merge(out, melted3)
}


# get literal distance (in km) between soil cores
if (spatial.distance) {
  if (is.null(phys.path)) {
    cat(
      "Please enter a filepath to NEON soil physical properties dataset with geolocations added."
    )
  } else {
    soil_phys_all <- readRDS(phys.path)
    sites <- sort(unique(metadata$siteID))
    out.list <- list()
    for (s in 1:length(sites)){
      site <- sites[[s]]
      soil_phys <- soil_phys_all[which(as.character(soil_phys_all$sampleID) %in% out$soilSample1),]
      dec.degrees.mat<-as.matrix(cbind(soil_phys$adjDecimalLongitude, soil_phys$adjDecimalLatitude))
      rownames(dec.degrees.mat)<-soil_phys$sampleID
      #create great circle matrix
      distance.matrix<-rdist.earth(dec.degrees.mat[,1:2],miles=FALSE)
      distance.mat <- reshape2::melt(distance.matrix)
      colnames(distance.mat) <- c("soilSample1","soilSample2","spatial_distance")
      distance.mat$soilSample1 <- as.character(distance.mat$soilSample1)
      distance.mat$soilSample2 <- as.character(distance.mat$soilSample2)
      distance.mat <- distance.mat[which(!is.na(distance.mat$spatial_distance)),]
      distance.mat <- distance.mat[!duplicated(distance.mat),]
      distance.mat <- distance.mat[distance.mat$spatial_distance != 0,]
      distance.mat$siteID <- sites[[s]]
      out.list[[s]] <- distance.mat
    }
    distance.mat.all <- do.call(rbind, out.list)
    out <- merge(out, distance.mat.all)
  }
} # end spatial distance section

return(out)

} # end function

