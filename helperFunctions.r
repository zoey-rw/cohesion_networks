# helper functions for cohesion scripts

####################create necessary functions######################

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

# evaluate SpiecEasi output
evalSpiecEasi <- function(se = se){
stability <- getStability(se) # stability should be close to .05
chosen.lambda <- getOptLambda(se) # the lambda the algorithm chose
lambda.path <- se$lambda # the lambda the algorithm was choosing from (the space between possible values is determined by lambda.min.ratio)
if (stability < .045) {
  cat(c(paste("Warning: stability is not close to .05 - revisit parameters:\nStability:", stability),
          paste("\nChosen lambda:", chosen.lambda),  
          ("\nLambda path:"), lambda.path),
      "\nSee explanations of these parameters by SpiecEasi developer: https://github.com/zdk123/SpiecEasi/issues/85")
}
}
