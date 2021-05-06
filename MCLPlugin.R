
################################################################
## This scripts shows how run the MCL clustering algorithm
## developed by Stijn van Dongen with R
## source('~/Documents/enseignement/bioinformatics_courses/statistics_bioinformatics/R-files/network_topology/sylvain/mcl.R')
##
## Authors: Sylvain Brohee
## Date:    October 2008

########################################
## Functions
#oldlog <- log
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
#log <- oldlog

# Inflation step of MCL
inflate <- function (M,
		     inf) {
  M <- M^(inf);
  return (M);
}

# Normalize the matrix by column
norm <- function (M) {
  colum.sum <- apply(M,2,sum)
  M <- t(M) / colum.sum
  return (t(M))
}

# MCL procedure
mcl <- function (M, 	# Matrix
		 inf, 	# Inflation value
		 iter, 	# Number of iterations
		 verbose = F,
		 heatmaps = F
		 ) { 
print("FXN")
	print(nrow(M))
      	for (i in 1:iter) {
    old.M <- M;
    M.norm <- norm(M);
    #print("NORM")
    #print(nrow(M))
    #print(ncol(M))
    #print(nrow(M.norm))
    #print(ncol(M.norm))
    M <- M.norm%*%M.norm;
    #print(nrow(M))
    M <- inflate(M, inf);
    M <- norm(M);
    if (sum(old.M == M) == dim(M)[1]*dim(M)[2]) {
      break;
    }
    if (verbose) {
      #print (paste ("iteration", i));
      log(paste ("iteration", i));
    } 
  }
  return (M);
}


# collect.mcl.clusters
# having a MCL output, makes a data frame having 
# in the first column the node id
# in the second column the cluster to which this node belongs
# run through the columns of the mcl result matrix and looks for all nodes that compose a cluster
# Only selects clusters with more than the indicated number of elements
collect.mcl.clusters2 <- function (M, Mrows, size 	# Matrix (mcl result)
		                 ) {
  cluster.list <- list();
  cluster.number <- 1;
  M.names <- Mrows;#row.names(M);
  clustered.nodes <- vector(mode = "logical", length = dim(M)[1])
  for (i in 1:dim(M)[1]) {
    nodes <- M.names[which(M[i,] != 0)];
    if (length(nodes) >= size && !clustered.nodes[which(M[i,] != 0)]) {
		print (nodes);
		#Add the nodes to the cluster list
		cluster.list[[cluster.number]] <- nodes;
		
		clustered.nodes[which(M[i,] != 0)] = T;
		cluster.number <- cluster.number + 1;
    }
  }
  return(cluster.list);
}


input <- function(inputfile) {
#pc <- read.csv("Never.pvalued.csv");
pc <<- read.csv(inputfile, check.names=FALSE);
print("XXX")
print(nrow(pc))
print(ncol(pc))
}

run <- function() {
#pc <<- pc[,-1] #Status is present
#print(nrow(pc))
#print(ncol(pc))

mcl.data <<- as.matrix(pc)

#########################################################
#Second step, calculate the first clustering
#########################################################
#Do mapping
#Map (values in the config file)
mcl.data.ext <<- mcl.data;

# Set negative edges to zero (don't think they should be used in the clustering), new TMC
for (i in 1:nrow(mcl.data.ext)) {
   for (j in 1:ncol(mcl.data.ext)) {
      if (mcl.data.ext[i, j] < 0) {
         mcl.data.ext[i, j] <<- 0;
      }
   }
}

#Absolute Values
mcl.data.abs <<- abs(mcl.data.ext);

#Remove columns and rows with only 0s, normalization cannot be done in mcl if a column value is 0
mcl.data.abs <<- mcl.data.abs[,which(colSums(mcl.data.abs[,])!=0)]
mcl.data.abs <<- mcl.data.abs[which(rowSums(mcl.data.abs[,])!=0),]

# Launch MCL, LEVEL 1
inf_lev1 <<- 1.5
inf_lev2 <<- 1.5
minClusSize_L1 <<- 5
minClusSize_L2 <<- 5
mcl.clusters <<- mcl(mcl.data.abs,inf_lev1,2000, verbose = T, heatmaps=F);
#print("HERE")
#print(nrow(mcl.clusters))
#print(ncol(mcl.clusters))
junk <- matrix(mcl.clusters, ncol(mcl.data), ncol(mcl.data));
mcl.list <<- collect.mcl.clusters2(junk,colnames(mcl.clusters),minClusSize_L1);
}

output <- function(outputfile) {
mcl.listFull <- list();
for(i in 1:length(mcl.list)){
     if (i == 1) {
	write.table(mcl.list[[i]], file=paste(outputfile,"1.clusters.csv",sep="."), sep=",", append=FALSE, col.names=NA, na="");
	write.table(mcl.clusters, file=paste(outputfile,"1.clusters.values.csv",sep="."), sep=",", append=FALSE, col.names=NA, na="");
     }
     else {
	write.table(mcl.list[[i]], file=paste(outputfile,"1.clusters.csv",sep="."), sep=",", append=TRUE, col.names=NA, na="");
	write.table(mcl.clusters, file=paste(outputfile,"1.clusters.values.csv",sep="."), sep=",", append=TRUE, col.names=NA, na="");
     }
        if(length(mcl.list[[i]]) >= minClusSize_L1 && length(mcl.list[[i]]) > 1){
           new.data <- mcl.data[mcl.list[[i]], mcl.list[[i]]];
           new.data[which(new.data<0)] <- 0
           mcl.clusters2 <- mcl(new.data,inf_lev2,2000, verbose = F);
           mcl.list2 <- collect.mcl.clusters2(mcl.clusters2, colnames(mcl.clusters2),minClusSize_L2);
           mcl.listFull <- c(mcl.listFull, mcl.list2);
           if(length(mcl.list2) > 0){
                        for(j in 1:length(mcl.list2)){
     if (j == 1) {
	write.table(mcl.list2[[j]], file=paste(outputfile,"2.cluster", i, "csv", sep="."), sep=",", append=FALSE, col.names=NA, na="");
	write.table(mcl.clusters2, file=paste(outputfile,"2.cluster", i, "values.csv", sep="."), sep=",", append=FALSE, col.names=NA, na="");
     }
     else {
	write.table(mcl.list2[[j]], file=paste(outputfile,"2.cluster", i, "csv", sep="."), sep=",", append=TRUE, col.names=NA, na="");
	write.table(mcl.clusters2, file=paste(outputfile,"2.cluster", i, "values.csv", sep="."), sep=",", append=TRUE, col.names=NA, na="");
     }
                        }
           }
        }
}

write.table(mcl.data[rle(unlist(mcl.listFull))$values, rle(unlist(mcl.listFull))$values], file=paste(outputfile,"sortedcorrelations.csv", sep="."), sep=",", append=FALSE, col.names=NA, na="");

for (j in 1:length(mcl.listFull)) {
   if (j == 1) {
      write.table(mcl.listFull[[j]], file=paste(outputfile, "sortedclusters.csv", sep="."), sep=",", append=FALSE, col.names=NA, na="");
   }
   else {
      write.table(mcl.listFull[[j]], file=paste(outputfile, "sortedclusters.csv", sep="."), sep=",", append=TRUE, col.names=NA, na="");
   }
}

}




