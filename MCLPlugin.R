
################################################################
## This scripts shows how run the MCL clustering algorithm
## developed by Stijn van Dongen with R
## source('~/Documents/enseignement/bioinformatics_courses/statistics_bioinformatics/R-files/network_topology/sylvain/mcl.R')
##
## Authors: Sylvain Brohee
## Date:    October 2008

########################################
## Functions

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
  for (i in 1:iter) {
    old.M <- M;
    M.norm <- norm(M);
    print(dim(M.norm));
    M <- M.norm%*%M.norm;
    M <- inflate(M, inf);
    M <- norm(M);
    if (sum(old.M == M) == dim(M)[1]*dim(M)[2]) {
      break;
    }
    if (verbose) {
      print (paste ("iteration", i));
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
		#print (nodes);
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
pc <<- read.csv(inputfile);
}

run <- function() {
pc <<- pc[,-1] #Status is present

mcl.data <- as.matrix(pc)

#########################################################
#Second step, calculate the first clustering
#########################################################
#Do mapping
#Map (values in the config file)
mcl.data.ext <- mcl.data;

# Set negative edges to zero (don't think they should be used in the clustering), new TMC
for (i in 1:nrow(mcl.data.ext)) {
   for (j in 1:ncol(mcl.data.ext)) {
      if (mcl.data.ext[i, j] < 0) {
         mcl.data.ext[i, j] <- 0;
      }
   }
}

#Absolute Values
mcl.data.abs <- abs(mcl.data.ext);

#Remove columns and rows with only 0s, normalization cannot be done in mcl if a column value is 0
mcl.data.abs <- mcl.data.abs[,which(colSums(mcl.data.abs[,])!=0)]
mcl.data.abs <- mcl.data.abs[which(rowSums(mcl.data.abs[,])!=0),]

# Launch MCL
inf_lev1 <- 1.5
minClusSize_L1 <- 3
mcl.clusters <- mcl(mcl.data.abs,inf_lev1,2000, verbose = T, heatmaps=F);
junk <- matrix(mcl.clusters, ncol(mcl.data), ncol(mcl.data));
mcl.list <<- collect.mcl.clusters2(junk,colnames(mcl.clusters),minClusSize_L1);
}


output <- function(outputfile) {
for(j in 1:length(mcl.list)){
	write.table(mcl.list[[j]], file=outputfile, sep=",", append=TRUE, col.names=NA, na="");
}
}




