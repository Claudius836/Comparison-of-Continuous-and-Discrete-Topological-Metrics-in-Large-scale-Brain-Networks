library(igraph)
library(mvtnorm)
set.seed(319)

#creating the number of simulations and dataframe
in_sim <- 10
data_frame <- data.frame("Weighted_Global", "Weighted_Local", "Unweighted_Global", "Unweighted_Local")

visualise <- function(){
  #change arrow size and edge color:
  E(G)$arrow.size <- 1
  E(G)$edge.color <- "gray80"
  E(G)$width <-1+E(G)$weight/12
  V(G)$size <- V(G)$audience.size*0.7
  plot(G)
  
  
  ceb <- cluster_edge_betweenness(G) 
  dendPlot(ceb, mode="hclust")
  
  plot(ceb, G)
  
  length(ceb) # number of communities 
  membership(ceb)   # community membership for each node
  modularity(ceb)   # how modular the graph partitioning is
}


for (i in 1:100){
  n <- 100 # number of subjects not subject but time 
  p <- 30 # number of ROIs (Vertices)
  X <- rmvnorm(n,mean=rep(0,p),sigma=diag(p))
  X <- abs(X)
  C <- cov(X); dim(C) #creating the matrix
  C[C<0.05] <- 0 #thresholding
  #colnames(C) = rownames(C) = LETTERS[1:10] #assigning the col and row of the matrix
  #Creating the weighted grap
  weighted <- C
  weighted_graph <- graph_from_adjacency_matrix(C, weighted=T, mode="undirected", diag=F)
  
  #Creating the unweighted graph
  C[C<0.01] <- 0
  C[C>0.01] <- 1
  unweighted <- C
  unweighted_graph <- graph_from_adjacency_matrix(C, weighted=NULL, mode="undirected", diag=F)
  
  #visualise()
  
  #Efficiencies
  weighted_global_efficiency <- global_efficiency(weighted_graph)
  #weighted_elocal_efficiencyfficiency <- local_efficiencyfficiency(weighted_graph)
  weighted_average_local_efficiency <- average_local_efficiency(weighted_graph)
  
  unweighted_global_efficiency <- global_efficiency(unweighted_graph)
  #unweighted_local_efficiency <- local_efficiency(unweighted_graph)
  unweighted_average_local_efficiency <- average_local_efficiency(unweighted_graph)
  
  #On each iteration, we create a new row for the variables
  data_frame <- rbind(data_frame, c(weighted_global_efficiency, weighted_average_local_efficiency, unweighted_global_efficiency, unweighted_average_local_efficiency))
  
  
}

  
global_efficiency(weighted_graph)
average_local_efficiency(unweighted_graph)

mean(unweighted)

unweighted_local_efficiency
weighted_graph
unweighted_graph
X
C
