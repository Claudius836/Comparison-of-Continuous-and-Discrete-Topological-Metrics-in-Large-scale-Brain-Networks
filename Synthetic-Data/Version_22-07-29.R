library(igraph)
library(mvtnorm)
set.seed(319)

#creating the number of simulations and dataframe
in_sim <- 10

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

n.sim <- 100 ### Number of (experimental) simulations. 
T <- 500  ### Number of time points. 
n <- 30   ### Number of subjects.
p <- 116  ### Number of ROIs, like in the real data set for the AAL parcellation scheme.
g <- 2
threshold <- 0.05
metric <- "Global_eff" ### flag or option (to be changed).
#metric <- "Local_eff" ### flag or option (to be changed).
#...

### Output object of dimension n.sim * 2 = 100*2.
data_frame <- data.frame("Weighted_p_value","Unweighted_p_value")
out_data

### Start sims.
for(s in 1:n.sim){
  ### Storage.
  #data <- list()
  ## Data frame with 4 colums. Exactly what you had before, with n rows = 30. 
  sim_data <- data.frame("Weighted_metric_grp1","Weighted_metric_grp2",
                           "Unweighted_metric_grp1","Unweighted_metric_grp2")
  for(j in 1:g){ ### Groups.
     for(i in 1:n){ ### Subjects.
        data[[j]][[i]] <- list()
        X <- rmvnorm(T,mean=rep(g-1,p),sigma=diag(p))
        X <- abs(X)
        C <- cov(X)
        
        ### Weighted graph
        C[C<threshold] <- 0
        weighted <- C
        weighted_graph <- graph_from_adjacency_matrix(C, weighted=T, mode="undirected", diag=F)
        data[[j]][[i]][[1]] <- weighted_graph
     
        ### Unweighted graph
        C[C<0.01] <- 0
        C[C>0.01] <- 1
        unweighted <- C
        unweighted_graph <- graph_from_adjacency_matrix(C, weighted=NULL, mode="undirected", diag=F)
  
        #Efficiencies for each group! 
        if(metric=="Global_eff") weighted_value <- global_efficiency(weighted_graph)
        if(metric=="Local_eff")  weighted_value <- average_local_efficiency(weighted_graph)

        ### Adapt this:
        data[[j]][[i]][[1]] <- weighted_value  ### Ideally, store it into a data_frame.
        data[[j]][[i]][[2]] <- unweighted_value  ### Ideally, store it into a data_frame.        
    }#n
  }#g
  ### Weighted TEST:
  data_frame with two columns for each group for the WEIGHTED metrics (e.g. global_eff), of dimension n*g = 30*2. 
  #t.test(rnorm(100),rnorm(100))$p.value
  #[1] 0.6971857
  Store p.value for that simulation, s. 
    
  ### Unweighted TEST:
  data_frame with two columns for each group for the UNWEIGHTED metrics (e.g. global_eff), of dimension n*g = 30*2. 
  #t.test(rnorm(100),rnorm(100))$p.value
  #[1] 0.6971857
  Store p.value for that simulation, s. 
    
  ### Final output of sims:
  #On each iteration, we create a new row for the variables for each simulation.
  data_frame <- rbind(data_frame, c(Weighted_p_value,Unweighted_p_value))   
}#n.sim

### Plot outputs
box plot of weighted vs. unweighted p-values.

### Redo everything choosing local eff. 
