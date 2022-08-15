install.packages("igraph")
install.packages("mvtnorm")

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
time <- 500  ### Number of time points. 
n <- 30   ### Number of subjects.
p <- 116  ### Number of ROIs, like in the real data set for the AAL parcellation scheme.
g <- 2
m <- 2 #changing between local and global
threshold <- 0.05
h <- 0 #Start of each group. Group seperator
#metric <- "Global_eff" ### flag or option (to be changed).
#metric <- "Local_eff" ### flag or option (to be changed).
#...

### Output object of dimension n.sim * 2 = 100*2.
#data_frame <- data.frame("Weighted_p_value","Unweighted_p_value")
#out_data

#data <- data.frame("global_weighted", "global_unweighted", "local_weighted", "local_unweighted")
#p_value_data <- data.frame("global_w_p-value", "local_w_p-value", "global_u_p-value", "local_u_p-value")

p_value_data <- data.frame(1,2,3,4)
colnames(p_value_data) <- c("global_w_p-value", "local_w_p-value", "global_u_p-value", "local_u_p-value")
p_value_data <- p_value_data[-1,]

#group1 <- data.frame("global_weighted", "global_unweighted", "local_weighted", "local_unweighted")
#grpup2 <- data.frame("global_weighted", "global_unweighted", "local_weighted", "local_unweighted")
### Start sims.
for(s in 1:n.sim){
  ### Storage.
  #data <- list()
  ## Data frame with 4 colums. Exactly what you had before, with n rows = 30. 
  #sum_data <- data.frame("Weighted_metric_grp1","Weighted_metric_grp2",
  #                         "Unweighted_metric_grp1","Unweighted_metric_grp2")
  
  data <- data.frame(1,2,3,4)
  colnames(data) <- c("global_weighted", "global_unweighted", "local_weighted", "local_unweighted")
  data <- data[-1,]
  for(j in 1:g){ ### Groups.
     for(i in 1:n){ ### Subjects.
        #data[[j]][[i]] <- data_frame()
        X <- rmvnorm(time,mean=rep(g-1,p),sigma=diag(p))
        X <- abs(X)
        C <- cov(X)
        
       ### Weighted graph
        C[C<threshold] <- 0
        weighted <- C
        weighted_graph <- graph_from_adjacency_matrix(C, weighted=T, mode="undirected", diag=F)
        #data[[j]][[i]][[1]] <- weighted_graph
     
        ### Unweighted graph
        C[C<0.01] <- 0
        C[C>0.01] <- 1
        unweighted <- C
        unweighted_graph <- graph_from_adjacency_matrix(C, weighted=NULL, mode="undirected", diag=F)
        
        global_weighted_value <- global_efficiency(weighted_graph)
        global_unweighted_value <- global_efficiency(unweighted_graph)
  
        local_weighted_value <- average_local_efficiency(weighted_graph)
        local_unweighted_value <- average_local_efficiency(unweighted_graph)
        
        #for (i in 1:m){ #metric change
          #Efficiencies for each group! 
        #  if (m==1){
            #calculating metrics
        #    global_weighted_value <- global_efficiency(weighted_graph)
        #    global_unweighted_value <- global_efficiency(unweighted_graph)
        #  }
        #  if (m==2){
        #    #calculating metrics
        #    local_weighted_value <- average_local_efficiency(weighted_graph)
        #    local_unweighted_value <- average_local_efficiency(unweighted_graph)
        #  }
        #} #m
        
    #if(metric=="Global_eff") weighted_value <- global_efficiency(weighted_graph)
        #if(metric=="Local_eff")  weighted_value <- average_local_efficiency(weighted_graph)
        
    #if(metric=="Global_eff") unweighted_value <- global_efficiency(unweighted_graph)
    #if(metric=="Local_eff")  unweighted_value <- average_local_efficiency(unweighted_graph)
        
    ### Adapt this:
    
    #data[[j]][[i]][[j*1]] <- global_weighted_value ### Ideally, store it into a data_frame.
    #data[[j]][[i]][[j*2]] <- global_unweighted_value  ### Ideally, store it into a data_frame.
    #data[[j]][[i]][[j*3]] <- local_unweighted_value
    #data[[j]][[i]][[j*4]] <- local_unweighted_value
    
    #subject placement
      i <- h+i
        
      data[nrow(data) + 1,] = c(global_weighted_value,global_unweighted_value,local_weighted_value,local_unweighted_value)
    
    #data[i][1] <- global_weighted_value
    #data[i][2] <- global_unweighted_value
    #data[i][3] <- local_weighted_value
    #data[i][4] <- local_unweighted_value

    }#n
    #group number * length of subject to start the next append below the last subject. 
    #e.g number of subjects = 30. This starts the next subject(from group 2) on like 30.
    # h subject start. J is group number. n is subject count. 
    h <- j*n
  }#g
  
  #calculate metrics after iterating between both groups. 
  #seperate the dataframe depending on group position
  
  
  group1 <- as.data.frame(data[c(0:n),])
  #group1 <- group1[-1,]
  top_number <- 2*n
  group2 <- as.data.frame(data[c(31:top_number),])
  #group2 <- data(c(n+1:2*n),)
  
  #group1 <- data[[1:n],[:]]
  #group2 <- data[]
  
  #Weighted
  #global p-value
  p_value_data[nrow(p_value_data)+1,] <- c(t.test(group1["global_weighted"], group2["global_weighted"])$p.value,
                                           t.test(group1["local_weighted"], group2["local_weighted"])$p.value,
                                           t.test(group1["global_unweighted"], group2["global_unweighted"])$p.value,
                                           t.test(group1["local_unweighted"], group2["local_unweighted"])$p.value)
  #p_value_data[s][1] <- t.test(group1["global_weighted"], group2["global_weighted"])$p.value
  #local
  #p_value_data[s][2] <- t.test(group1, group2["global_u_p-value"]$p.value
  
  #unweighted
  #p_value_data[s][3] <- t.test(group1(,2), group2(,2))$p.value
  
  #p_value_data[s][4] <- t.test(group1(,4), group2(,4))$p.value
  
  rm(data)
  
  ### Weighted TEST:
  #data_frame with two columns for each group for the WEIGHTED metrics (e.g. global_eff), of dimension n*g = 30*2. 
  #t.test(rnorm(100),rnorm(100))$p.value
  #[1] 0.6971857
  #Store p.value for that simulation, s. 
    
  ### Unweighted TEST:
  #data_frame with two columns for each group for the UNWEIGHTED metrics (e.g. global_eff), of dimension n*g = 30*2. 
  #t.test(rnorm(100),rnorm(100))$p.value
  #[1] 0.6971857
  #Store p.value for that simulation, s. 
    
  ### Final output of sims:
  #On each iteration, we create a new row for the variables for each simulation.
  #data_frame <- rbind(data_frame, c(Weighted_p_value,Unweighted_p_value))   
}#n.sim



boxplot(p_value_data)


### Plot outputs
#box plot of weighted vs. unweighted p-values.

### Redo everything choosing local eff. 
