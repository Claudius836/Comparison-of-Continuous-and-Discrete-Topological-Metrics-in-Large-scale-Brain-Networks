library(igraph)
set.seed(319)

matrix_demo <- dget("matrix_demo.R")

demographics <- as.data.frame(matrix_demo[2])

demographics["demographics.sex"]

demographics[4]
demo_sex <- demographics[4]
summary(demo_age)

male <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("global_weighted", "global_unweighted", "local_weighted", "local_unweighted"))))
female <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("global_weighted", "global_unweighted", "local_weighted", "local_unweighted"))))
p_value_data <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("global_w_p-value", "local_w_p-value", "global_u_p-value", "local_u_p-value"))))


for (i in 1:nrow(demo_sex)){
  gender <- demo_sex[i, 'demographics.sex']
  c <- matrix_demo[[c(1, i)]] #current matrix
  c <- abs(c)
  
  #create weighted graph 
  weighted_graph <- graph_from_adjacency_matrix(c, weighted=T, mode="undirected", diag=F)
  
  #turn the data into binary
  c[c<0.5] <- 0
  c[c>0.5] <- 1
  #create unweighted graph. 
  unweighted_graph <- graph_from_adjacency_matrix(c, weighted=NULL, mode="undirected", diag=F)
  
  #creating global values
  global_weighted_value <- global_efficiency(weighted_graph)
  global_unweighted_value <- global_efficiency(unweighted_graph)
  #creating local values
  local_weighted_value <- average_local_efficiency(weighted_graph)
  local_unweighted_value <- average_local_efficiency(unweighted_graph)
  
  if (gender=='m'){
    male[nrow(male) + 1,] = c(global_weighted_value,global_unweighted_value,local_weighted_value,local_unweighted_value)
  print('male')
  }
  else if(gender == 'f'){
    female[nrow(female) + 1,] = c(global_weighted_value,global_unweighted_value,local_weighted_value,local_unweighted_value)
    print('female')
  }
}# patient gender

p_value_data[nrow(p_value_data)+1,] <- c(t.test(male["global_weighted"], female["global_weighted"])$p.value,
                                         t.test(male["local_weighted"], female["local_weighted"])$p.value,
                                         t.test(male["global_unweighted"], female["global_unweighted"])$p.value,
                                         t.test(male["local_unweighted"], female["local_unweighted"])$p.value)


boxplot(p_value_data)




#combine two dataframes. 



demographics <- as.data.frame(matrix_demo[2])
demo_age <- demographics[3]
median_age <- median(as.numeric(unlist(demo_age)))

old <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("global_weighted", "global_unweighted", "local_weighted", "local_unweighted"))))
young <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("global_weighted", "global_unweighted", "local_weighted", "local_unweighted"))))
p_value_data_age <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("global_w_p-value", "local_w_p-value", "global_u_p-value", "local_u_p-value"))))



for (i in 1:nrow(demo_sex)){
  age <- demo_age[i, 'demographics.age']
  c <- matrix_demo[[c(1, i)]] #current matrix
  c <- abs(c)
  
  #create weighted graph 
  weighted_graph <- graph_from_adjacency_matrix(c, weighted=T, mode="undirected", diag=F)
  
  #turn the data into binary
  c[c<0.5] <- 0
  c[c>0.5] <- 1
  #create unweighted graph. 
  unweighted_graph <- graph_from_adjacency_matrix(c, weighted=NULL, mode="undirected", diag=F)
  
  #creating global values
  global_weighted_value <- global_efficiency(weighted_graph)
  global_unweighted_value <- global_efficiency(unweighted_graph)
  #creating local values
  local_weighted_value <- average_local_efficiency(weighted_graph)
  local_unweighted_value <- average_local_efficiency(unweighted_graph)
  
  if (age >= median(as.numeric(unlist(demo_age)))){
    old[nrow(old) + 1,] = c(global_weighted_value,global_unweighted_value,local_weighted_value,local_unweighted_value)
    print('old')
  }
  else if(age < median(as.numeric(unlist(demo_age)))){
    young[nrow(young) + 1,] = c(global_weighted_value,global_unweighted_value,local_weighted_value,local_unweighted_value)
    print('young')
  }
}# patient gender

p_value_data_age[nrow(p_value_data_age)+1,] <- c(t.test(old["global_weighted"], young["global_weighted"])$p.value,
                                         t.test(old["local_weighted"], young["local_weighted"])$p.value,
                                         t.test(old["global_unweighted"], young["global_unweighted"])$p.value,
                                         t.test(old["local_unweighted"],young["local_unweighted"])$p.value)


boxplot(p_value_data_age)


summary(lm(demographics.age~demographics.site,demographics))




# Scatter plot
plot(x, y,
     pch = 19,
     col = factor(group))

# Legend
legend("topleft",
       legend = levels(factor(group)),
       pch = 19,
       col = factor(levels(factor(group))))




