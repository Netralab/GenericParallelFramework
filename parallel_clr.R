library(minet)
library(igraph)
library(doParallel)
library(ppcor)

#Specify the file path to obtain dataset
file_path = "address of dataset"
data_main = read.table(file = file_path, sep = '\t', header = TRUE)
col_name_ori = colnames(data_main[,-1]) #get column names
tot_no_of_col = ncol(data_main) #get number of columns
# done

#       Inference algorithms
#2. CLR
func_clr <- function(ds)
{
  result <- minet( ds, method="clr", estimator="spearman" )
  res <- result[,-1][-1,]
  for(row in 1:nrow(res))
  {
    for(col in 1:ncol(res)) 
    {
      if(res[row, col] > 0) 
      {res[row, col] <- 1}
      else {res[row, col] <- 0}
    }
  }
  return(res)
}
# CLR done #

#Recording start time
START =Sys.time()

# Setting up Parallel Environment 
numCores <- detectCores()
registerDoParallel(cores = numCores)
# done

# Splitting dataset
no_of_div = 10 #specify no. of divisions
col_in_each_div = (tot_no_of_col-1)/no_of_div #deciding the number of columns in each divisions
time_col = 1
count = 2

for(i in 1:no_of_div)
{
  temp_store = c(data_main[time_col]) #storing time column
  for(j in count:(count+col_in_each_div-1))
  {
    temp_store = cbind(temp_store,data_main[j])
  }
  count = count+col_in_each_div
  assign(paste("div", i, sep = "_"), as.data.frame(c(temp_store)))
}
list_of_div <- lapply( paste0('div_',1:no_of_div), get )
#done

# Initial inference # Note: function call has to be varied depending on inference algorithm #
# generating initial sub-matrices
list_of_submatrices = foreach(i=1:no_of_div)%dopar% {
  submatrix = as.data.frame(list_of_div[i])
  func_clr(submatrix)
}
# done

# Finding Hub genes in sub-matrices
hub_gene = foreach(i=1:no_of_div)%dopar% {
  sub_matrix = as.data.frame(as.data.frame(list_of_submatrices[i])) ##check
  row_sum_degree <- as.array(rowSums(sub_matrix))
  x=max(row_sum_degree)
  list_of_hub_genes = c()
  for(j in 1:nrow(row_sum_degree))
  {
    if(row_sum_degree[j]>=(x-1) & row_sum_degree[j] !=0 )
    {
      list_of_hub_genes = c(list_of_hub_genes,row_sum_degree[j])
    }
    hub_matrix = as.data.frame(list_of_hub_genes)
  }
  return(hub_matrix)
}
#done

# Appending hub genes data to the divisions
mod_div = foreach(i=1:no_of_div)%dopar%
  {
    re_div = as.data.frame(list_of_div[i])
    for(x in 1:no_of_div)
    {
      if(x==i)# excluding the div's hub genes
      {
        next
      }
      else
      {
        read_hub_genes = as.data.frame(hub_gene[x])
        if(nrow(read_hub_genes) == 0) # if hub gene file is empty then exclude
        {
          next
        }
        else
        {
          array_of_hub_genes=as.array(rownames(read_hub_genes)) #extract the name of the hub gene
          for(y in 1:length(array_of_hub_genes))#merging
          {
            gene_name = as.character(array_of_hub_genes[y])
            gene_data = c(data_main[gene_name])
            re_div[gene_name]=gene_data
          }
        }
      }#----else
    }
    #write appended data chunks  
    mod_div = as.data.frame(re_div)
    return(mod_div)
  }
# done

# Re-infering to get final sub-matrices # varies acc to inference algorithm
mod_matrix =foreach(i=1:no_of_div)%dopar%{
  re_div = as.data.frame(mod_div[i])
  res_matrix = func_clr(re_div)
  return(res_matrix) 
}
#done

# Merging the sub-matrices
dim = length(col_name_ori)
adj_matrix = matrix(0L, nrow = dim, ncol = dim)
colnames(adj_matrix) = col_name_ori
rownames(adj_matrix) = col_name_ori
for(a in 1:no_of_div)# Necessary Packages
{
  mod = as.data.frame(mod_matrix[a])
  mat = as.matrix(mod)
  node_name = colnames(mat)
  width = length(node_name)
  for(i in 1:width)
  {
    for(j in 1:width)
    {
      if(mat[node_name[i],node_name[j]] != 0)
      {
        adj_matrix[node_name[i],node_name[j]] = mat[node_name[i],node_name[j]]
        adj_matrix[node_name[j],node_name[i]] = mat[node_name[i],node_name[j]]
      }
    }
  }
}

END =Sys.time() # stop recording time
write.table(adj_matrix, file=paste("output folder address","clr_4000", sep = "_"), quote=FALSE, sep='\t')
# done

# Calculating time taken       
time_taken = END-START
print(time_taken)
#done
