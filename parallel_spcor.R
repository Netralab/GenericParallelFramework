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
#SPCOR
spcor <- function(x, method = c("pearson", "kendall", "spearman"))
{
  # correlation method
  method <- match.arg(method)
  
  # check the data
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x)) 
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  
  # sample number
  n <- dim(x)[1]
  
  # given variables' number
  gp <- dim(x)[2]-2
  
  # covariance matrix
  cvx <- cov(x,method=method)
  
  # inverse covariance matrix
  if(det(cvx) < .Machine$double.eps){
    warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
    icvx <- ginv(cvx)
  }else
    icvx <- solve(cvx)
  
  # semi-partial correaltion
  spcor <- -cov2cor(icvx)/sqrt(diag(cvx))/sqrt(abs(diag(icvx)-t(t(icvx^2)/diag(icvx))))
  diag(spcor) <- 1
  
  # p-value
  if(method == "kendall"){
    statistic <- spcor/sqrt(2*(2*(n-gp)+5)/(9*(n-gp)*(n-1-gp)))
    p.value <- 2*pnorm(-abs(statistic))
    
  }else{
    statistic <- spcor*sqrt((n-2-gp)/(1-spcor^2))
    p.value <- 2*pt(-abs(statistic),(n-2-gp))
    #p.value <- 2*pnorm(-abs(statistic))
  }
  
  diag(statistic) <- 0
  diag(p.value) <- 0
  
  list(estimate=spcor,p.value=p.value,statistic=statistic,n=n,gp=gp,method=method)
}

func_spcor <- function(ds)
{
  #ds <- ds[,-1]
  col_name_ds = colnames(ds[,-1]) #get column names
  res1 <- spcor( ds, method="pearson" )
  result <- as.data.frame(res1[1])
  res <- result[,-1][-1,]
  for(row in 1:nrow(res))
  {
    for(col in 1:ncol(res)) 
    {
      if (res[row, col] == 0)
      {res[row, col] <- 0}
      else if ((res[row, col] > 0.05) || (res[row, col] < -0.05))
      {res[row, col] <- 1}
      else {res[row, col] <- 0}
    }
  }
  colnames(res) = col_name_ds
  rownames(res) = col_name_ds
  return(res)
}
# SPCOR DONE#######


#Recording start time
START =Sys.time()

# Setting up Parallel Environment 
numCores <- detectCores()
registerDoParallel(cores = numCores)
# done

# Splitting dataset
no_of_div = 100 #specify no. of divisions
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
  func_spcor(submatrix)
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
  res_matrix = func_spcor(re_div)
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
write.table(adj_matrix, file=paste("output folder address","spcor_4000", sep = "_"), quote=FALSE, sep='\t')
# done

# Calculating time taken       
time_taken = END-START
print(time_taken)
#done
