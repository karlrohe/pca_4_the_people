# reserved column names for features in this package:
# outcome, row_id, col_id, degree, weighted_degree, leverage, outcome_unweighted_1
# anything that starts with "pc_" or "vpc_"
# Avoid these as variable names in your data tibble.

library(tidyverse)
library(irlba)
library(Matrix)
library(gdim)
library(softImpute)

# pca_count / pca_sum 
# 

#
# diagnose(fo, tib)
# this looks at the degree distribution. 
# perhaps run this function before anything to make sure that there is enough data.
# 
# pick_dim(fo, tib, dimMax=20, num_bootstraps=2)
# Call Eigcv to get z-scores for each successive dimension. 
# 
# pca_count(fo, tib, k) / pca_sum(fo, tib, k)
# computes svd embeddings for normalized and regularized laplacian. 
#  this returns an "embed" object
#  row_features, 
#  column_features, 
#  middle_B, 
#  settings (a list with various details)
# 
# plot.embed computes multiple plots
# 
# Localization(embeddings) 
# diagnostic plot for Deg vs leverage. For both modes.  
# 
# streaks(embeddings,mode = "row", plot_columns= NULL)
# this makes a pairs plot
# 
# rotate(embeddings)
#  this returns an "embed" object after varimax rotation of rows and columns and constructing middle_B

############################################
##### The three core internal functions ####
############################################

# Here are the 3 "core internal functions" that you might like, but I do not imagine the typical user using. 
#  Given a formula and a tibble, we want to make a matrix (and other stuff)
#  fo = y~var1*(var2 & var3)
# 1) parse_variables(fo,tib); returns a list of length three. each element is a character vector of variable names.
#       first element is outcome variable.  second element is row variables.  third element is column variables.
#       "y"
#       "var1"
#       "var2" "var3"
# 2) make_edge_tib(fo, tib); this makes a tibble that contains all the necessary columns of tib. 
#       importantly, this function gives a row_id and col_id that become row/column indexes to define a matrix. 
#       this outputs a list of three elements, all tibbles: (edge_tib, row_index, col_index)
#       row/col indexes are kinda like row and column names, but allow for the fact that rows/columns might be defined by more than one element, e.g. month&day. 
# 3) make_sparse_matrix_raw(fo, tib); this outputs the sparse matrix to perform the svd upon. 
#       this outputs a list of three elements (A, row_universe, column_universe)
#       A is the sparse matrix.
#       row_universe and column_universe are identical to row_index and col_index in make_edge_tib


# parse variables calls the two functions:
# extract_interaction
# extract_variables

extract_interaction <- function(str) {
  
  pattern <- "\\(([^)]+)\\)|<([^>]+)>|\\b([a-zA-Z0-9_]+)\\b"
  #   This line defines a regular expression pattern to match specific parts of the input string.
  #   The pattern is composed of three parts, separated by the | (OR) operator, meaning it will match any one of these three patterns:
  #     \\(([^)]+)\\): Matches anything inside parentheses.
  # \\( and \\): Escaped parentheses to match literal parentheses in the string.
  # ([^)]+): A capturing group that matches any character except a closing parenthesis ()) one or more times (+).
  # <([^>]+)>: Matches anything inside angle brackets (< and >).
  # < and >: Matches literal angle brackets.
  # ([^>]+): A capturing group that matches any character except a closing angle bracket (>) one or more times.
  # \\b([a-zA-Z0-9_]+)\\b: Matches a standalone variable name.
  # \\b: Word boundary, ensuring we match whole words only.
  # [a-zA-Z0-9_]+: A capturing group that matches one or more alphanumeric characters or underscores.
  
  matches <- regmatches(str, gregexpr(pattern, str))
  # gregexpr(pattern, str): This function searches for all occurrences of the pattern in the input string str. It returns a list of match positions.
  # regmatches(str, ...): This function extracts the actual matched strings based on the positions found by gregexpr. The ... represents the list of match positions.
  return(matches[[1]])
}

extract_variables <- function(str) {
  # Remove parentheses if present
  str <- gsub("[()]", "", str)
  
  # Split the string by '|' or '*' and trim whitespace
  variables <- strsplit(str, "\\s*&\\s*")
  
  return(variables)
}


parse_variables <- function(fo, tib) {
  # this makes a list of length three. each element is a character vector of variable names.
  #      first element is outcome variable.  second element is row variables.  third element is column variables.
  
  # fo should be like this:
  # fo = val~(var1&var2) * (var4&var5 ) 
  # or
  # fo = val~ v1 * var2
  # tib should be a tibble the contains all the variable names. 
  
  # first, a simple check to make sure fo is a formula  check to make sure for is of the right form... 
  
  if(!is_formula(fo)){
    stop("first argument must be a formula; outcome ~ user_id * product_id")
  }
  
  # "left side" is the outcome...
  
  left_var <- all.vars(fo[[2]])
  # if it is an unweighted graph, then should 
  #  should have 1 on left hand side... 1 ~ userid*productid
  #  this makes all.vars(fo[[2]]) as character(0)
  #  also, we need to make a column in tib for it...
  if(length(left_var)==0){
    left_var="1"
    if("outcome_unweighted_1" %in% colnames(tib)){
      stop("outcome_unweighted_1 is 'reserved' column name for parse_variables with '1~' formulas. 
            this code is removing your column outcome_unweighted_1 and defining it as all 1s. this might cause an error later.")
    }
    tib = tib %>% mutate("outcome_unweighted_1" = 1)
    left_var = "outcome_unweighted_1"
  }
  
  # right side are the variables that define the features.
  
  
  right_expr <- fo[3]
  rhs = deparse(right_expr[[1]])
  rhs_terms = extract_interaction(rhs)
  if(length(rhs_terms)!=2){
    stop("Your input formula \n\n    ", 
         deparse(fo), 
         "\n\nhas been parsed to have these terms on the right hand side:\n\n    ", 
         paste0(rhs_terms,sep="\n    "), 
         "\nCurrently, this function only handles one interaction term on the right hand side of the formula. 
for example, these are valid: 

    y ~ v1*v3
    y ~ (v1 & v2) * v3
    y ~ (v1 & v2) * (v3 & v4 & v5)
    y ~ v1 * (v3 & v4).

However, these are not yet accepted:

    y ~ v1*v3 + v4 
    y ~ v1*v3 + v2*v4 
    y ~ v1*v2*v3")
  }
  right_vars = extract_variables(rhs_terms)
  
  
  if(left_var=="1") tib= tib %>% mutate("")
  
  # Combine left and right variables
  vars <- c(left_var, right_vars)
  
  # Check if the variables exist in the tibble
  
  missing_vars <- setdiff(unlist(vars), names(tib))
  if (length(missing_vars) > 0) {
    stop("The following variables are not in the tibble: ", paste(missing_vars, collapse = ", "))
  }
  
  return(vars)
  
}



make_edge_tib = function(fo, tib){
  vars = parse_variables(fo, tib)
  # Select the columns from the tibble
  
  
  outcome_column = vars[[1]]
  row_column = vars[[2]]
  column_column = vars[[3]]
  
  # row_index = tib %>% distinct(.data[[row_column]]) %>% mutate(row_id = row_number())
  # col_index = tib %>% distinct(.data[[column_column]]) 
  # 
  row_index = tib %>% 
    select(all_of(row_column)) %>% 
    drop_na() %>% 
    distinct() %>% 
    mutate(row_id = row_number())
  col_index = tib %>% 
    select(all_of(column_column)) %>% 
    drop_na() %>% 
    distinct() %>% 
    mutate(col_id = row_number())
  
  # if a weighted graph...
  if(outcome_column != "outcome_unweighted_1"){
    edge_list_tib = tib %>% select(all_of(unlist(vars)))
    
    edge_tib = edge_list_tib %>% 
      left_join(row_index, by = row_column) %>% 
      left_join(col_index, by = column_column) %>% 
      mutate(outcome= .data[[outcome_column]])
  }else{
    # if an un-weighted graph...
    edge_list_tib = tib %>% select(all_of(unlist(vars[-1]))) 
    
    edge_tib = edge_list_tib %>% 
      left_join(row_index, by = row_column) %>% 
      left_join(col_index, by = column_column) %>% 
      mutate(outcome= 1)
  }
  list(edge_tib = edge_tib %>% drop_na(), 
       row_index = row_index, 
       col_index = col_index)
}

make_sparse_matrix_raw = function(fo, tib, dropNA = TRUE){
  # This returns a list with elements
  #  A: sparse matrix for formula fo on tibble tib.
  #  row_universe: a tibble of distinct row-variable values, with the row_id (i.e. corresponding row number in A)
  #  col_universe: a tibble of distinct column-variable values, with the col_id (i.e. corresponding column number in A)
  
  vars = parse_variables(fo, tib)
  outcome_column = vars[[1]]
  row_column = vars[[2]]
  column_column = vars[[3]]
  
  edge_tib_list = make_edge_tib(fo,tib)
  
  edge_tib = edge_tib_list$edge_tib 
  if(dropNA) edge_tib = edge_tib %>% drop_na()
  row_universe = edge_tib_list$row_index
  column_universe = edge_tib_list$col_index
  A = sparseMatrix( 
    i = edge_tib$row_id, 
    j = edge_tib$col_id, 
    x = edge_tib$outcome,
    dims = c(nrow(row_universe), nrow(column_universe)))
  
  
  sp_A_dat = list(A=A, 
                  row_universe=row_universe, 
                  column_universe = column_universe)
  return(sp_A_dat)
}

make_incomplete_matrix_raw = function(fo, tib){
  
  vars = parse_variables(fo, tib)
  outcome_column = vars[[1]]
  row_column = vars[[2]]
  column_column = vars[[3]]
  
  edge_tib_list = make_edge_tib(fo,tib) 
  
  edge_tib = edge_tib_list$edge_tib
  mean_outcome = edge_tib %>% 
    group_by(row_id, col_id) %>% 
    summarize(outcome = mean(outcome, na.rm=T)) %>% ungroup
  
  column_means = mean_outcome %>% 
    group_by(col_id) %>% 
    summarize(col_mean = mean(outcome, na.rm=T)) %>% ungroup
  
  edge_tib2 = mean_outcome %>% 
    left_join(column_means) %>% 
    mutate(outcome = outcome - col_mean) %>% 
    select(-col_mean) %>% 
    left_join(
      edge_tib %>% 
        select(-outcome) %>% 
        select(all_of(c(row_column, column_column)), 
               row_id, col_id) %>% 
        distinct()
    )
  
  row_universe = edge_tib_list$row_index
  column_universe = edge_tib_list$col_index %>% 
    left_join(column_means) %>% 
    drop_na()
  
  Imat = softImpute::Incomplete(i = edge_tib$row_id, 
                                j = edge_tib$col_id,
                                x = edge_tib$outcome)
  
  incomplete_A_dat = list(A=Imat, 
                          row_universe=row_universe, 
                          column_universe = column_universe)
  return(incomplete_A_dat)
}



##### First, check that degrees are ok. 
### diagnose(fo, tib)
# two functions to help diagnose

itty_pivot = function(itty_tibby){
  # if you have a column like this:
  #   # A tibble: 162,263 × 2
  #   userId        degree
  #         <dbl>     <int>
  #         1           7
  #         2           36
  #         3          108
  # and you want it like this:
  #   # A tibble: 162,263 × 3
  #   type      id degree
  #   <chr>  <dbl>  <int>
  #   userId     1      7
  #   userId     2     36
  
  itty_tibby %>% 
    pivot_longer(1, names_to = "type", values_to = "id") %>% 
    relocate(type, id) %>% select(-id)
}
list_2_tib = function(list_of_data){
  lapply(list_of_data, itty_pivot) %>% bind_rows 
}
make_degrees = function(fo, tib){
  sparse_matrix_data = make_sparse_matrix_raw(fo, tib)
  row_degrees = tibble(row_id = 1:nrow(sparse_matrix_data$A), degree = rowSums(sparse_matrix_data$A!=0))
  col_degrees = tibble(col_id = 1:ncol(sparse_matrix_data$A), degree = colSums(sparse_matrix_data$A!=0))
  # edge_tib_list = make_edge_tib(fo,tib)
  # row_degrees = edge_tib_list$edge_tib %>% count(row_id) %>% rename(degree=n)
  # col_degrees = edge_tib_list$edge_tib %>% count(col_id) %>% rename(degree=n)
  list(row_degrees, col_degrees)
}
transpose_tibble <- function(data) {
  # Convert the first column to row names
  data <- tibble::column_to_rownames(data, var = names(data)[1])
  
  # Transpose the data frame
  data_transposed <- as.data.frame(t(data))
  
  # Convert row names back to a column
  data_transposed <- tibble::rownames_to_column(data_transposed, var = "measurement") 
  
  return(as_tibble(data_transposed))
}


diagnose = function(fo, tib, make_plot = TRUE){
  
  degrees_list = make_degrees(fo, tib)
  degrees_data = list_2_tib(degrees_list)
  # model_variables = parse_variables(fo,tib) %>% lapply(function(x) str_glue(x,sep = " & "))
  model_variables = parse_variables(fo,tib) %>% lapply(function(x) paste0(x,collapse = " & "))
  degrees_data = degrees_data %>% left_join(tibble(type=c("row_id", "col_id"), 
                                                   type_label = c(model_variables[[2]],model_variables[[3]]))) %>% 
    select(-type)
  
  if (make_plot) {
    p =  ggplot(degrees_data, aes(x = degree)) + 
      geom_histogram()+
      # geom_density(aes(y = after_stat(count))) + 
      scale_x_log10()+
      scale_y_log10()+
      facet_wrap(~type_label, scales="free")
    print(p)
  }
  
  dat = degrees_data %>% 
    group_by(type_label) %>% 
    summarize(number_of_items = n(),
              average_degree = round(mean(degree)),
              median_degree = median(degree),
              percent_le_1 = round(mean(degree<=1),2)*100,
              percent_le_2 = round(mean(degree<=2),2)*100,
              percent_le_3 = round(mean(degree<=3),2)*100)
  
  # dat %>% left_join(tibble(type=c("row_id", "col_id"), 
  #                          type_label = c(model_variables[[2]],model_variables[[3]]))) %>% 
  #   select(-type) %>% 
  #   relocate(type_label) %>% 
  dat %>%   transpose_tibble
  # transpose_tibble(dat)
  
}



# this computes the Z-scores for the cross-validated eigenvalues in 
warning_message_pick_dim = 
"the left side of the formula for pick_dim should be 1. 
            
Other variables on the left-hand-side are likely to create unreliable p-values. In many settings it likely still makes a good-bit of sense to have 1.  This code will still run, but you have been warned!"
pick_dim = function(fo, tib, dimMax=20, num_bootstraps=2){
  parsed_model =  parse_variables(fo, tib)
  outcome_variables = parsed_model[[1]]
  if(outcome_variables!="1"){
    warning(warning_message_pick_dim)
  }
  
  
  A = make_sparse_matrix_raw(fo, tib)$A
  gdim::eigcv(A = A, k_max = dimMax, laplacian= TRUE,num_bootstraps)
}


###############################
############  pca_count  ###########
###############################
# This is the key function that computes the "graph principal components"
#  in particular, for the rectangular matrix specified by fo and tibs_svd
#  it computes the leading k singular vectors (left and right)
#  of the normalized and regularized "graph Laplacian"

# glaplacian: this normalizes and regularizes the data matrix.

glaplacian <- function(A, regularize = TRUE) {
  # deg_row <- Matrix::rowSums(abs(A))
  # deg_col <- Matrix::colSums(abs(A))
  deg_row <- Matrix::rowSums(A!=0)
  deg_col <- Matrix::colSums(A!=0)
  
  tau_row <- mean(deg_row)
  tau_col <- mean(deg_col)
  
  D_row = Diagonal(nrow(A), 1 / sqrt(deg_row + tau_row))
  D_col = Diagonal(ncol(A), 1 / sqrt(deg_col + tau_col))
  L = D_row %*% A %*% D_col
  rownames(L) = rownames(A)
  colnames(L) = colnames(A)
  return(L)
}


generate_colnames <- function(matrix_data, prefix) {
  num_cols = ncol(matrix_data)
  num_digits = nchar(as.character(num_cols))
  
  formatted_colnames = sprintf(paste0(prefix, "_%0", num_digits, "d"), 1:num_cols)
  return(formatted_colnames)
}

s_2_embedding = function(sparse_matrix_data, s, factor_prefix){
  # sparse_matrix_data is a list with:
  A = sparse_matrix_data$A
  row_universe = sparse_matrix_data$row_universe
  column_universe = sparse_matrix_data$column_universe
  
  #  s is output of svd with s$u, s$d, s$v.
  # factor_prefix is a character string that will be put at the front (e.g. "pc" for pc_1,...)
  
  #  we want to make a tidy output.
  
  
  # sometimes first dimension is all negative and that's annoying...
  #  make it mostly positive...
  if(median(s$u[,1])<0){
    s$u[,1]= -s$u[,1]
    s$v[,1]= -s$v[,1]
  }
  u = s$u  * sqrt(nrow(A))
  v = s$v  * sqrt(ncol(A))
  
  #  I currently naming dimensions/factors pc_1, pc_2, ... for both row and column embedding.
  #    this has some risks. 
  
  colnames(u) = generate_colnames(u, factor_prefix)
  colnames(u) = paste0(colnames(u), "_rows")
  colnames(v) = generate_colnames(v, factor_prefix)
  colnames(v) = paste0(colnames(v), "_columns")
  
  
  keep_these_rows = rowSums(abs(A))>0
  row_features = bind_cols(row_universe, 
                           degree= rowSums(A!=0)[keep_these_rows],
                           weighted_degree = rowSums(abs(A))[keep_these_rows],
                           as_tibble(u[keep_these_rows,]))
  keep_these_cols = colSums(abs(A))>0
  column_features = bind_cols(column_universe, 
                              degree= colSums(A!=0)[keep_these_cols],
                              weighted_degree = colSums(abs(A))[keep_these_cols],
                              as_tibble(v[keep_these_cols,]))
  
  
  middle_B= bind_cols(row_factors = colnames(u), 
                      column_factors = colnames(v), 
                      value = s$d / (sqrt(nrow(A)*ncol(A))))
  list(row_features = row_features, column_features= column_features, middle_B=middle_B)
}



pca_count = function(fo, tib, k){
  
  # this is the key user function to generate embeddings. 
  # the output is a list of:
  # row_features (for whichever term is first on the right hand side of the formula)
  # column_features (for whichever term is second on the right hand side of the formula)
  # middle_B (this is a diagonal matrix in triplet for, stored in a tibble)
  # settings (this is a list of details)
  
  
  
  
  sp_A_dat = make_sparse_matrix_raw(fo, tib)
  A = sp_A_dat$A
  A@x = sqrt(A@x)
  L = glaplacian(A)
  s_svd = irlba(L,nu = k, nv = k)
  factor_prefix = "pc"
  
  
  
  embeds = s_2_embedding(sparse_matrix_data = sp_A_dat, s = s_svd, factor_prefix=factor_prefix)
                           
  
  parsed_model =  parse_variables(fo, tib)
  
  settings = list(fit_method = "pca_count", 
                  prefix_for_embeddings = str_glue(factor_prefix, "_"), 
                  k = k,
                  normalized = TRUE,
                  reguarlized = TRUE,
                  outcome_variables = parsed_model[[1]],
                  row_variables  = parsed_model[[2]],
                  column_variables  = parsed_model[[3]])
  
  embeds[[4]] = settings
  names(embeds)[4] = "settings"

  class(embeds) = "embed"
  embeds
  
  
}
pca_sum = pca_count







remove_L_normalization = function(s_svd, A, orthogonalize= FALSE){
  # if using L in  s_svd(L) ,
  #  then we might want to remove that normalization in the output. 
  #  this function does that.
  deg_row <- Matrix::rowSums(abs(A))
  deg_col <- Matrix::colSums(abs(A))
  
  tau_row <- mean(deg_row)
  tau_col <- mean(deg_col)
  
  D_row = Diagonal(nrow(A), sqrt(deg_row + tau_row))
  D_col = Diagonal(ncol(A), sqrt(deg_col + tau_col))
  s_svd$u = D_row %*% s_svd$u
  s_svd$v = D_col %*% s_svd$v
  s_svd$u = as.matrix(s_svd$u)
  s_svd$v = as.matrix(s_svd$v)
  
  if(!orthogonalize) return(s_svd)
  # if we want to re-orthogonalize s$u and s$v....
  #  then you need to do some algebra.
  #  I think this is right, but has not been tested:
  su= svd(s_svd$u)
  sv= svd(s_svd$v)
  b= diag(su$d) %*% t(su$v) %*% diag(s_svd$d) %*% sv$v %*% diag(sv$d)
  sb = svd(b)
  u = su$u%*% sb$u
  d = sb$d
  v = sv$u %*% sb$v
  return(list(u = as.matrix(u), d = d, v = as.matrix(v)))
}

# pca_average does low-rank matrix completion


pca_average = function(fo, tib, k){
  
  # this is the second user function to generate embeddings. 
  # the output is an embed object... a list of:
  # row_features (for whichever term is first on the right hand side of the formula)
  # column_features (for whichever term is second on the right hand side of the formula)
  # middle_B (this is a diagonal matrix, but in a triplet form in a tibble)
  # settings (this is a list of details)
  
  
  sp_A_dat = make_incomplete_matrix_raw(fo, tib)
  A = sp_A_dat$A
  L = glaplacian(A)
  s_svd = softImpute::softImpute(L,rank.max = k)
  s_svd = remove_L_normalization(s_svd,A)
  
  factor_prefix = "na_pc"
  
  embeds = s_2_embedding(sparse_matrix_data = sp_A_dat, s = s_svd, factor_prefix=factor_prefix)

  
  parsed_model =  parse_variables(fo, tib)
  
  settings = list(fit_method = "pca_average", 
                  prefix_for_embeddings = str_glue(factor_prefix, "_"), 
                  k = k,
                  normalized = TRUE,
                  reguarlized = TRUE,
                  outcome_variables = parsed_model[[1]],
                  row_variables  = parsed_model[[2]],
                  column_variables  = parsed_model[[3]])
  
  embeds[[4]] = settings
  names(embeds)[4] = "settings"
  
  class(embeds) = "embed"
  embeds
  
}
pca_mean = pca_average



# make_leverage is called in "localization"
make_leverage = function(embeddings){
  # types_of_modes = embeddings %>% distinct(type) %>% pull(type) 
  # this_type = types_of_modes[1]
  lev_list=list()
  
  row_embedding_mat = embeddings$row_features  %>% 
    select(starts_with(embeddings$settings$prefix_for_embeddings)) %>% 
    as.matrix 
  
  row_leverage = rowSums(row_embedding_mat^2)
  
  lev_list[[1]] = embeddings$row_features  %>% 
    select(-starts_with(embeddings$settings$prefix_for_embeddings)) %>% 
    bind_cols(tibble(leverage = row_leverage))
  
  col_embedding_mat = embeddings$column_features  %>% 
    select(starts_with(embeddings$settings$prefix_for_embeddings)) %>% 
    as.matrix 
  
  col_leverage = rowSums(col_embedding_mat^2)
  
  lev_list[[2]] = embeddings$column_features  %>% 
    select(-starts_with(embeddings$settings$prefix_for_embeddings)) %>% 
    bind_cols(tibble(leverage = col_leverage))
  
  
  names(lev_list) = names(embeddings)[1:2]
  lev_list
  
}
make_deg_lev_resid = function(lev_tib){
  fit = lm(I(log(leverage+.00000001))~I(log(degree+1)),  data = lev_tib)
  bind_cols(lev_tib , tibble(residuals = fit$resid))
}
localization = function(embeddings){
  # plot degree vs leverage...
  #   actually, the residual of log(leverage)~log(degree)
  #   against log(degree).
  
  
  leverages_data = make_leverage(embeddings) 
  
  
  # make mode labels
  my_row_names = embeddings$settings$row_variables %>% paste(collapse  = " & ")
  my_col_names = embeddings$settings$column_variables %>% paste(collapse  = " & ")
  
  
  row_dat = leverages_data[[1]] %>% 
    select(degree,leverage) %>% 
    mutate(mode = my_row_names) %>% 
    make_deg_lev_resid
  
  col_dat = leverages_data[[2]] %>% 
    select(degree,leverage) %>% 
    mutate(mode = my_col_names) %>% 
    make_deg_lev_resid
  
  deg_lev_resid = bind_rows(row_dat, col_dat)
  
  
  deg_lev_resid %>% 
    ggplot(aes(x=degree+1, y= residuals)) +
    geom_point()+
    # ggtitle(paste("   ", round(fit$coefficients[1],2),  " x deg ^", round(fit$coefficients[2],2)))+
    geom_smooth(se = F, method = "gam", data=subset(deg_lev_resid, degree >=2))+
    facet_wrap(~mode, scales = "free")+
    scale_x_log10()   
}


# these next two functions help the function streaks
pair = function(u, n = 1000){
  if(nrow(u)>1000){
    lev = rowSums(u^2)
    samp = sample(nrow(u), n, prob = lev)
  }else{samp= 1:nrow(u)}
  par(bty="n")
  pairs(u[samp,], panel =my_panel)
  par(bty="o")
}
my_panel = function(x,y){
  
  points(x,y,pch ='.')
  abline(h=0, col='red'); abline(v=0, col='red')
  points(0,0, pch = "O", col = "red")
  
}
streaks = function(embeddings, mode = "rows",plot_columns= NULL){
  #  pairs plot
  #  type_mode plots the 
  # make plot to look for radial streaks
  # if(is.null(type_mode)) type_mode = default_mode(embeddings)
  #   if(!is_tibble(embedding)){
  #     stop(
  #       "error: you need to give it a tibble. 
  # this would happen if you gave it an embed object
  # Instead, you need to pick either $row_features or $column_features.
  # ")
  #   }
  if(mode == "rows"){
    embedding = embeddings$row_features %>% select(starts_with(embeddings$settings$prefix_for_embeddings))
  }else{
    embedding = embeddings$column_features %>% select(starts_with(embeddings$settings$prefix_for_embeddings))
  }
  
  embedding_mat = embedding %>% 
    as.matrix 
  # scale(center=F)
  
  k = ncol(embedding_mat)
  
  # identify which columns to plot if plot_columns is NULL
  if(is.null(plot_columns)){
    if(k<=10){plot_columns = 1:k}else{
      plot_columns = c(1:5, (k-4):k)}
  }
  embedding_mat = embedding_mat[,plot_columns]
  nn = nrow(embedding_mat)
  if(nn>1000){
    levs = rowSums(embedding_mat^2)
    samp = sample(nrow(embedding_mat), size = 1000,prob = levs)
  }else{samp = 1:nn}
  
  pair(embedding_mat[samp,])
}



plot.embed = function(embeddings){
  B = get_middle_matrix(embeddings)
  k = nrow(B)
  par(mfrow = c(1,1))
  plot(diag(B), main = "Singular values of L (biased)")
  # , ylim = c(0, max(diag(B))))
  # cutoff1 = (sqrt(1/mean(embeddings$row_features$degree)) + sqrt(1/mean(embeddings$column_features$degree)))/2
  # cutoff2 = (sqrt(1/mean(embeddings$row_features$weighted_degree)) + sqrt(1/mean(embeddings$column_features$weighted_degree)))/2
  # lines(x = c(1,max(dim(B))), y = cutoff1*c(1,1))
  # lines(x = c(1,max(dim(B))), y = cutoff2*c(1,1))
  readline(prompt="Press [Enter] to continue to the next plot...")
  
  skip_first_singular_values = diag(B)[-1]
  plot(2:k, skip_first_singular_values, main = "2:k singular values of L (biased)")
  # , ylim = c(0, max(diag(B)[-1])))
  readline(prompt="Press [Enter] to continue to the next plot...")
  
  
  localization_plot = localization(embeddings)
  print(localization_plot)
  
  readline(prompt="Press [Enter] to continue to the next plot...")
  
  streaks(embeddings) 
  readline(prompt="Press [Enter] to continue to the next plot...")
  
  streaks(embeddings, "cols")
  # readline(prompt="Press [Enter] to continue to the next plot...")
  # 
  # image(B, main = "middle B matrix")
}

rotate = function(embeddings){
  # rotate both modes with variamax.  return an embedding object.
  #  we want the axes for rows and columns to roughly align with each other 
  #    (makes other interpretation easier... tends to make B have strong diagonal).
  #  so, rotate the smaller one first... then use that rotation to pre-rotate the bigger one...
  #           .... then take that rotation back to pre-rotate the smaller one.
  #
  if(nrow(embeddings[[1]]) < nrow(embeddings[[2]])) mode_order_vector = c(1,2)
  if(nrow(embeddings[[1]]) > nrow(embeddings[[2]])) mode_order_vector = c(2,1)
  
  
  # rotate the smaller one...
  
  embedding_mat_mode1 = embeddings[[mode_order_vector[1]]] %>% 
    select(starts_with(embeddings$settings$prefix_for_embeddings)) %>% 
    as.matrix()
  
  varimax_rotation_mode1 = varimax(embedding_mat_mode1,normalize = F)$rot
  
  # pre-rotate the bigger one...
  
  embedding_mat_mode2 = embeddings[[mode_order_vector[2]]] %>% 
    select(starts_with(embeddings$settings$prefix_for_embeddings)) %>% 
    as.matrix()
  
  
  # # pre-rotate with last varimax rotation:
  v_mode_2 = varimax_with_pre_rotation(embedding_mat_mode2, varimax_rotation_mode1)
  sparse_embedding_mat_mode2 = v_mode_2$data_after_rotation
  
  # pre-rotate first mode with both of the previous varimax rotations.
  v_mode_1 = varimax_with_pre_rotation(embedding_mat_mode1, v_mode_2$varimax_rotation)
  sparse_embedding_mat_mode1 = v_mode_1$data_after_rotation
  
  old_prefix = embeddings$settings$prefix_for_embeddings %>% str_remove("_")
  mode1_output = make_sparse_output(sparse_embedding_mat_mode1,
                                    v_mode_1$varimax_rotation,
                                    old_prefix = old_prefix) 
  mode2_output = make_sparse_output(sparse_embedding_mat_mode2,
                                    v_mode_2$varimax_rotation,
                                    old_prefix = old_prefix) 
  
  
  #if we did 2 first, then 1... then we need to flop mode1_output and mode2_output.
  if(mode_order_vector[1] == 2){
    row_output = mode2_output
    column_output = mode1_output
  }
  if(mode_order_vector[1] == 1){ # if we did 1 first, then no swapping
    row_output = mode1_output
    column_output = mode2_output
  }
  row_features = bind_cols(
    embeddings$row_features %>% select(-starts_with(embeddings$settings$prefix_for_embeddings)),
    row_output$sparse_tib)
  
  column_features = bind_cols(
    embeddings$column_features %>% select(-starts_with(embeddings$settings$prefix_for_embeddings)),
    column_output$sparse_tib)
  
  
  # construct the new middle B matrix by passing through rotations from:
  #  this should likely be it's own function at some point.... would be useful
  #  to abstract away this construction of B 
  # row_output$rot_mat
  # mode2_output$rot_mat  
  # embeddings$middle_B
  
  # first, $middle_B is stored as "edge list form".  So, convert to a matrix....
  
  old_B = get_middle_matrix(embeddings)
  
  new_B = t(row_output$rot_mat) %*% old_B %*% column_output$rot_mat
  
  
  
  B_tib= make_middle_B_tibble(new_B, factor_prefix = "")
  
  new_settings = c(
    list(
      fit_method = str_glue(embeddings$settings$fit_method, "varimax",.sep = " + "),
      prefix_for_embeddings = str_glue("v",old_prefix,"_")
    ),
    embeddings$settings[-(1:2)])
  
  out_list = list(row_features = row_features,
                  column_features = column_features,
                  middle_B = B_tib,
                  settings = new_settings)
  
  class(out_list) = "embed"
  out_list
  
}

get_middle_matrix = function(embeddings){
  b_el = embeddings$middle_B
  B_list = make_sparse_matrix_raw(value~row_factors*column_factors, tib = b_el)
  B_matrix = B_list$A
  rownames(B_matrix) = B_list$row_universe$row_factors
  colnames(B_matrix) = B_list$column_universe$column_factors
  return(B_matrix)
}
make_middle_B_tibble = function(B_matrix, factor_prefix){
  colnames(B_matrix) = paste0(factor_prefix,1:ncol(B_matrix))
  B_tib = as_tibble(as.matrix(B_matrix)) %>% 
    mutate(row_factors = paste0(factor_prefix,1:nrow(B_matrix))) %>% 
    relocate(row_factors) %>% 
    pivot_longer(-1, names_to="column_factors", values_to = "value")
  return(B_tib)
  
}

# these next two functions help with the function rotate.

make_sparse_output = function(sparse_embed_mat, rot_mat, old_prefix){
  #  used in rotate.
  # given the output of varimax rotated matrix... we want to
  # 1) make skew positive.
  # 2) give nice column names.
  # 3) make it a tibble
  # 4) output the rotation matrix (that make it skew positive)
  
  # get columns skew sign:
  skew_sign = apply(sparse_embed_mat, 2, function(x) sign(sum(x^3)))
  #  multiply each column of sparse_embed_mat by corresponding skew_sign element.
  sparse_embed_mat_output = sweep(sparse_embed_mat, MARGIN = 2, STATS = skew_sign, FUN = "*")
  # same for the rotation
  rot_mat_output = sweep(rot_mat, MARGIN = 2, STATS = skew_sign, FUN = "*")
  
  # nice column names
  colnames(sparse_embed_mat_output) = generate_colnames(sparse_embed_mat_output, str_glue("v",old_prefix))
  list(sparse_tib = as_tibble(sparse_embed_mat_output), rot_mat = rot_mat_output)
}

varimax_with_pre_rotation = function(matrix_to_rotate, pre_rotation){
  # used in function rotate.
  #  varimax has multiple optima.  it is nice for rotation of both modes to find
  #  ~analogous~ modes of varimax optima.  we do that we pre-rotation.
  # before computing varimax on an embedding, pre-rotate with an optima from the other mode.
  #  then, iterate from there.  empirically, this should speed up convergence and improve interpretability.
  matrix_to_rotate_pre_rotated = matrix_to_rotate %*% pre_rotation
  next_varimax_rotation = varimax(matrix_to_rotate_pre_rotated, normalize = F)$rot
  list(data_after_rotation = matrix_to_rotate_pre_rotated %*% next_varimax_rotation,
       varimax_rotation = pre_rotation%*%next_varimax_rotation)
}



# 
# 
# TODO:
# predict = function(embeddings, tib_test){
#   # Add  implementation details here
# }
# 
# 
# 
# Interpret = function(embeddings, features){
#   # Add implementation details here
#   bff
# }
# 
# tree = function(sparse_embeddings){
#   # Add implementation details here
# }
# bff_with_words = function(loadings, strings, num_best = 10){
#   text_df <- tibble(id = 1:length(strings),
#                     text = strings)
#   
#   # this does a lot of processing!
#   #  to lower, remove @ # , .
#   #  often these make sense on a first cut.
#   #  worth revisiting before "final results"!
#   tt  = text_df %>% unnest_tokens(word, text)
#   dt = cast_sparse(tt, id, word)
#   cs = colSums(dt)
#   dt = dt[,cs>3]
#   these_have_words = tt$id %>% unique %>% sort
#   vsp::bff(as.matrix(loadings)[these_have_words,], dt, num_best = num_best)
# }
# 
# 
# 
# top_elements = function(embeddings,dimension,mode="columns"){
#   if(mode == "columns"){
#     embeddings$column_features %>% arrange(abs())
#     apply(loadings,2, function(x) names[order(-x)[1:10]])
#   }
# }
# 
# top_loadings
# 

