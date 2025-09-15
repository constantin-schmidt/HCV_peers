# Number of CV data sets created
M <- 50

# Load the full list of panel data
load('./03_Data/ZZ_Temp/panel_het_blue.RData')

# Define a function to create data sets with random missing
generate_random_matrices <- function(mat, G, M) {
  
  # Predefine results list
  list_of_results <- vector("list", M)
  
  # Define function generating a CV data set
  set_missings_CV <- function(mat, G) {
    modified_mat <- mat
    g_col <- rep(NA, nrow(mat))
    
    for (i in 1:nrow(mat)) {
      if (G[i] < ncol(mat)) { # Check if G is less than the number of columns
        g_col[i] <- sample(1:ncol(mat), 1) # Randomly select a column from the whole row
        modified_mat[i, g_col[i]] <- NA # Set the selected column to missing
      }
    }
    
    return(list(modified_mat = modified_mat, g_col = g_col))
  } 
  
  # Generate CV data sets
  for (i in 1:M) {
    list_of_results[[i]] <- set_missings_CV(mat, G)
  }
  
  return(list_of_results)
}

# Generate a list with data sets
set.seed(42)
CV_het_data_blue <- generate_random_matrices(panel_list$Y, 
                                             panel_list$G, 
                                             M)

# Save the generated data sets
save(CV_het_data_blue, file = './03_Data/ZZ_Temp/CV_het_data_blue.RData')
