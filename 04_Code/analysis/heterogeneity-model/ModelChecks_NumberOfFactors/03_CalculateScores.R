####################
##  Load the data ##
####################

scores_list <- list()
for (i in 1:50) {
  file_name <- paste0("./03_Data/zz_CVoutput/scores_", i, ".rds")
  load(file_name)
  scores_list[[i]] <- get("score_matrix")
}

############################
##  Calculate overall MSE ##
############################

MSE_df <- do.call(rbind, lapply(scores_list, function(x) x[, 1])) |>
  as.data.frame()
MSE_df[] <- lapply(MSE_df, as.numeric)

MSE <- colMeans(MSE_df)
names(MSE) <- c("0", "1", "2", "3", "4", "5")
print("Mean squared errors by number of factors (k)")
print(MSE)

#######################################
##  Calculate overall interval score ##
#######################################

interval_df <- do.call(rbind, lapply(scores_list, function(x) x[, 2])) |>
  as.data.frame()
interval_df[] <- lapply(interval_df, as.numeric)

interval_scores <- colMeans(interval_df)
names(interval_scores) <- c("0", "1", "2", "3", "4", "5")
print("Interval scores by number of factors (k)")
print(interval_scores)
