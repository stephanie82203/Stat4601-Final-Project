# Author: Stephanie Cheng

######### PCA ###############
pca <- function(df,df_name) {
  # Apply PCA
  df_pca <- prcomp(df, scale. = TRUE)
  
  # Print PCA summary and the PCs' variables name
  cat("\n===== PCA Summary for", df_name, "=====\n")
  print(summary(df_pca))
  cat("Contributing variable for each PC:\n")
  variables <- apply(df_pca$rotation, 2, function(x) {
    names(x)[which.max(abs(x))]
  })
  print(variables)
  colnames(df_pca$x) <- variables
  return(df_pca)
}

######### k-Mean ############

# This code is from Section 6 - vq.helpers.r
Davies.Bouldin <- function(A, SS, m) {
  # A  - the centres of the clusters
  # SS - the within sum of squares
  # m  - the sizes of the clusters
  N <- nrow(A)   # number of clusters
  # intercluster distance
  S <- sqrt(SS/m)
  # Get the distances between centres
  M <- as.matrix(dist(A))
  # Get the ratio of intercluster/centre.dist
  R <- matrix(0, N, N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      R[i,j] <- (S[i] + S[j])/M[i,j]
      R[j,i] <- R[i,j]
    }
  }
  return(mean(apply(R, 1, max)))
}

# Calculate DBI and WSS for PCA
calculate_k_stats_PCA <- function(pca_result, max_k) {
  X.syn <- pca_result$x[, 1:1:ncol(pca_result$x)] 
  DBI <- numeric(max_k - 1)
  errs <- numeric(max_k - 1)

  for (k in 2:max_k) {
    set.seed(12345)
    KM <- kmeans(X.syn, centers = k, nstart = 15, iter.max = 100)
    errs[k - 1] <- sum(KM$withinss)
    DBI[k - 1] <- Davies.Bouldin(KM$centers, KM$withinss, KM$size)
  }
  return(list(DBI = DBI, errs = errs, X.syn = X.syn))
}

# Calculate DBI and WSS for standard dataframe
calculate_k_stats <- function(data, max_k) {
  DBI <- numeric(max_k - 1)
  errs <- numeric(max_k - 1)

  for (k in 2:max_k) {
    set.seed(12345)
    KM <- kmeans(data, centers = k, nstart = 15, iter.max = 100)
    errs[k - 1] <- sum(KM$withinss)
    DBI[k - 1] <- Davies.Bouldin(KM$centers, KM$withinss, KM$size)
  }
  return(list(DBI = DBI, errs = errs))
}

# Show the DBI and the Sum of Squares 
plot_kmeans <- function(errs, DBI) {
  oldpar <- par(mfrow = c(1, 2))
  err_k <- 2:(length(errs) + 1)
  dbi_k <- 2:(length(DBI) + 1)

  ## ---- WSS Plot ----
  plot(err_k, errs, type = "b", pch = 17, main = "Within-Cluster Sum of Squares", 
       xlab = "k", ylab = "WSS")
  
  wss_lm <- lm(errs ~ err_k)
  abline(wss_lm, col = "blue", lty = 2)

  ## ---- DBI Plot ----
  plot(dbi_k, DBI, type = "b", pch = 17, main = "Davies-Bouldin Index", 
       xlab = "k", ylab = "DBI")
  
  par(oldpar)
}

# Clustering Plot
plot_clusters <- function(X.syn, min_k , max_k) {
  set.seed(12345)

  oldpar <- par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
  colors <- rainbow(max_k + 1)
  
  for (k in min_k:max_k) {
    set.seed(12345)
    KM <- kmeans(X.syn, k, 15)
    
    # Ensure pch values are within valid range (Ref: ChatGpt)
    valid_pch <- 0:25
    pch_vals <- valid_pch[(KM$cluster %% length(valid_pch)) + 1]
    
    plot(X.syn[, 1:2], col = colors[KM$cluster], pch = pch_vals,
         main = paste(k, "clusters"), 
         xlab = "PC1", 
         ylab = "PC2"[2])
  } 
  par(oldpar)
}