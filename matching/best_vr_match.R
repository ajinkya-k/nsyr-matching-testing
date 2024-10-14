# best_block_vr_match

# Loop over the blocks getting the matched_set_vectors
# pull out the vector indicating matched_set id. this is currently saved as tmp_matched_set
# paste in the block number to these elements


best_vr_match_block <- function(prop_score_df, out_data, X_mat_full, max_K = 5, block_subset = list(),cont_vars = c(), cat_vars = c(), mix_vars = c(), mix_values = list(), verbose = FALSE){
  # make a matrix to count the imbalances
  imbalance_counts <- matrix(nrow = max_K-1, ncol = 6, dimnames = list(c(), c("K","Pre-match", "Post-match", "Pre-match.missing", "Post-match.missing", "n")))
  imbalance_counts[,"K"] <- 2:max_K
  z_all <- prop_score_df$treated
  names(z_all) <- rownames(prop_score_df)
  for(K in 2:max_K){
    print(paste("Starting K = ", K, "at", Sys.time()))
    matched_sets_all <- c()
    for(block in unique(prop_score_df[,"block"])){
      #print(paste("Starting block", block, "at", Sys.time()))
      # subset to only those members of the block
      tmp_prop_score_df <- prop_score_df[prop_score_df[,"block"] == block,]
      z <- tmp_prop_score_df[,"treated"]
      prop_score <- tmp_prop_score_df[,"Propensity"]
      logit_prop <- tmp_prop_score_df[,"Logit(Propensity)"]
      names(z) <- rownames(tmp_prop_score_df)
      names(prop_score) <- rownames(tmp_prop_score_df)
      names(logit_prop) <- rownames(tmp_prop_score_df)
      # compute exclude cs
      exclude_cs_treated <- names(z)[z == 1 & prop_score > max(prop_score[z == 0])]
      exclude_cs_control <- names(z)[z == 0 & prop_score < min(prop_score[z == 1])]
      exclude_cs <- c(exclude_cs_treated, exclude_cs_control)
      
      # pull out the appropriate submatrix of X_mat
      X_mat <- X_mat_full[rownames(tmp_prop_score_df),]
      X_mat <- X_mat[!rownames(X_mat) %in% exclude_cs,]
      z <- z[!names(z) %in% exclude_cs]
      prop_score <- prop_score[!names(prop_score) %in% exclude_cs]
      logit_prop <- logit_prop[!names(logit_prop) %in% exclude_cs]
      
      subject_index <- seq(1, length(z))
      
      drop_cols <- c()
      for(x in colnames(X_mat)){
        if(var(X_mat[,x]) == 0){
          drop_cols <- c(drop_cols,x)
        }
      }
      X_mat <- X_mat[,!colnames(X_mat) %in% drop_cols]
      tmp_dist_mat <- smahal(z, X_mat)
      dist_mat <- addcaliper(tmp_dist_mat, z, logit_prop) # adds  propensity score caliper
      rownames(dist_mat)<- subject_index[z == 1]
      subject_id <- names(z)
      
      # do all the matching
      match_vec_list <- list()
      for(kk in 1:K){
        if(kk == 1){
          treated_ind <- which(z == 1 & prop_score > 1/3)
          control_ind <- which(z == 0 & prop_score > 1/3)
        } else if(kk == K){
          treated_ind <- which(z == 1 & prop_score < 1/(K+1))
          control_ind <- which(z == 0 & prop_score < 1/(K+1))
        } else{
          treated_ind <- which(z == 1 & prop_score > 1/(kk+2) & prop_score < 1/(kk+1))
          control_ind <- which(z == 0 & prop_score > 1/(kk+2) & prop_score < 1/(kk+1))
        }
        row_ind <- rownames(dist_mat) %in% treated_ind
        col_ind <- colnames(dist_mat) %in% control_ind
        n_treated <- sum(row_ind)
        n_control <- sum(col_ind)
        
        #print(paste("kk = ", kk, "n.treated =", n.treated, "n.control =", n.control))
        if(n_control == 1 & n_treated == 1){
          # matching software doesn't handle this well so let's ignore for now
        } else if(n_control >= n_treated & n_treated > 0){
          num_controls <- min(floor(n_control/n_treated), kk)
          tmp_dist_mat <- matrix(dist_mat[row_ind, col_ind], nrow = n_treated, ncol = n_control, dimnames = list(rownames(dist_mat)[row_ind], colnames(dist_mat)[col_ind]))
          match_vec_list[[paste0("kk.",kk)]] <- pairmatch(tmp_dist_mat, controls = num_controls)
        } else if(n_control < n_treated & n_control > 0){ # relabel controls as treated and optimally exclude treated units
          tmp_dist_mat <- t(matrix(dist_mat[row_ind, col_ind], nrow = n_treated, ncol = n_control, dimnames = list(rownames(dist_mat)[row_ind], colnames(dist_mat)[col_ind])))
          match_vec_list[[paste0("kk.",kk)]] <- pairmatch(tmp_dist_mat)
        } else if(n_treated == 0 | n_control == 0){
          #print("No controls or no treated. No match possible in this stratum")
        }
      } # closes loop over kk
      tmp_matched_sets <- stratified.matched.set(match_vec_list, z)
      tmp_matched_sets[exclude_cs] <- NA
      tmp_matched_sets[!is.na(tmp_matched_sets)] <- paste0(block, ".", tmp_matched_sets[!is.na(tmp_matched_sets)])
      matched_sets_all <- c(matched_sets_all, tmp_matched_sets)
    } # closes loop over the blocks
    
    # make the data frame
    #print("    making the data frame")
    tmp_match <- stratified.df(prop_score_df, out_data, matched_sets_all, z_all)
    # assess balance
    tmp_balance <- assess_balance(tmp_match$cov, cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
    # need to loop over the appropriate subsets of blocks to assess match
    if(length(block_subset) == 0){
      # we just need balance overall and not additionally on any particular subset
      tmp_balance <- assess_balance(tmp_match$cov, cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
      imbalance_counts[imbalance_counts[,"K"] == K, "Pre-match"] <- sum(abs(tmp_balance$full[,"Std.Diff.Before"]) >= 0.2, na.rm = T)
      imbalance_counts[imbalance_counts[,"K"] == K, "Post-match"] <- sum(abs(tmp_balance$full[,"Std.Diff.After"]) >= 0.2, na.rm = T)
      if (!is.null(tmp_balance$missing)){
        imbalance_counts[imbalance_counts[,"K"] == K, "Pre-match.missing"] <- sum(abs(tmp_balance$missing[,"Std.Diff.Before"] >= 0.2), na.rm = TRUE)
        imbalance_counts[imbalance_counts[,"K"] == K, "Post-match.missing"] <- sum(abs(tmp_balance$missing[,"Std.Diff.After"] >= 0.2), na.rm = TRUE)
      }
      imbalance_counts[imbalance_counts[,"K"] == K, "n"] <- sum(!is.na(tmp_match$cov[,"stratum"]))
    } else{
      pre_match_count <- 0
      post_match_count <- 0
      pre_match_missing_count <- 0
      post_match_missing_count <- 0
      for(bs_ix in 1:length(block_subset)){
        #print(paste("bs_ix = ", bs_ix))
        tmp_balance <- assess_balance(tmp_match$cov[tmp_match$cov[,"block"] %in% block_subset[[bs_ix]],], cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
        pre_match_count <- pre_match_count + sum(abs(tmp_balance$full[,"Std.Diff.Before"]) >= 0.2, na.rm = T)
        post_match_count <- post_match_count + sum(abs(tmp_balance$full[,"Std.Diff.After"]) >= 0.2, na.rm = T)
        if (!is.null(tmp_balance$missing)){
          pre_match_missing_count <- pre_match_missing_count + sum(abs(tmp_balance$missing[,"Std.Diff.Before"]) >= 0.2, na.rm = TRUE)
          post_match_missing_count <- post_match_missing_count + sum(abs(tmp_balance$missing[,"Std.Diff.After"]) >= 0.2, na.rm = TRUE)
        }
      }
      imbalance_counts[imbalance_counts[,"K"] == K, "Pre-match"] <- pre_match_count
      imbalance_counts[imbalance_counts[,"K"] == K, "Post-match"] <- post_match_count
      if (!is.null(tmp_balance$missing)){
        imbalance_counts[imbalance_counts[,"K"] == K, "Pre-match.missing"] <- pre_match_missing_count
        imbalance_counts[imbalance_counts[,"K"] == K, "Post-match.missing"] <- post_match_missing_count
      }
      imbalance_counts[imbalance_counts[,"K"] == K, "n"] <- sum(!is.na(tmp_match$cov[,"stratum"]))
    }
    if(verbose) {
      print(imbalance_counts)
    }
  } # closes loop over K

  match_index <- which(imbalance_counts[,"Post-match"] == 0)
  if(length(match_index) != 0){
    vr_n <- imbalance_counts[match_index,"n"]
    K <- imbalance_counts[which.max(vr_n), "K"]
    n_vr <- max(vr_n)
    print(paste("Best vr match uses", n_vr, "subjects and K =", K))
  } else{
    print("None of the variable ratio matches were adequately balanced on all covariates")
    K <- NULL
    match_index <- which(imbalance_counts[,"Post-match"] < 2)
    if(length(match_index) != 0){
      vr_n <- imbalance_counts[match_index, "n"]
      K <- match_index[which.max(vr_n)]
      n_vr <- max(vr_n)
      print(paste("Will attempt a vr matching using", n_vr, "subjects and K = ", K))
    } else{
      print(paste("All attempted vr matches leave at least", min(imbalance_counts[,"Post-match"], na.rm = TRUE), "covariates imbalanced"))
      K <- NULL
    }
  }
  if(!is.null(K)){
    matched_sets_all <- c()
    exclude_cs_all <- c()
    for(block in unique(prop_score_df[,"block"])){
      #print(paste("Starting block", block, "at", Sys.time()))
      # subset to only those members of the block
      tmp_prop_score_df <- prop_score_df[prop_score_df[,"block"] == block,]
      z <- tmp_prop_score_df[,"treated"]
      prop_score <- tmp_prop_score_df[,"Propensity"]
      logit_prop <- tmp_prop_score_df[,"Logit(Propensity)"]
      names(z) <- rownames(tmp_prop_score_df)
      names(prop_score) <- rownames(tmp_prop_score_df)
      names(logit_prop) <- rownames(tmp_prop_score_df)
      # compute exclude cs
      exclude_cs_treated <- names(z)[z == 1 & prop_score > max(prop_score[z == 0])]
      exclude_cs_control <- names(z)[z == 0 & prop_score < min(prop_score[z == 1])]
      exclude_cs <- c(exclude_cs_treated, exclude_cs_control)
      exclude_cs_all <- c(exclude_cs_all, exclude_cs)
      # pull out the appropriate submatrix of X_mat
      X_mat <- X_mat_full[rownames(tmp_prop_score_df),]
      X_mat <- X_mat[!rownames(X_mat) %in% exclude_cs,]
      z <- z[!names(z) %in% exclude_cs]
      prop_score <- prop_score[!names(prop_score) %in% exclude_cs]
      logit_prop <- logit_prop[!names(logit_prop) %in% exclude_cs]
      
      subject_index <- seq(1, length(z))
      
      drop_cols <- c()
      for(x in colnames(X_mat)){
        if(var(X_mat[,x]) == 0){
          drop_cols <- c(drop_cols,x)
        }
      }
      X_mat <- X_mat[,!colnames(X_mat) %in% drop_cols]
      tmp_dist_mat <- smahal(z, X_mat)
      dist_mat <- addcaliper(tmp_dist_mat, z, logit_prop) # adds  propensity score caliper
      rownames(dist_mat)<- subject_index[z == 1]
      subject_id <- names(z)
      
      # do all the matching
      match_vec_list <- list()
      for(kk in 1:K){
        if(kk == 1){
          treated_ind <- which(z == 1 & prop_score > 1/3)
          control_ind <- which(z == 0 & prop_score > 1/3)
        } else if(kk == K){
          treated_ind <- which(z == 1 & prop_score < 1/(K+1))
          control_ind <- which(z == 0 & prop_score < 1/(K+1))
        } else{
          treated_ind <- which(z == 1 & prop_score > 1/(kk+2) & prop_score < 1/(kk+1))
          control_ind <- which(z == 0 & prop_score > 1/(kk+2) & prop_score < 1/(kk+1))
        }
        row_ind <- rownames(dist_mat) %in% treated_ind
        col_ind <- colnames(dist_mat) %in% control_ind
        n_treated <- sum(row_ind)
        n_control <- sum(col_ind)
        
        #print(paste("kk = ", kk, "n.treated =", n.treated, "n.control =", n.control))
        if(n_control == 1 & n_treated == 1){
          # matching software doesn't handle this well so let's ignore for now
        } else if(n_control >= n_treated & n_treated > 0){
          num_controls <- min(floor(n_control/n_treated), kk)
          tmp_dist_mat <- matrix(dist_mat[row_ind, col_ind], nrow = n_treated, ncol = n_control, dimnames = list(rownames(dist_mat)[row_ind], colnames(dist_mat)[col_ind]))
          match_vec_list[[paste0("kk.",kk)]] <- pairmatch(tmp_dist_mat, controls = num_controls)
        } else if(n_control < n_treated & n_control > 0){ # relabel controls as treated and optimally exclude treated units
          tmp_dist_mat <- t(matrix(dist_mat[row_ind, col_ind], nrow = n_treated, ncol = n_control, dimnames = list(rownames(dist_mat)[row_ind], colnames(dist_mat)[col_ind])))
          match_vec_list[[paste0("kk.",kk)]] <- pairmatch(tmp_dist_mat)
        } else if(n_treated == 0 | n_control == 0){
          #print("No controls or no treated. No match possible in this stratum")
        }
      } # closes loop over kk
      tmp_matched_sets <- stratified.matched.set(match_vec_list, z)
      tmp_matched_sets[exclude_cs] <- NA
      tmp_matched_sets[!is.na(tmp_matched_sets)] <- paste0(block, ".", tmp_matched_sets[!is.na(tmp_matched_sets)])
      matched_sets_all <- c(matched_sets_all, tmp_matched_sets)
    } # closes loop over the blocks
    
    # make the data frame
    tmp_match <- stratified.df(prop_score_df, out_data, matched_sets_all, z_all)
    # assess balance
    if(length(block_subset) == 0){
      tmp_balance <- assess_balance(tmp_match$cov, cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
      return(list("cov" = tmp_match$cov, "out" = tmp_match$out, "balance" = tmp_balance, "imbalance_counts" = imbalance_counts, "exclude_cs_all" = exclude_cs_all))
    } else{
      balance_list <- list()
      for(bs_ix in 1:length(block_subset)){
        balance_list[[bs_ix]] <- assess_balance(tmp_match$cov[tmp_match$cov[,"block"] %in% block_subset[[bs_ix]],], cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
      }
      return(list("cov" = tmp_match$cov, "out" = tmp_match$out, "balance" = balance_list, "imbalance_counts" = imbalance_counts, "exclude_cs_all" = exclude_cs_all))

    }
  } else{
    return(list("imbalance_counts" = imbalance_counts))
  }
}


