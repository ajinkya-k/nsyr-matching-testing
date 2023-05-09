  # best_block_vr_match

# Loop over the blocks getting the matched_set_vectors
# pull out the vector indicating matched_set id. this is currently saved as tmp_matched_set
# paste in the block number to these elements


best_full_match_block <- function(prop_score_df, out_data, X_mat_full, max_K = 5, block_subset = list(),cont_vars = c(), cat_vars = c(), mix_vars = c(), mix_values = list(), verbose = TRUE){
  # make a matrix to count the imbalances
  imbalance_counts <- matrix(nrow = max_K, ncol = 7, dimnames = list(c(), c("K","Pre-match", "Post-match", "Pre-match.missing", "Post-match.missing", "n", "num.weak.bal")))
  imbalance_counts[,"K"] <- 1:max_K
  z_all <- prop_score_df$treated
  names(z_all) <- rownames(prop_score_df)
  for(K in 1:max_K){
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
      
      fm <- fullmatch(dist_mat, min.controls = 1/K, max.controls = K)
      fm <- fm[!is.na(fm)]
      # do all the matching
      tmp_matched_sets <- matched.set(fm, z)
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
      wb <- abs(tmp_balance$full[,"Std.Diff.After"]) < 0.2 &  abs(tmp_balance$full[,"Std.Diff.After"]) > 0.1 # wb means weakly balanced
      imbalance_counts[imbalance_counts[,"K"] == K, "num.weak.bal"] <- sum(wb, na.rm = T)
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
      wb <- abs(tmp_balance$full[,"Std.Diff.After"]) < 0.2 &  abs(tmp_balance$full[,"Std.Diff.After"]) > 0.1 # wb means weakly balanced
      imbalance_counts[imbalance_counts[,"K"] == K, "num.weak.bal"] <- sum(wb, na.rm = T)
      if (!is.null(tmp_balance$missing)){
        imbalance_counts[imbalance_counts[,"K"] == K, "Pre-match.missing"] <- pre_match_missing_count
        imbalance_counts[imbalance_counts[,"K"] == K, "Post-match.missing"] <- post_match_missing_count
      }
      imbalance_counts[imbalance_counts[,"K"] == K, "n"] <- sum(!is.na(tmp_match$cov[,"stratum"]))
    }
  } # closes loop over K

  if(verbose) {
      print(imbalance_counts)
  }
  match_index <- which(imbalance_counts[,"Post-match"] == 0 & imbalance_counts[,"n"] > 0)
  if(length(match_index) != 0){
    # print("Zero imbalances: ")
    # print(match_index)
    vr_n <- imbalance_counts[match_index,"num.weak.bal"]
    # print(vr_n)
    K <- imbalance_counts[match_index[which.min(vr_n)], "K"]
    n_vr <- min(vr_n)
    print(paste("Best vr match uses", n_vr, "subjects and K =", K))
  } else{
    print("None of the variable ratio matches were adequately balanced on all covariates")
    K <- NULL
    match_index <- which(imbalance_counts[,"Post-match"] < 2)
    if(length(match_index) != 0){
      vr_n <- imbalance_conts[match_index, "n"]
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
      fm <- fullmatch(dist_mat, min.controls = 1/K, max.controls = K) # closes loop over kk
      fm <- fm[!is.na(fm)]
      tmp_matched_sets <- matched.set(fm, z)
      tmp_matched_sets[exclude_cs] <- NA
      tmp_matched_sets[!is.na(tmp_matched_sets)] <- paste0(block, ".", tmp_matched_sets[!is.na(tmp_matched_sets)])
      matched_sets_all <- c(matched_sets_all, tmp_matched_sets)
    } # closes loop over the blocks
    
    # make the data frame
    tmp_match <- stratified.df(prop_score_df, out_data, matched_sets_all, z_all)
    # assess balance
    if(length(block_subset) == 0){
      tmp_balance <- assess_balance(tmp_match$cov, cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
      return(list("cov" = tmp_match$cov, "out" = tmp_match$out, "balance" = tmp_balance, "imbalance_counts" = imbalance_counts, "exclude_cs_all" = exclude_cs_all, "exclude_cs_treated" = exclude_cs_treated, "exclude_cs_control" = exclude_cs_control))
    } else{
      balance_list <- list()
      for(bs_ix in 1:length(block_subset)){
        balance_list[[bs_ix]] <- assess_balance(tmp_match$cov[tmp_match$cov[,"block"] %in% block_subset[[bs_ix]],], cont.vars = cont_vars, cat.vars = cat_vars, mix.vars = mix_vars, mix.values = mix_values)
      }
      return(list("cov" = tmp_match$cov, "out" = tmp_match$out, "balance" = balance_list, "imbalance_counts" = imbalance_counts, "exclude_cs_all" = exclude_cs_all, "exclude_cs_treated" = exclude_cs_treated, "exclude_cs_control" = exclude_cs_control))

    }
  } else{
    return(list("imbalance_counts" = imbalance_counts))
  }
}


