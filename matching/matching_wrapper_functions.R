
# 1. dt_covars_treat: dataframe with IDS, the covariates, treatments, and prop scores for treatments
# for each trt say tr_name, the propensity score column should be named prop_tr_name
# dt_outcome: dataframe containing IDS and outomces
# matched_covars: Names of covariate columns to be used for matching
# tr_name: name of the treatment to be used for matching
# out_name: name of the outcome variable in dt_outcome
# incl_rows: rows to explicitly include
# excl_rows: rows to explicilty exclude
fit_vr_match <- function(
    dt_covars_treat,
    dt_outcome,
    matched_covars,
    cont_covars,
    cat_covars,
    tr_name,
    out_name,
    incl_rows,
    excl_rows
) {
    dt_xt <- dt_covars_treat[incl_rows, ]
    ps_val <- dt_xt[[paste0("prop_", tr_name)]]
    ps_df <- data.frame(
            "Propensity" = ps_val,
            "Logit(Propensity)" = log(ps_val / (1 - ps_val)),
            "treated" = dt_xt[[tr_name]],
            "block" = rep(1, length(ps_val)), # assign all to block 1
            row.names = rownames(dt_xt)
        )

    setnames(ps_df, "Logit.Propensity.", "Logit(Propensity)")

    ps_df <- cbind(
        ps_df,
        dt_xt[, matched_covars]
    )

    vr_fit <- best_vr_match_block(
        ps_df,
        dt_outcome[incl_rows, out_name],
        ps_df[, matched_covars],
        cont_vars = cont_covars,
        cat_vars = cat_covars
    )

    return(vr_fit)
}


# 1. dt_covars_treat: dataframe with IDS, the covariates, treatments, and prop scores for treatments
# for each trt say tr_name, the propensity score column should be named prop_tr_name
# 2. dt_outcome: dataframe containing IDS and outomces
# 3. matched_covars: Names of covariate columns to be used for matching
# 4. tr_name: name of the treatment to be used for matching
# 5. out_name: name of the outcome variable in dt_outcome
# 6. incl_rows: rows to explicitly include
# 7. excl_rows: rows to explicilty exclude
# 8. max_K: max ratio of control:trt or vice versa, i.e. ratio is between 1:K and K:1 at most
fit_full_match <- function(
    dt_covars_treat,
    dt_outcome,
    matched_covars,
    cont_covars,
    cat_covars,
    tr_name,
    out_name,
    incl_rows,
    excl_rows,
    max_K = 10
) {
    dt_xt <- dt_covars_treat[incl_rows, ]
    ps_val <- dt_xt[[paste0("prop_", tr_name)]]
    ps_df <- data.frame(
            "Propensity" = ps_val,
            "Logit(Propensity)" = log(ps_val / (1 - ps_val)),
            "treated" = dt_xt[[tr_name]],
            "block" = rep(1, length(ps_val)), # assign all to block 1
            row.names = rownames(dt_xt)
        )

    setnames(ps_df, "Logit.Propensity.", "Logit(Propensity)")

    ps_df <- cbind(
        ps_df,
        dt_xt[, matched_covars]
    )

    vr_fit <- best_full_match_block(
        ps_df,
        dt_outcome[incl_rows, out_name],
        ps_df[, matched_covars],
        cont_vars = cont_covars, # c("age"),
        cat_vars = cat_covars, # matched_covars[-1],
        max_K = max_K
    )

    return(vr_fit)
}

# Makes the Love plot with arrows for ONE fit object obtained from `best_vr_match_block`
love_plot_arrow <- function(vr_fit) {
    baltabl <- vr_fit$balance$full

    baltabl <- baltabl[which(baltabl$X.val != 0) , ]
    baltabl$Diff.Before <- with(baltabl, abs(Before.Treated.Mean - Before.Control.Mean))
    baltabl$Diff.After <- with(baltabl, abs(After.Treated.Mean - After.Control.Mean))
    bal_love <- baltabl[,
            c(
                "X",
                "Std.Diff.Before",
                "Std.Diff.After",
                "Diff.Before",
                "Diff.After"
            )
        ]

    bal_love$"Std.Diff.Before" <- abs(bal_love$"Std.Diff.Before")
    bal_love$"Std.Diff.After" <- abs(bal_love$"Std.Diff.After")

    bal_love <- bal_love[order(bal_love$"Std.Diff.Before", decreasing = TRUE), ]
    bal_love$X <-  as.factor(bal_love$X)
    bal_love$sgn <- as.factor(with(bal_love, ifelse(Std.Diff.After < Std.Diff.Before, "Good", "Bad")))

    # plot(x = as.factor(bal_love$X), y = bal_love$"Std.Diff.Before")


    p1 <- bal_love %>%
        mutate(X = fct_reorder(X, Std.Diff.Before)) %>%
        ggplot() +
        geom_point(aes(x=X, y = Std.Diff.Before, colour="before")) +
        geom_line(aes(x=X, y = Std.Diff.Before, colour="before", group=1)) + 
        geom_point(aes(x=X, y = Std.Diff.After, colour = "after")) +
        geom_line(aes(x=X, y = Std.Diff.After, colour="after", group=1)) + 
        geom_hline(yintercept=0.1) +
        coord_flip()

    #
    p2 <- bal_love %>%
        mutate(X = fct_reorder(X, Std.Diff.Before)) %>%
        ggplot() +
        geom_point(aes(x=X, y = Std.Diff.Before, colour = sgn)) +
        geom_point(aes(x=X, y = Std.Diff.After, colour = sgn)) +
        geom_segment(
            aes(x=X, xend=X, y = Std.Diff.Before, yend=Std.Diff.After, colour = sgn),
            arrow = arrow(length = unit(0.25, "cm"))
            ) +
        geom_hline(yintercept=c(0.1, 0.2), color="gray") +
        labs(x = "covariates", y = "standardized difference") +
        scale_x_discrete(labels = lbls) +
        scale_colour_manual(
          values = c("red","black"),
          aesthetics = c("colour")
        ) +
        theme_bw() +
        theme(legend.position="none") +
        coord_flip()
    #
    return(list(p1 = p1, p2 = p2))
}



# 1. dt_covars_treat: dataframe with IDS, the covariates, treatments, and prop scores for treatments
# dt_outcome: dataframe containing outcomes
# matched_covars: names of pre-exposure covariates to be used for matching
# cont_covars: names of continuous covarites,
# cat_covars: names of categorical covariates,
# treat_vars: names of the treatment columns,
# treat_hierachy: named list with key = treatment name and value = parent treatment name, (value="root" for root treatment),
# out_name: name of the outcome variable,
# keep_ids: IDS of the rows to keep in the dataset,
# excl_rows: row nums to explicity exclude,
# max_K = NULL


full_match_for_tree <- function(
    dt_covars_treat,
    dt_outcome,
    matched_covars,
    cont_covars,
    cat_covars,
    treat_vars,
    treat_hierachy,
    out_name,
    keep_ids,
    excl_rows,
    max_K = NULL
) {
    fullmatch_objects <- list()
    # top level treatment
    trt_name <- treat_vars[1]

    # if max K not set, use max_K = 10
    if (is.null(max_K)) {
        print("max K not provided, using max_K = 10")
        max_K = 10
    }
    vr_fit <- fit_full_match(
        dt_covars_treat,
        dt_outcome,
        matched_covars,
        cont_covars, # c("age"),
        cat_covars, # matched_covars[-1],
        trt_name,
        out_name,
        keep_ids,
        excl_rows = NULL, max_K = max_K)

    fullmatch_objects[[trt_name]] <- vr_fit

    # matching for the non-root treatments
    for (tr_name in treat_vars[-1]) {
        print(
            paste0("Computing matches for ", tr_name)
        )
        tr_parent <- treat_hierachy[[tr_name]]
        print(
            paste0("Parent treatment: ", tr_parent)
        )
        vr_fit <- fullmatch_objects[[tr_parent]]
        keep2 <- rownames(vr_fit$cov[which(!is.na(vr_fit$cov$stratum)), ])

        print(
            paste0("For ", tr_name, " we start with sample size: ", length(keep2) )
        )
        dt_xt_only_sports <- dt_covars_treat[keep2, ]
        keep2 <- rownames(dt_xt_only_sports[which(!is.na(dt_xt_only_sports[, tr_name])) , ])

        vr_fit2 <- fit_full_match(
            dt_covars_treat,
            dt_outcome,
            matched_covars,
            cont_covars,
            cat_covars,
            tr_name,
            out_name,
            as.character(keep2),
            excl_rows = NULL)
        fullmatch_objects[[tr_name]] <- vr_fit2
        print(paste0(rep("=", 78), collapse = ""))
    }

    return(fullmatch_objects)
}
