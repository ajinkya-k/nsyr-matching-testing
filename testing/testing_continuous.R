# Testing the effect on continuous outcomes of a hierarchy of exposures
# using the `senfm` package by Rosenbaum et al.

hier_continuous_test <- function(
    dt_outcome_trt,  # dataframe with observed outcomes, treatment, and stratum
    treat_names,     # treatment column names
    strat_names,     # named list mapping trt colname -> stratum colname
    treat_hierarchy, # named list mapping trt colname -> parent trt colname
    p_budget         # named list mapping trt colname -> type I budget
) {

    n_treat <- length(treat_hierarchy)

    # create named empty vector for p-values
    p_vals <- rep(NA, n_treat)
    names(p_vals) <- names(treat_hierarchy)

    # create named empty vector for the decision for each hyp
    # - DNT: null hypothesis not tested (default)
    # - ACC: null hypothesis tested and NOT rejected
    # - REJ: null hypothesis tested and WAS rejected
    rej_v <- rep("DNT", n_treat)
    names(rej_v) <- names(treat_hierarchy)

    for (tvar in treat_names) {
        strat_colname <- strat_names[[tvar]]
        trt_colname <-  tvar
        trt_parent <- treat_hierarchy[[trt_colname]]

        if (is.na(trt_parent) || rej_v[trt_parent] == "REJ") {
            df_tmp <- dt_outcome_trt %>%
                filter(!is.na(get(strat_colname)))
            rownames(df_tmp) <- df_tmp$IDS

            # find max size of matched set
            max_sz <- df_tmp %>%
                count(get(strat_colname)) %>%
                pull(n) %>%
                max()

            # create the y matrix needed for `senfm` for one treatment
            ymat <- matrix(NA, length(unique(df_tmp[[strat_colname]])), max_sz)
            # create the indicator vector for senfm
            first_trt <- rep(NA, length(unique(df_tmp[[strat_colname]])))

            # fill in the ymat
            selvarlist <- c("IDS", trt_colname)
            indexlst <- list()
            first_trt <- rep(NA, length(unique(df_tmp[[strat_colname]])))
            for (strat in unique(df_tmp[[strat_colname]])) {
                tmp_trt_y <- df_tmp %>%
                    filter(get(strat_colname) == !!strat) %>%
                    select(all_of(selvarlist)) %>%
                    group_by(get(trt_colname)) %>%
                    mutate(n_per_trt = n()) %>%
                    ungroup() %>%
                    arrange(n_per_trt)
                first_trt[strat] <- tmp_trt_y[[1, trt_colname]]
                indexlst[[strat]] <- tmp_trt_y %>% pull(IDS)
                yv <- df_tmp[as.character(indexlst[[strat]]), "y"]
                ymat[strat, 1:length(yv)] <- yv
            }
            # m-test for one side
            sentest_gt <- senfm(
                ymat, as.logical(first_trt), # data
                gamma = 1, inner = 0, trim = 3, lambda = 1/2, #parameters
                tau = 0, # null value
                alternative = "greater" # type of test
            )

            # m-test for the other side
            sentest_lt <- senfm(
                ymat, as.logical(first_trt), # data
                gamma = 1, inner = 0, trim = 3, lambda = 1/2, # parameters
                tau = 0, # null
                alternative = "less" # type of test
            )

            pval_mtest <- min(1, 2 * min(sentest_gt$pval, sentest_lt$pval))

            p_vals[[trt_colname]] <- pval_mtest
            rej_v[[trt_colname]] <- ifelse(
                pval_mtest < p_budget[[trt_colname]],
                "REJ", "ACC"
            )
        }

    }
    return(
        list(p_values = p_vals, decisions = rej_v)
    )
}
