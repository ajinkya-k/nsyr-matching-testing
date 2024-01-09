# Testing the effect on continuous outcomes of a hierarchy of exposures...
# using the test of composite null for difference...
# by Fogharty et al.

source("https://raw.githubusercontent.com/colinbfogarty/SensitivityCompositeBinary/main/compositeBinary.R") # nolint: line_length_linter.

# Return value:
# - p_vals: list of p-values for the individual tests
# - decisions: was each hyp was accepted(ACC), rejected(REJ) or not tested(DNT)
#NOTE the term ``accepted'' is only used for simplicty, it means fail to reject

hier_binary_tests <- function(
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
            # genereate subset for matched samples
            df_tmp <- dt_outcome_trt %>%
                filter(!is.na(get(strat_colname)))

            # run test
            comp_test <- compositeBinary( # test of coposite null
                df_tmp[[strat_colname]],
                df_tmp$y,
                df_tmp[[trt_colname]]
            )

            p_vals[[trt_colname]] <- comp_test$pval
            rej_v[[trt_colname]] <- ifelse(
                comp_test$pval < p_budget[[trt_colname]],
                "REJ", "ACC"
            )
        }
    }

    return(
        list(p_values = p_vals, decisions = rej_v)
    )
}