# Testing the effect on continuous outcomes of a hierarchy of exposures...
# using the test of composite null for difference...
# by Fogharty et al.

source("https://raw.githubusercontent.com/colinbfogarty/SensitivityCompositeBinary/main/compositeBinary.R") # nolint: line_length_linter.

# Return value:
# - p_vals: list of p-values for the individual tests
# - decisions: was each hyp was accepted(ACC), rejected(REJ) or not tested(DNT)
#NOTE the term ``accepted'' is only used for simplicty, it means fail to reject

#TODO: Needs an IDS column!!!!!!

hier_binary_tests <- function(
    dt_outcome_trt,  # dataframe with observed outcomes, treatment, and stratum
    treat_names,     # treatment column names
    strat_names,     # named list mapping trt colname -> stratum colname
    outcome_name,    # name of the outcome column
    treat_hierarchy, # named list mapping trt colname -> parent trt colname
    p_budget,        # named list mapping trt colname -> type I budget,
    fwer,            # Family Wise Error Rate (for uncorrected intervals)
    gamma            # sensitivity parameter
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

    unc_tau_est <- rep(NA, n_treat)
    names(unc_tau_est) <- names(treat_hierarchy)

    cor_tau_est <- rep(NA, n_treat)
    names(cor_tau_est) <- names(treat_hierarchy)

    unc_confints = list()
    cor_confints = list()
    unc_conf_objects = list()
    cor_conf_objects = list()
    unc_conf_levels = list()
    cor_conf_levels = list()

    for (tvar in treat_names) {
        print(paste0("Binary test for ", tvar))
        strat_colname <- strat_names[[tvar]]
        trt_colname <-  tvar
        trt_parent <- treat_hierarchy[[trt_colname]]

            # genereate subset for matched samples
            df_tmp <- dt_outcome_trt %>%
                filter(!is.na(get(strat_colname)))

            # run test
            unc_comp_test <- compositeBinary( # test of coposite null
                df_tmp[[strat_colname]],
                df_tmp[[outcome_name]],
                df_tmp[[trt_colname]],
                alpha = fwer,
                Gamma = gamma
            )

        unc_tau_est[[trt_colname]] <- unc_comp_test$estimate
        p_vals[[trt_colname]] <- unc_comp_test$pval
        unc_confints[[trt_colname]] <- unc_comp_test$confint
        unc_conf_levels[[trt_colname]] <- unc_comp_test$confidencelevel
        unc_conf_objects[[trt_colname]] <- unc_comp_test

        cor_comp_test <- compositeBinary( # test of coposite null
                df_tmp[[strat_colname]],
                df_tmp[[outcome_name]],
                df_tmp[[trt_colname]],
                alpha = p_budget[[trt_colname]],
                Gamma = gamma
        )

        cor_tau_est[[trt_colname]] <- cor_comp_test$estimate
        cor_confints[[trt_colname]] <- cor_comp_test$confint
        cor_conf_levels[[trt_colname]] <- cor_comp_test$confidencelevel
        cor_conf_objects[[trt_colname]] <- cor_comp_test


        # run test only if root or parent rejected
        if (trt_parent == "root" | rej_v[trt_parent] == "REJ") {
            rej_v[[trt_colname]] <- ifelse(
                unc_comp_test$pval < p_budget[[trt_colname]],
                "REJ", "ACC"
            )
        }
    }

    return(
        list(
            p_values = p_vals,
            decisions = rej_v,
            unc_tau_est = unc_tau_est,
            unc_confints = unc_confints,
            unc_confint_objects = unc_conf_objects,
            unc_conf_levels = unc_conf_levels,
            cor_tau_est = cor_tau_est,
            cor_confints = cor_confints,
            cor_confint_objects = cor_conf_objects,
            cor_conf_levels = cor_conf_levels
        )
    )
}