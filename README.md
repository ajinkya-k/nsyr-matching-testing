# Matching and Testing for Tree structured Hypotheses

Repository with code for the matching, testing, and simulation studies for the NSYR sports participation study.
Pre-analysis protocol is available on arXiv: [arXiv:2211.02104](https://arxiv.org/abs/2211.02104).


# Usage

To use the code given in this repository, three things are required:

1. A dataframe with the covariates and the treatments.
2. A list of the names of the treatment variables
3. A dictionary (i.e. named list) with the keys being the names of the treatments and the values being the names of the parent treatment. For the root node treatment, the value must be "root"

The following code chunk shows an example of the treatment names and the associated hierarchy, for the exposure tree depicted in Figure 1 of the pre-analysis protocol[^1].

```r
treat_vars <- c(
    "tr_any_activity",
    "tr_any_sports",
    "tr_no_sports",
    "tr_any_contact",
    "tr_no_contact",
    "tr_any_collision",
    "tr_no_collision"
)

trt_hierarcy <- list(
    "tr_any_activity" = "root",
    "tr_any_sports" = "tr_any_activity",
    "tr_no_sports" = "tr_any_activity",
    "tr_any_contact" = "tr_any_sports",
    "tr_no_contact" = "tr_any_sports",
    "tr_any_collision" = "tr_any_contact",
    "tr_no_collision" = "tr_any_contact"
)
```
Once all of these are defined run using the following template: 
```r
source("matching/matching_wrapper_functions.R")

full_matches <- full_match_for_tree(
    dt_covars_treat, # dataframe with IDS, the covariates, treatments, and prop scores for treatments
    dt_outcome, # dataframe containing outcomes
    matched_covars, # names of pre-exposure covariates to be used for matching
    cont_covars, # names of continuous pre-exposure covariates
    cat_covars, # names of categorical pre-exposure covariates
    treat_vars, # names of the treatment columns in dt_covars_treat
    treat_hierachy, # treatment hierarchy named list (as described in README)
    out_name, # name of the outcome (i.e. response) variable
    keep_ids, # IDS of the rows to keep in the dataset,
    excl_rows, # rows of the dataset to exclude (functionality not checked)
    max_K = NULL # max K to use in full matching with control-to-exposure is between 1:K and K:1
)
```

[^1]: Ajinkya H. Kokandakar, Yuzhou Lin, Steven Jin, Jordan Weiss, Amanda R. Rabinowitz, Reuben A. Buford May, Sameer K. Deshpande, Dylan Small. _Protocol for an observational study on the effects of adolescent sports participation on health in early adulthood_. [arXiv:2211.02104](https://arxiv.org/abs/2211.02104)
