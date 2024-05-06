library(tidyverse)
library(sensitivityfull)



df_test <- data.frame(
    IDS = 1:12,
    y = c(
        0, 0, 1, 0,
        1, 1, 1, 0,
        1, 0, 0, 0),
    z = c(
        0, 0, 1, 0,
        1, 1, 1, 0,
        1, 0, 0, 0),
    strat = rep(c(1, 2, 3), each = 4)
)


df_tmp <- df_test
strat_colname <- "strat"
trt_colname <- "z"
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

ymat
