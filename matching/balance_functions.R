# Functions used to assess balance

# assess_balance_cont: for continuous measurements
# assess_balance_cat: for categorical measurements
# assess_balance_mix: for "mixed" measurements. Some values are categorical and others are continuous
library(Hmisc)
assess_balance_mix <- function(x,cov.data, cat.values = 0){
  
  # Get the covariates
  # Check whether paste0(x, ".missing")
  # if it is and there are missing values, then we need 2 + length(cat.values) rows in the match.table
  if(x %in% colnames(cov.data)){

    if(paste0(x, ".missing") %in% colnames(cov.data)){
      match.table <- data.frame("X" = rep(x, times = 2 + length(cat.values)), 
  	                            "X.val" = c(sort(cat.values), "missing", "continuous"),
  	                            "Before.Treated.Mean" = rep(NA, times = 2 + length(cat.values)),
  	                            "After.Treated.Mean" = rep(NA, times = 2 + length(cat.values)),
  	                            "Before.Control.Mean" = rep(NA, times = 2 + length(cat.values)),
  	                            "After.Control.Mean" = rep(NA, times = 2 + length(cat.values)),
  	                            "Std.Diff.Before" = rep(NA, times = 2 + length(cat.values)),
  	                            "Std.Diff.After" = rep(NA, times = 2 + length(cat.values)),
  	                            "Before.Treated.Count" = rep(NA, times = 2 + length(cat.values)),
  	                            "Before.Control.Count" = rep(NA, times = 2 + length(cat.values)),
  	                            "After.Treated.Count" = rep(NA, times = 2 + length(cat.values)),
  	                            "After.Control.Count" = rep(NA, times = 2 + length(cat.values)))
  	 
  	  # Deal with the categorical variables
  	  before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	  after.control.index <- which(cov.data[,"treated"] == 0 & !is.na(cov.data[,"stratum"]))

  	 
      for(x.val in cat.values){
  	    before.treated.mean <- mean(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE)
  	    before.control.mean <- mean(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE)
  	   
  	    after.treated.mean <- wtd.mean(1*(cov.data[after.treated.index,x] == x.val), weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
  	    after.control.mean <- wtd.mean(1*(cov.data[after.control.index,x] == x.val), weight = cov.data[after.control.index, "weight"], na.rm = TRUE)
  	   
  	    combined.sd <- sqrt(0.5*var(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE) + 0.5*var(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE))
  	   
  	    std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	    std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	   
  	    before.treated.count <- paste0(sum(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE), "/", length(before.treated.index))
  	    before.control.count <- paste0(sum(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE), "/", length(before.control.index))
  	    after.treated.count <- paste0(sum(cov.data[after.treated.index,x] == x.val, na.rm = TRUE), "/", length(after.treated.index))
  	    after.control.count <- paste0(sum(cov.data[after.control.index,x] == x.val, na.rm = TRUE), "/", length(after.control.index))
  	   
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Count"] <- before.treated.count
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Count"] <- before.control.count
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Count"] <- after.treated.count
        match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Count"] <- after.control.count
  	  }
  	  # Now deal with missingness
  	  before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	  after.control.index <- which(cov.data["treated"] == 0 & !is.na(cov.data[,"stratum"]))
  	 
  	  before.treated.mean <- mean(cov.data[before.treated.index, paste0(x, ".missing")], na.rm = TRUE)
  	  before.control.mean <- mean(cov.data[before.control.index, paste0(x, ".missing")], na.rm = TRUE)
  	  after.treated.mean <- wtd.mean(cov.data[after.treated.index, paste0(x, ".missing")], weight = cov.data[after.treated.index, "weight"], na.rm = TRUE)
  	  after.control.mean <- wtd.mean(cov.data[after.control.index, paste0(x, ".missing")], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
  	 
  	  combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,paste0(x, ".missing")], na.rm = TRUE) + 0.5*var(cov.data[before.control.index,paste0(x, ".missing")], na.rm = TRUE))
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	 
  	  before.treated.count <- paste0(sum(cov.data[before.treated.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(before.treated.index))
  	  before.control.count <- paste0(sum(cov.data[before.control.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(before.control.index))
  	  after.treated.count <- paste0(sum(cov.data[after.treated.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(after.treated.index))
  	  after.control.count <- paste0(sum(cov.data[after.control.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(after.control.index))
 
  	  match.table[match.table[,"X.val"] == "missing","Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "missing","Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "missing","Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "missing","Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "missing","After.Treated.Count"] <- after.treated.count
  	  match.table[match.table[,"X.val"] == "missing","After.Control.Count"] <- after.control.count

  	  # Now deal with the continuous
  	  # we want to take averages only of the treated and controls that have continuous values and are not missing
  	  before.treated.index <- which(cov.data[,"treated"] == 1 & !cov.data[,x] %in% cat.values & cov.data[,paste0(x, ".missing")] == 0)
  	  before.control.index <- which(cov.data[,"treated"] == 0 & !cov.data[,x] %in% cat.values & cov.data[,paste0(x, ".missing")] == 0)
  	 
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !cov.data[,x] %in% cat.values & cov.data[,paste0(x, ".missing")] == 0 & !is.na(cov.data[,"stratum"]))
      after.control.index <- which(cov.data[,"treated"] == 0 & !cov.data[,x] %in% cat.values & cov.data[,paste0(x, ".missing")] == 0 & !is.na(cov.data[,"stratum"]))
      
      before.treated.mean <- mean(cov.data[before.treated.index,x], na.rm = TRUE)
      before.control.mean <- mean(cov.data[before.control.index,x], na.rm = TRUE)
      
      after.treated.mean <- wtd.mean(cov.data[after.treated.index,x], weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
      after.control.mean <- wtd.mean(cov.data[after.control.index,x], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
      
      combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,x], na.rm = TRUE) + 0.5*var(cov.data[before.control.index,x], na.rm = TRUE))
  	  
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	  
  	  before.treated.count <- NA
  	  before.control.count <- NA
  	  after.treated.count <- NA
  	  after.control.count <- NA
  	  
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Mean"] <- round(before.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Mean"] <- round(before.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Mean"] <- round(after.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Mean"] <- round(after.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Count"] <- after.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Count"] <- after.control.count
    } else{   # if not, match.table table requires 1 + length(cat.values) rows in the match.table

  	  match.table <- data.frame("X" = rep(x, times = 1 + length(cat.value)),
  	                          "X.val" = c(cat.value, "continuous"),
  	                          "Before.Treated.Mean" = rep(NA, times = 1 + length(cat.value)),
  	                          "After.Treated.Mean" = rep(NA, times = 1 + length(cat.value)),
  	                          "Before.Control.Mean" = rep(NA, times = 1 + length(cat.value)),
  	                          "After.Control.Mean" = rep(NA, times = 1 + length(cat.value)),
  	                          "Std.Diff.Before" = rep(NA, times = 1 + length(cat.value)),
  	                          "Std.Diff.After" = rep(NA, times = 1 + length(cat.value)),
  	                          "Before.Treated.Count" = rep(NA, times = 1 + length(cat.value)),
  	                          "Before.Control.Count" = rep(NA, times = 1 + length(cat.value)),
  	                          "After.Treated.Count" = rep(NA, times = 1 + length(cat.value)),
  	                          "After.Control.Count" = rep(NA, times = 1 + length(cat.value)))
  	
      # Deal with the categorical variables
      before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	  after.control.index <- which(cov.data[,"treated"] == 0 & !is.na(cov.data[,"stratum"]))

   	  for(x.val in cat.values){
  	    before.treated.mean <- mean(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE)
  	    before.control.mean <- mean(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE)
  	   
  	    after.treated.mean <- wtd.mean(1*(cov.data[after.treated.index,x] == x.val), weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
  	    after.control.mean <- wtd.mean(1*(cov.data[after.control.index,x] == x.val), weight = cov.data[after.control.index, "weight"], na.rm = TRUE)
  	   
  	    combined.sd <- sqrt(0.5*var(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE) + 0.5*var(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE))
  	   
  	    std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	    std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	   
  	    before.treated.count <- paste0(sum(cov.data[before.treated.index,x] == x.val, na.rm = TRUE), "/", length(before.treated.index))
  	    before.control.count <- paste0(sum(cov.data[before.control.index,x] == x.val, na.rm = TRUE), "/", length(before.control.index))
  	    after.treated.count <- paste0(sum(cov.data[after.treated.index,x] == x.val, na.rm = TRUE), "/", length(after.treated.index))
  	    after.control.count <- paste0(sum(cov.data[after.control.index,x] == x.val, na.rm = TRUE), "/", length(after.control.index))
  	   
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Count"] <- before.treated.count
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Count"] <- before.control.count
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Count"] <- after.treated.count
        match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Count"] <- after.control.count 	    
      }
  	  # Now deal with the continuous variables. We want to restrict attention to the people without categorical values
  	  before.treated.index <- which(cov.data[,"treated"] == 1 & !cov.data[,x] %in% cat.values)
  	  before.control.index <- which(cov.data[,"treated"] == 0 & !cov.data[,x] %in% cat.values)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !cov.data[,x] %in% cat.values & !is.na(cov.data[,"stratum"]))
      after.control.index <- which(cov.data[,"treated"] == 0 & !cov.data[,x] %in% cat.values & !is.na(cov.data[,"stratum"]))

      before.treated.mean <- mean(cov.data[before.treated.index,x], na.rm = TRUE)
      before.control.mean <- mean(cov.data[before.control.index,x], na.rm = TRUE)
      
      after.treated.mean <- wtd.mean(cov.data[after.treated.index,x], weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
      after.control.mean <- wtd.mean(cov.data[after.control.index,x], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
      
      combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,x], na.rm = TRUE) + 0.5*var(cov.data[before.control.index,x], na.rm = TRUE))
  	  
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	  
  	  before.treated.count <- NA
  	  before.control.count <- NA
  	  after.treated.count <- NA
  	  after.control.count <- NA
  	  
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Mean"] <- round(before.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Mean"] <- round(before.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Mean"] <- round(after.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Mean"] <- round(after.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Count"] <- after.treated.count
      match.table[match.table[,"X.val"] == "continuous", "After.Control.Count"] <- after.control.count
    } # closes else statement seeing whether x.missing is in colnames(cov.data)
    return(match.table)
  } else { # closes if statement checking whether x is in colnames(cov.data)
  	return(NULL)
  }
}

assess_balance_cat <- function(x,cov.data){
  
  # Get the covariates
  # Check whether paste0(x, ".missing")
  # if it is and there are missing values, then we need 2 + length(cat.value) rows in the match.table
  if(x %in% colnames(cov.data)){

    if(paste0(x, ".missing") %in% colnames(cov.data)){
      unik.values <- sort(unique(cov.data[cov.data[,paste0(x, ".missing")] == 0,x])) # unique values that are NOT imputed
      n.values <- length(unik.values) + 1
      match.table <- data.frame("X" = rep(x, times = n.values), 
  	                            "X.val" = c(unik.values, "missing"),
  	                            "Before.Treated.Mean" = rep(NA, times = n.values),
  	                            "After.Treated.Mean" = rep(NA, times = n.values),
  	                            "Before.Control.Mean" = rep(NA, times = n.values),
  	                            "After.Control.Mean" = rep(NA, times = n.values),
  	                            "Std.Diff.Before" = rep(NA, times = n.values),
  	                            "Std.Diff.After" = rep(NA, times = n.values),
  	                            "Before.Treated.Count" = rep(NA, times = n.values),
  	                            "Before.Control.Count" = rep(NA, times = n.values),
  	                            "After.Treated.Count" = rep(NA, times = n.values),
  	                            "After.Control.Count" = rep(NA, times = n.values))
  	 
  	  # Deal with the categorical variables
  	  before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	  after.control.index <- which(cov.data[,"treated"] == 0 & !is.na(cov.data[,"stratum"]))

  	 
      for(x.val in unik.values){
  	    before.treated.mean <- mean(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE)
  	    before.control.mean <- mean(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE)
  	   
  	    after.treated.mean <- wtd.mean(1*(cov.data[after.treated.index,x] == x.val), weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
  	    after.control.mean <- wtd.mean(1*(cov.data[after.control.index,x] == x.val), weight = cov.data[after.control.index, "weight"], na.rm = TRUE)
  	   
  	    combined.sd <- sqrt(0.5*var(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE) + 0.5*var(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE))
  	   
  	    std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	    std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	   
  	    before.treated.count <- paste0(sum(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE), "/", length(before.treated.index))
  	    before.control.count <- paste0(sum(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE), "/", length(before.control.index))
  	    after.treated.count <- paste0(sum(cov.data[after.treated.index,x] == x.val, na.rm = TRUE), "/", length(after.treated.index))
  	    after.control.count <- paste0(sum(cov.data[after.control.index,x] == x.val, na.rm = TRUE), "/", length(after.control.index))
  	   
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Count"] <- before.treated.count
  	    match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Count"] <- before.control.count
  	    match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Count"] <- after.treated.count
          match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Count"] <- after.control.count 	   
  	   
  	  }
  	  # Now deal with missingness
  	  before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	  after.control.index <- which(cov.data["treated"] == 0 & !is.na(cov.data[,"stratum"]))
  	 
  	  before.treated.mean <- mean(cov.data[before.treated.index, paste0(x, ".missing")], na.rm = TRUE)
  	  before.control.mean <- mean(cov.data[before.control.index, paste0(x, ".missing")], na.rm = TRUE)
  	  after.treated.mean <- wtd.mean(cov.data[after.treated.index, paste0(x, ".missing")], weight = cov.data[after.treated.index, "weight"], na.rm = TRUE)
  	  after.control.mean <- wtd.mean(cov.data[after.control.index, paste0(x, ".missing")], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
  	 
  	  combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,paste0(x, ".missing")], na.rm = TRUE) + 0.5*var(cov.data[after.treated.index,paste0(x, ".missing")], na.rm = TRUE))
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	 
  	  before.treated.count <- paste0(sum(cov.data[before.treated.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(before.treated.index))
  	  before.control.count <- paste0(sum(cov.data[before.control.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(before.control.index))
  	  after.treated.count <- paste0(sum(cov.data[after.treated.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(after.treated.index))
  	  after.control.count <- paste0(sum(cov.data[after.control.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(after.control.index))
 
  	  match.table[match.table[,"X.val"] == "missing","Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "missing","Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "missing","Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "missing","Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "missing","After.Treated.Count"] <- after.treated.count
  	  match.table[match.table[,"X.val"] == "missing","After.Control.Count"] <- after.control.count
    } else{   # if not, match.table table requires 1 + length(cat.values) rows in the match.table
      unik.values <- sort(unique(cov.data[,x]))
      n.values <- length(unik.values)
  	match.table <- data.frame("X" = rep(x, times = n.values), 
  	                          "X.val" = unik.values,
  	                          "Before.Treated.Mean" = rep(NA, times = n.values),
  	                          "After.Treated.Mean" = rep(NA, times = n.values),
  	                          "Before.Control.Mean" = rep(NA, times = n.values),
  	                          "After.Control.Mean" = rep(NA, times = n.values),
  	                          "Std.Diff.Before" = rep(NA, times = n.values),
  	                          "Std.Diff.After" = rep(NA, times = n.values),
  	                          "Before.Treated.Count" = rep(NA, times = n.values),
  	                          "Before.Control.Count" = rep(NA, times = n.values),
  	                          "After.Treated.Count" = rep(NA, times = n.values),
  	                          "After.Control.Count" = rep(NA, times = n.values))
  	
      # Deal with the categorical variables
      before.treated.index <- which(cov.data[,"treated"] == 1)
  	before.control.index <- which(cov.data[,"treated"] == 0)
  	after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	after.control.index <- which(cov.data[,"treated"] == 0 & !is.na(cov.data[,"stratum"]))
  	 
  	 
  	for(x.val in unik.values){
  	  before.treated.mean <- mean(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE)
  	  before.control.mean <- mean(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE)
  	   
  	  after.treated.mean <- wtd.mean(1*(cov.data[after.treated.index,x] == x.val), weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
  	  after.control.mean <- wtd.mean(1*(cov.data[after.control.index,x] == x.val), weight = cov.data[after.control.index, "weight"], na.rm = TRUE)
  	   
  	  combined.sd <- sqrt(0.5*var(1*(cov.data[before.treated.index,x] == x.val), na.rm = TRUE) + 0.5*var(1*(cov.data[before.control.index,x] == x.val), na.rm = TRUE))
  	   
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	   
  	  before.treated.count <- paste0(sum(cov.data[before.treated.index,x] == x.val, na.rm = TRUE), "/", length(before.treated.index))
  	  before.control.count <- paste0(sum(cov.data[before.control.index,x] == x.val, na.rm = TRUE), "/", length(before.control.index))
  	  after.treated.count <- paste0(sum(cov.data[after.treated.index,x] == x.val, na.rm = TRUE), "/", length(after.treated.index))
  	  after.control.count <- paste0(sum(cov.data[after.control.index,x] == x.val, na.rm = TRUE), "/", length(after.control.index))
  	   
  	  match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == as.character(x.val),"Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == as.character(x.val),"Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == as.character(x.val),"Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == as.character(x.val),"After.Treated.Count"] <- after.treated.count
        match.table[match.table[,"X.val"] == as.character(x.val),"After.Control.Count"] <- after.control.count 	    
      }
	
    } # closes else statement seeing whether x.missing is in colnames(cov.data)
  
    return(match.table)
  # Deal with the catergorical values
  } else { # closes if statement checking whether x is in colnames(cov.data)
  	return(NULL)
  }
}


assess_balance_cont <- function(x, cov.data){

  if(x %in% colnames(cov.data)){
  	
  	if(paste0(x, ".missing") %in% colnames(cov.data)){
  	  match.table <- data.frame("X" = rep(x, times = 2), 
  	                            "X.val" = c("continuous", "missing"),
  	                            "Before.Treated.Mean" = rep(NA, times = 2),
  	                            "After.Treated.Mean" = rep(NA, times = 2),
  	                            "Before.Control.Mean" = rep(NA, times = 2),
  	                            "After.Control.Mean" = rep(NA, times = 2),
  	                            "Std.Diff.Before" = rep(NA, times = 2),
  	                            "Std.Diff.After" = rep(NA, times = 2),
  	                            "Before.Treated.Count" = rep(NA, times = 2),
  	                            "Before.Control.Count" = rep(NA, times = 2),
  	                            "After.Treated.Count" = rep(NA, times = 2),
  	                            "After.Control.Count" = rep(NA, times = 2))
  	  # Now deal with missingness
  	  before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  	  after.control.index <- which(cov.data["treated"] == 0 & !is.na(cov.data[,"stratum"]))
  	 
  	  before.treated.mean <- mean(cov.data[before.treated.index, paste0(x, ".missing")], na.rm = TRUE)
  	  before.control.mean <- mean(cov.data[before.control.index, paste0(x, ".missing")], na.rm = TRUE)
  	  after.treated.mean <- wtd.mean(cov.data[after.treated.index, paste0(x, ".missing")], weight = cov.data[after.treated.index, "weight"], na.rm = TRUE)
  	  after.control.mean <- wtd.mean(cov.data[after.control.index, paste0(x, ".missing")], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
  	 
  	  combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,paste0(x, ".missing")], na.rm = TRUE) + 0.5*var(cov.data[after.treated.index,paste0(x, ".missing")], na.rm = TRUE))
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	 
  	  before.treated.count <- paste0(sum(cov.data[before.treated.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(before.treated.index))
  	  before.control.count <- paste0(sum(cov.data[before.control.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(before.control.index))
  	  after.treated.count <- paste0(sum(cov.data[after.treated.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(after.treated.index))
  	  after.control.count <- paste0(sum(cov.data[after.control.index,paste0(x, ".missing")], na.rm = TRUE), "/", length(after.control.index))
 
  	  match.table[match.table[,"X.val"] == "missing","Before.Treated.Mean"] <- round(100*before.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","Before.Control.Mean"] <- round(100*before.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","After.Treated.Mean"] <- round(100*after.treated.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","After.Control.Mean"] <- round(100*after.control.mean, digits = 2)
  	  match.table[match.table[,"X.val"] == "missing","Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "missing","Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "missing","Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "missing","Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "missing","After.Treated.Count"] <- after.treated.count
  	  match.table[match.table[,"X.val"] == "missing","After.Control.Count"] <- after.control.count

  	  # Now deal with the continuous
  	  # we want to take averages only of the treated and controls that have continuous values and are not missing
  	  before.treated.index <- which(cov.data[,"treated"] == 1 & cov.data[,paste0(x, ".missing")] == 0)
  	  before.control.index <- which(cov.data[,"treated"] == 0 & cov.data[,paste0(x, ".missing")] == 0)
  	 
  	  after.treated.index <- which(cov.data[,"treated"] == 1 &  cov.data[,paste0(x, ".missing")] == 0 & !is.na(cov.data[,"stratum"]))
        after.control.index <- which(cov.data[,"treated"] == 0 &  cov.data[,paste0(x, ".missing")] == 0 & !is.na(cov.data[,"stratum"]))
      
        before.treated.mean <- mean(cov.data[before.treated.index,x], na.rm = TRUE)
        before.control.mean <- mean(cov.data[before.control.index,x], na.rm = TRUE)
      
        after.treated.mean <- wtd.mean(cov.data[after.treated.index,x], weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
        after.control.mean <- wtd.mean(cov.data[after.control.index,x], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
      
        combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,x], na.rm = TRUE) + 0.5*var(cov.data[before.control.index,x], na.rm = TRUE))
  	  
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	  
  	  before.treated.count <- NA
  	  before.control.count <- NA
  	  after.treated.count <- NA
  	  after.control.count <- NA
  	  
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Mean"] <- round(before.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Mean"] <- round(before.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Mean"] <- round(after.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Mean"] <- round(after.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Count"] <- after.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Count"] <- after.control.count  

  		
  	} else{ # closes if statement checking whether x.missing is in colnames(cov.data)
  	  match.table <- data.frame("X" = rep(x, times = 1), 
  	                            "X.val" = c("continuous"),
  	                            "Before.Treated.Mean" = rep(NA, times = 1),
  	                            "After.Treated.Mean" = rep(NA, times = 1),
  	                            "Before.Control.Mean" = rep(NA, times = 1),
  	                            "After.Control.Mean" = rep(NA, times = 1),
  	                            "Std.Diff.Before" = rep(NA, times = 1),
  	                            "Std.Diff.After" = rep(NA, times = 1),
  	                            "Before.Treated.Count" = rep(NA, times = 1),
  	                            "Before.Control.Count" = rep(NA, times = 1),
  	                            "After.Treated.Count" = rep(NA, times = 1),
  	                            "After.Control.Count" = rep(NA, times = 1))
  	                            
  	  before.treated.index <- which(cov.data[,"treated"] == 1)
  	  before.control.index <- which(cov.data[,"treated"] == 0)
  	 
  	  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
        after.control.index <- which(cov.data[,"treated"] == 0 & !is.na(cov.data[,"stratum"]))
        before.treated.mean <- mean(cov.data[before.treated.index,x], na.rm = TRUE)
        before.control.mean <- mean(cov.data[before.control.index,x], na.rm = TRUE)
      
        after.treated.mean <- wtd.mean(cov.data[after.treated.index,x], weight = cov.data[after.treated.index,"weight"], na.rm = TRUE)
        after.control.mean <- wtd.mean(cov.data[after.control.index,x], weight = cov.data[after.control.index,"weight"], na.rm = TRUE)
      
        combined.sd <- sqrt(0.5*var(cov.data[before.treated.index,x], na.rm = TRUE) + 0.5*var(cov.data[before.control.index,x], na.rm = TRUE))
  	  
  	  std.diff.before <- (before.treated.mean - before.control.mean)/combined.sd
  	  std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
  	  
  	  before.treated.count <- NA
  	  before.control.count <- NA
  	  after.treated.count <- NA
  	  after.control.count <- NA
  	  
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Mean"] <- round(before.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Mean"] <- round(before.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Mean"] <- round(after.treated.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Mean"] <- round(after.control.mean, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.Before"] <- round(std.diff.before, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Std.Diff.After"] <- round(std.diff.after, digits = 3)
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Treated.Count"] <- before.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "Before.Control.Count"] <- before.control.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Treated.Count"] <- after.treated.count
  	  match.table[match.table[,"X.val"] == "continuous", "After.Control.Count"] <- after.control.count	
  	}
  	return(match.table)
  } else{ # closes if statement checking whether x is in colnames(cov.data)
  	return(NULL)
  }


}


# mix.list: names of the variables that have categorical and continuous values
# mix.values is a list whose elements are given in mix.list, containing the categorical values
assess_balance <- function(cov.data, cont.vars, cat.vars = c(), mix.vars = c(), mix.values = list()){

  full.table <- data.frame("X" = NA, "X.val" = NA,
                           "Before.Treated.Mean" = NA,
                           "After.Treated.Mean" = NA,
                           "Before.Control.Mean" = NA,
                           "After.Control.Mean" = NA,
                           "Std.Diff.Before" = NA,
                           "Std.Diff.After" = NA,
                           "Before.Treated.Count" = NA,
                           "Before.Control.Count" = NA,
                           "After.Treated.Count" = NA,
                           "After.Control.Count" = NA)
  for(x in cont.vars){
    #print(x)
    tmp <- assess_balance_cont(x, cov.data)
    if(!is.null(tmp)){
      full.table <- rbind(full.table, tmp)
    }
  }
  for(x in cat.vars){
    #print(x)
    tmp <- assess_balance_cat(x, cov.data)
    if(!is.null(tmp)){
      full.table <- rbind(full.table, tmp)
    }
  }
  for(x in mix.vars){
    #print(x)
    tmp <- assess_balance_mix(x, cov.data, mix.values[[x]])
    if(!is.null(tmp)){
      full.table <- rbind(full.table, tmp)
    }
  }
  full.table <- full.table[-1,]

  
  
  # Now return three tables:
  missing.index <- which(full.table[,"X.val"] == "missing")
  if(length(missing.index > 1)){
  	table.red <- full.table[-missing.index,]
  	table.missing <- full.table[missing.index,]
  } else{
  	table.red <- NULL
  	table.missing <- NULL
  }

  results <- list("full" = full.table, "reduced" = table.red, "missing" = table.missing)
  return(results)
}



