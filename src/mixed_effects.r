##############################################
##    Run Mixed Effects as part of grmpy    ##
##############################################


### Put Get coefficients functions here ###
#############
#####



### Run Mixed Effects ###
r_mixedeffects <- function(datacsv, ...) {
	
	## Libraries:
	library(lme4)

	## Read in file:
	db <- read.table(datacsv)
	
	## Set the V1 and V6 columns, columns of station names and event numbers, as a factor
	db$evnum <- as.factor(db$evnum)
	db$sta <- as.factor(db$sta)
	
	## Subset the data so no 0 or negative PGA:
	db <- subset(db, pga>0)
	
	##			  ##
	## Run Models ##
	##			  ##
	
	model <- lmer(pga ~ m + m2 + lnR + rrup + (1|evnum) + (1|sta), data=db)
	
	## Print to screen:
	print(model)
	
	## Also save coefficients to a csv file, ##NAME##
	fixed_coefs <- get_coefs(
**********	filename <- file.path(".","my_results_random.csv")
	readr::write_csv(fixed_coefs, path=filename)
	message("Fixed coefficients saved to: ", filename)
	
	site_coefs <- get_site(
	
	event_coefs <- get_event(
	
	
### Call Function ###

print("Running mixed effects model...")

