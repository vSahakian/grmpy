##############################################
##    Run Mixed Effects as part of grmpy    ##
##############################################


############## Functions needed for mixed effects ##############

# Method to get fixed effects

get_fixed <- function(model, ...){ UseMethod("get_fixed") }

get_fixed.merMod <- function(model, ...){
	require(lme4)
	# get a data.frame of the fixed effects
	coef.tbl <- {
		modsum <- summary(model) #Can also use 'fixef' but this doesn't return std errors
		coefs <- modsum[['coefficients']]
		colnames(coefs) <- c('Estimate','Std.error','t.value')
		coefs <- cbind(data.frame(fixed.effect=rownames(coefs)), coefs)
		rownames(coefs) <- NULL
		return(coefs)
	}
	return(coef.tbl)
}

####################################################################

# Method to get the site random effects

get_site <- function(model, ...){ UseMethod("get_site") }

get_site.merMod <- function(model, ...){
	require(lme4)
	# get a data.frame of the random site effects
	modran <- lme4::ranef(model, condVar=TRUE)
	# internal function to add the std errors to the output
	.add_stderr <- function(x){
		id <- rownames(x)
		bias <- x$`(Intercept)`
		vars <- as.vector(attr(x, 'postVar'))
		std.err <- sqrt(vars)
		return(data.frame(ID=id, Bias=bias, Std.error=std.err))
	}
	
	coefs <- .add_stderr(modran[['sta']])
	return(coefs)
}
	
####################################################################

# Method to get the event random effects

get_event <- function(model, ...){ UseMethod("get_event") }

get_event.merMod <- function(model, ...){
	require(lme4)
	# get a data.frame of the random event effects
	modran <- lme4::ranef(model, condVar=TRUE)
	# internal function to add the std errors to the output
	.add_stderr <- function(x){
		id <- rownames(x)
		bias <- x$`(Intercept)`
		vars <- as.vector(attr(x, 'postVar'))
		std.err <- sqrt(vars)
		return(data.frame(ID=id, Bias=bias, Std.error=std.err))
	}
	
	coefs <- .add_stderr(modran[['evnum']])
	return(coefs)
}

	

###############################################################################
###############################################################################


### Run Mixed Effects ###
r_mixedeffects <- function(datacsv, home, database, ...) {
	
	# Input:
	#		datacsv: 		Full path to the csv file written by python
	#		home: 			String with path to home, no slash at end (i.e., /Users/vsahakian/anza or /home/vsahakian/katmai/anza)
	# 		database: 		String with database name (i.e., "test2013", or "abdb"
	#
	
	## Libraries:
	require(lme4)

	## Read in file:
	print(datacsv)
	db <- read.csv(datacsv)
	
	print("read in data")
	
	## Set the evnum and sta columns, columns of station names and event numbers, as a factor
	db$evnum <- as.factor(db$evnum)
	db$sta <- as.factor(db$sta)
	
	print("factored evnum and sta")
	
	## Subset the data so no 0 or negative PGA:
	#db <- subset(db, pga>0)
	
	##			  ##
	## Run Models ##
	##			  ##
	
	# Make
	
	
	model <- lme4::lmer(pga ~ m + m2 + lnR + rrup + (1|evnum) + (1|sta), data=db)
	
	# Make a new dataframe with just data to see if the predict funtion
	# 	always includes the random effects in the prediction:
	
	newdatframe <- data.frame(db$m, db$m2, db$lnR, db$rrup)
	
	prediction_vector <- predict(model)
	prediction <- data.frame(prediction_vector)
	
	## Print to screen:
	print(model)
	
	## Also save coefficients to a csv file, ##NAME##
	fixed_coefs <- get_fixed(model)
    filename <- file.path(home,"models/pckl",database,"r/results_fixed.csv")
	readr::write_csv(fixed_coefs, path=filename)
	message("Fixed coefficients saved to: ", filename)
	
	
	site_coefs <- get_site(model)
	filename <- file.path(home,"models/pckl",database,"r/results_site.csv")
	readr::write_csv(site_coefs, path=filename)
	message("Site effects saved to: ", filename)
	
	
	event_coefs <- get_event(model)
	filename <- file.path(home,"models/pckl",database,"r/results_event.csv")
	readr::write_csv(event_coefs, path=filename)
	message("Event effects saved to: ", filename)
	
	# Save predictions to a file:
	filename <- file.path(home,"models/pckl",database,"r/results_prediction.csv")
	readr::write_csv(prediction, path=filename)
	message("Wrote predicted values to: ", filename)
}	
	
	
### Call Function ###

# Echo...
print("Running mixed effects model...")

# Read in arguments:
args <- commandArgs(TRUE)

home_path <- eval(parse(text=args[1]))
database_name <- eval(parse(text=args[2]))

csvfile <- file.path(home_path,"models/pckl",database_name,"r","tmp_mixed.csv")

r_mixedeffects(csvfile,home_path,database_name)


