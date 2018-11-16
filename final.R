# Final Project - Older Korean Depression Data
library(readxl)
korea <- read_excel("~/WORKING_DIRECTORIES/biostat.653/Older_Korean_Depression_Data.xls")
# Dependent Variable: Depression Score
# Covariates:         Age
#                     Gender
#                     High School Graduate
#                     College Graduate
#                     Employed
#                     Income
#                     BMI

# Look into missing data
library(mice)
library(VIM)

# Examine missing data patterns in the full dataset
md.pattern(korea)

# Examine missing data patterns in specific years
korea06 <- korea[1:9]
korea08 <- korea[c(1,10:17)]
korea10 <- korea[c(1,18:25)]
md.pattern(korea06)
md.pattern(korea08)
md.pattern(korea10)

# Removes the observations missing full years of data
korea_miss <- is.na(korea)
korea2 <- korea[rowSums(korea_miss) < 7,]

# Examine missing data patterns in the reduced dataset
md.pattern(korea2)

# Performing the single imputation with K-nearest neighbors
seed(11162018)
korea_Imp <- kNN(korea2,
                 # Variables with missing values to impute
                 variable = c("hsgrad2006collgrad2006","hsgrad2008","hsgrad2010",
                 "collgrad2006","collgrad2008","collgrad2010",
                 "yrincome2006","yrincome2008","yrincome2010",
                 "cesd10_total2006","cesd10_total2008","cesd10_total2010",
                 "obese42006","obese42008","obese42010"),
                 # Variables to use in the distance measure calculations
                 dist_var = c("age2006","age2008","age2010",
                 "female2006","female2008","female2010",
                 "employed2006","employed2008","employed2010",
                 "hsgrad2006collgrad2006","hsgrad2008","hsgrad2010",
                 "collgrad2006","collgrad2008","collgrad2010",
                 "yrincome2006","yrincome2008","yrincome2010",
                 "cesd10_total2006","cesd10_total2008","cesd10_total2010",
                 "obese42006","obese42008","obese42010"), 
                 # Number of neighbors
                 k = 5,
                 # numFun = median
                 # catFun = maxCat # chooses the level with the most
                 #   occurences and random if the maximum is not unique
                 useImputedDist = FALSE)

# The imputed dataset to be used for analysis
library(foreign)
write.foreign(korea_Imp,
              "~/WORKING_DIRECTORIES/biostat.653/korea_Imp.txt",
              "~/WORKING_DIRECTORIES/biostat.653/korea_Imp.sas",
              package = "SAS")