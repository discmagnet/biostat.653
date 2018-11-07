# Homework 3 Supporting R Code
library(haven)
rat <- read_sas("~/WORKING_DIRECTORIES/biostat.653/rat.sas7bdat")

# 6.1.1 - On a single graph, construct a time plot that displays the
#         mean weight versus time (in weeks) for the three groups.
#         Describe the general characteristics of the time trends
#         for the three groups.


# 6.1.2 - Read the data into "long" format, with five "records"
#         per subject.
library(reshape2)
library(dplyr)
rat_long <- melt(rat[3:7])
rat_long <- mutate(rat_long,
                   group = rep(c(rep(1,10),rep(2,10),rep(3,7)),5),
                   id = rep(c(1:27),5))
rat_long <- arrange(rat_long, id)

# 6.1.3 - Assume that the rate of increase in each group is approx-
#         imately constant throughout the duration of the study.
#         Assuming an unstructured covariance matrix, construct a
#         test of whether the rate of increase differs in the groups.


# 6.1.4 - On a single graph, construct a time plot that displays the
#         estimated mean weight versus time (in weeks) for the three
#         treatment groups from the results generated from 6.1.3.


# 6.1.5 - Based on the results from 6.1.3, what is the estimated
#         rate of increase in mean weight in the control group?

#         What is the estimated rate of increase in mean weight in
#         the thiouracil group (group 2)?

#         What is the estimated rate of increase in mean weight in
#         the thyroxin group (group 3)?

# 6.1.6 - The study investigators conjectured that there would be an
#         increase in weight, but that the rate of increase would
#         level off toward the end of the study. They also conjectured
#         that this pattern of change may differ in the three treat-
#         ment groups. Assuming an unstructured covariance matrix,
#         construct a test of this hypothesis.
