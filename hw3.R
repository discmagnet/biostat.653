# Homework 3 Supporting R Code
library(haven)
rat <- read_sas("~/WORKING_DIRECTORIES/biostat.653/rat.sas7bdat")

# 6.1.1 - On a single graph, construct a time plot that displays the
#         mean weight versus time (in weeks) for the three groups.
#         Describe the general characteristics of the time trends
#         for the three groups.
library(dplyr)
group1 <- subset(rat, GROUP == 1)
group2 <- subset(rat, GROUP == 2)
group3 <- subset(rat, GROUP == 3)
summary1 <- summarise(group1, mean_y1 = mean(Y1),
                              mean_y2 = mean(Y2),
                              mean_y3 = mean(Y3),
                              mean_y4 = mean(Y4),
                              mean_y5 = mean(Y5))
summary2 <- summarise(group2, mean_y1 = mean(Y1),
                      mean_y2 = mean(Y2),
                      mean_y3 = mean(Y3),
                      mean_y4 = mean(Y4),
                      mean_y5 = mean(Y5))
summary3 <- summarise(group3, mean_y1 = mean(Y1),
                      mean_y2 = mean(Y2),
                      mean_y3 = mean(Y3),
                      mean_y4 = mean(Y4),
                      mean_y5 = mean(Y5))
library(ggplot2)
data611 <- data.frame(rbind(t(summary1),t(summary2),t(summary3)),
                      c(rep(1,5),rep(2,5),rep(3,5)),c(0:4,0:4,0:4))
colnames(data611) <- c("mean","group","index")
plot611 <- ggplot(data = data611,
                  aes(x = index, y = mean, color = as.factor(group))) +
  geom_line() +
  geom_point() +
  xlab("Follow-up Time (in weeks)") +
  ylab("Weight (in grams)")
plot611

# 6.1.2 - Read the data into "long" format, with five "records"
#         per subject.
library(reshape2)
rat_long <- melt(rat[3:7])
rat_long <- mutate(rat_long,
                   group = rep(c(rep(1,10),rep(2,10),rep(3,7)),5),
                   id = rep(c(1:27),5))
rat_long <- arrange(rat_long, id)
rat_long <- mutate(rat_long,
                   time = rep(c(0:4),27),
                   knot1 = (time-1)*(time > 0))

# 6.1.3 - Assume that the rate of increase in each group is approx-
#         imately constant throughout the duration of the study.
#         Assuming an unstructured covariance matrix, construct a
#         test of whether the rate of increase differs in the groups.
library(nlme)
res.gls <- gls(model = value ~ factor(group)*time,
               data = rat_long,
               # Covariance structure: unstructured
               correlation = corSymm(form = ~ as.numeric(factor(variable))|id),
               weights = varIdent(form = ~ 1 | factor(variable)),
               method = "ML")
summary(res.gls)
drop1(res.gls, test = "Chisq")

# 6.1.4 - On a single graph, construct a time plot that displays the
#         estimated mean weight versus time (in weeks) for the three
#         treatment groups from the results generated from 6.1.3.
data614 <- data.frame(c(res.gls$fitted[1:5],res.gls$fitted[51:55],
                      res.gls$fitted[101:105]),
                      c(rep(1,5),rep(2,5),rep(3,5)),c(0:4,0:4,0:4))
colnames(data614) <- c("mean","group","index")
plot614 <- ggplot(data = data614,
                  aes(x = index, y = mean, color = as.factor(group))) +
  geom_line() +
  geom_point() +
  xlab("Follow-up Time (in weeks)") +
  ylab("Weight (in grams)")
plot614

# 6.1.5 - Based on the results from 6.1.3, what is the estimated
#         rate of increase in mean weight in the control group?
#         What is the estimated rate of increase in mean weight in
#         the thiouracil group (group 2)?
#         What is the estimated rate of increase in mean weight in
#         the thyroxin group (group 3)?
res.gls$coefficients[4] # control group (group 1)
res.gls$coefficients[4] + res.gls$coefficients[5] # thiouracil group (group 2)
res.gls$coefficients[4] + res.gls$coefficients[6] # thyroxin group (group 3)


# 6.1.6 - The study investigators conjectured that there would be an
#         increase in weight, but that the rate of increase would
#         level off toward the end of the study. They also conjectured
#         that this pattern of change may differ in the three treat-
#         ment groups. Assuming an unstructured covariance matrix,
#         construct a test of this hypothesis.
res.gls <- gls(model = value ~ factor(group)*time + factor(group)*knot1,
               data = rat_long,
               # Covariance structure: unstructured
               correlation = corSymm(form = ~ as.numeric(factor(variable))|id),
               weights = varIdent(form = ~ 1 | factor(variable)),
               method = "ML")
summary(res.gls)
drop1(res.gls, test = "Chisq")