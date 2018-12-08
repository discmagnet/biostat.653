# Homework 4 Supporting R Code

# Question 1

# import data
library(readr)
exercise <- read_csv("~/WORKING_DIRECTORIES/biostat.653/exercise-data.csv")

# plot mean strength vs time for the two treatment groups
treatA <- subset(exercise, Treatment == 1)
treatB <- subset(exercise, Treatment == 2)
meanA <- vector()
meanB <- vector()
for(i in 1:7){
  meanA[i] <- colMeans(treatA[,i+2], na.rm = T)
  meanB[i] <- colMeans(treatB[,i+2], na.rm = T)
}
time <- c(0,2,4,6,8,10,12)
plot_data <- data.frame(c(meanA,meanB),
                        c(time,time),
                        c(rep(1,7),rep(2,7)))
colnames(plot_data) <- c("mean","time","group")
library(ggplot2)
plot01 <- ggplot(plot_data, aes(x=time,y=mean,col=as.factor(group))) +
  geom_point() +
  geom_line() +
  xlab("Time (in days)") +
  ylab("Measure of Strength") +
  ggtitle("Mean Strength vs Time") +
  scale_color_discrete(name="Program",labels=c("1","2"))
plot01

# convert to long form
library(reshape2)
library(dplyr)
exercise_long <- melt(exercise[3:9], value.name = "strength")
exercise_long <- mutate(exercise_long,
                        program = rep(c(rep(1,16),rep(2,21)),7),
                        id = rep(c(1:37),7))
exercise_long <- arrange(exercise_long, id)
exercise_long <- mutate(exercise_long, time = rep(2*0:6,37))

# fit a model with random intercepts and slopes
library(lme4)
library(broom)
model01a <- lmer(data = exercise_long,
                strength ~ factor(program) + 
                  time + factor(program):time + 
                  (1+time|id))
summary(model01a)

# fit a model with just a random intercept
model01b <- lmer(data = exercise_long,
                 strength ~ factor(program) + 
                   time + factor(program):time + 
                   (1|id))
summary(model01b)

# use LRT to test if a model with just random intercepts is adequate
lik_01a <- glance(model01a)$logLik
lik_01b <- glance(model01b)$logLik
LRT_stat <- -2*(lik_01b-lik_01a)
pchisq(LRT_stat,1,lower.tail = FALSE)

# test the effect of treatment on changes in strength
drop1(model01, test = "Chisq")

# obtain BLUP for each subject
ranef(model01)

# fit a simple linear model on subject 24
model02 <- lm(data = exercise_long,
              subset = (id == 24),
              strength ~ time)
summary(model02)

# Question 2

# import data
toenail <- read_csv("~/WORKING_DIRECTORIES/biostat.653/toenail-data.csv")

# fit a marginal model for the log odds of MOD or SEV onycholysis
library(geepack)
model03 <- geeglm(data = toenail,
                  id = ID,
                  corstr = "exchangeable",
                  family = binomial(link = "logit"),
                  Y ~ Month + factor(Treatment):Month)
summary(model03)

# fit a model with random intercepts
library(lme4)
model04 <- glmer(data = toenail,
                 family = binomial(link = "logit"),
                 nAGQ = 50, # sets the # of quadrature points
                 Y ~ Month + factor(Treatment):Month +
                   (1|ID)) # random intercept
summary(model04)
