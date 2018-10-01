# Homework 1 Supporting R Code

library(reshape2)
library(ggplot2)
library(stats)
library(car)
library(dplyr)

hw1 <- read.csv("~/WORKING_DIRECTORIES/biostat.653/hw1.csv")
hw1 <- hw1[complete.cases(hw1),]

plot01 <- ggplot(data = hw1, aes(x = wt, y = arm)) + 
  geom_line(aes(group = factor(id))) +
  geom_smooth() +
  facet_grid(. ~ sex) +
  ylab("Arm Circumference") +
  xlab("Weight") +
  theme_light()
plot01

model_2 <- lm(arm ~ wt + sex + age, data = hw1)
avPlots(model_2, terms = ~.)

hw1 <- mutate(hw1, age_sq = age^2)

model_3 <- lm(arm ~ wt + sex + age + age_sq, data = hw1)
data_3 <- data.frame(cbind(model_3$residuals,hw1$wt,hw1$sex))
plot03 <- ggplot(data = data_3, aes(x = X2, y = X1)) +
  geom_point() +
  geom_smooth() +
  ylab("Residuals") +
  xlab("Weight")
plot03

plot04 <- plot03 + facet_grid(. ~ X3)
plot04