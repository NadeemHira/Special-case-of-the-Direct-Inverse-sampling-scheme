# reading the CYP-GUIDES Trial data from Harford Hospital
data1 <- read.csv("C://Users/Hira/Desktop/PHD Thesis/clinical data.csv", header = TRUE)

# Looking at the summary of data
summary(data1)

# Subsetting the data based on the variables of interest
clinical.dat <- data.frame(data1$Assignment, data1$RAR)

# Summary for the clinical data
summary(clinical.dat)

# Subsetting the data based on the course of treatments in Assignment 
# variable for "G" and "S"

# The RCT recruited 1500 patients, genotyped CYP2D6 in 1459, and 
# randomized 477 to standard therapy (Group S), for whom 
# treatment-as-usual guidance was delivered without consideration of 
# patient CYP2D6 genotype, and 982 to genetically-guided therapy (Group G)
#where CYP2D6-based treatment recommendations were provided 
#via EMR to physicians.

assG <- subset(clinical.dat, data1.Assignment=="G")
assS <- subset(clinical.dat, data1.Assignment=="S")

# Confidence Interval for the special case for treatments G and S
# Reduced data of "G" treatment by fixing the number of successes 
# in the "S" treatment, based on the RAR value of 1.
# Therefore, by looking at the number of successes in S, we reduce the 
# data set by 548 observations
red1 <- head(assG, -548)

# Counting the number of successes for the reduced data set of G treatment,
# where the success is RAR = 1
count1 <- length(which(red1$data1.RAR == 1))
count1

# T is the number of successes in the first sample
T <- count1

# Number of observations in the reduced data set
n <- nrow(red1)

# Number of rows with RAR = 1 in the second sample with "S" treatment
count2 <- length(which(assS$data1.RAR == 1))
count2

# Observations in the "S" treatment groups
nuT <- nrow(assS)
nuT

# Computing the rho hat values
rhohat <- (nuT - T)/(n+1-T)

# Tau-squared values based on the values in the first and second samples
tausquare <- ((nuT-T)/(n+1-T))*((((n-T)/(T+1))*((n+1)/(T+1)))+
                                  (((nuT/T)-1)*((n+1)/(T+1))^2))*((T/(n+1-T))^2)
# Standard error for tau squared values
s2 <- tausquare/n

# 95% confidence intervals 
lower <- rhohat - qnorm(0.975)*sqrt(s2)
upper <- rhohat + qnorm(0.975)*sqrt(s2)
CI <- c(lower, upper)
CI

# Gender-wise study of treatment for Readmission rate
clinical.dat2 <- data.frame(data1$GENDER,data1$Assignment, data1$RAR)

# Summary of the clinical data with three columns
summary(clinical.dat2)

# Dividing data based on gender = female and gender = male
fem.data <- subset(clinical.dat2, data1.GENDER == "F")
male.data <- subset(clinical.dat2, data1.GENDER == "M")

# Dividing the gender-wise data based on treatments "G" and "S"
fem.data.assG <- subset(fem.data, data1.Assignment == "G")
fem.data.assS <- subset(fem.data, data1.Assignment == "S")

male.data.assG <- subset(male.data, data1.Assignment == "G")
male.data.assS <- subset(male.data, data1.Assignment == "S")

# Female reduced data based on treatments
# Reducing the first data set for "G" treatment based on the RAR = 1 
# for the second sample
# Reducing the female data by 351 observations for "G" treatment
red2 <- head(fem.data.assG, -351)

# Counting the observations for the reduced data set for "G" treatment in
# females for which the RAR = 1
count3 <- length(which(red2$data1.RAR == 1))
count3

# number of successes
T <- count3

# number of observations in the reduced data
n <- nrow(red2)

# Number of successes in the second sample with "S" treatment and RAR = 1
count4 <- length(which(fem.data.assS$data1.RAR == 1))
count4

# Number of observations in the second sample
nuT <- nrow(fem.data.assS)
nuT

# Rho hat estimate for female receiving "G" and "S" treatments
rhohat <- (nuT - T)/(n+1-T)

# Tau^2 estimate for female
tausquare <- ((nuT-T)/(n+1-T))*((((n-T)/(T+1))*((n+1)/(T+1)))+
                                  (((nuT/T)-1)*((n+1)/(T+1))^2))*((T/(n+1-T))^2)

# Standard error for tau^2
s2 <- tausquare/n

# 95% confidence interval for females receiving "G" and "S" treatments
lower <- rhohat - qnorm(0.975)*sqrt(s2)
upper <- rhohat + qnorm(0.975)*sqrt(s2)
CI <- c(lower, upper)
CI


# Male reduced data based on treatments
# Reducing the first data set for "G" treatment based on the RAR = 1 
# for the second sample
# Reducing the male data by 141 observations for "G" treatment
red3 <- head(male.data.assG, -141)

# Counting the observations for the reduced data set for "G" treatment in
# males for which the RAR = 1
count5 <- length(which(red3$data1.RAR == 1))
count5

# number of successes
T <- count5

# number of observations in the reduced data
n <- nrow(red3)

# Number of successes in the second sample with "S" treatment and RAR = 1
count6 <- length(which(male.data.assS$data1.RAR == 1))
count6

# Number of observations for the second sample
nuT <- nrow(male.data.assS)
nuT

# Rho hat estimate for males receiving "G" and "S" treatments
rhohat <- (nuT - T)/(n+1-T)

# Tau^2 estimate for males
tausquare <- ((nuT-T)/(n+1-T))*((((n-T)/(T+1))*((n+1)/(T+1)))+
                                  (((nuT/T)-1)*((n+1)/(T+1))^2))*((T/(n+1-T))^2)

# Standard error estimates for males
s2 <- tausquare/n

# 95% confidence intervals for males receiving "G" and "S" treatments
lower <- rhohat - qnorm(0.975)*sqrt(s2)
upper <- rhohat + qnorm(0.975)*sqrt(s2)
CI <- c(lower, upper)
CI

# Age wise study of treatments for the Readmission rates based on 
# people under 40 and over 40 years of age (based on median age)
clinical.dat3 <- data.frame(data1$AGE, data1$Assignment, data1$RAR)
summary(clinical.dat3)

# median Age
median(data1$AGE)

# subsetting data based on  the age of people under and over 40 years
age1 <- subset(clinical.dat3, data1.AGE <= 40)
age2 <- subset(clinical.dat3, data1.AGE > 40)

# Dividing the age-wise data based on treatments "G" and "S"
age1.assG <- subset(age1, data1.Assignment == "G")
age1.assS <- subset(age1, data1.Assignment == "S")

age2.assG <- subset(age2, data1.Assignment == "G")
age2.assS <- subset(age2, data1.Assignment == "S")

# Age reduced data based on treatments
# Reducing the first data set for "G" treatment based on the RAR = 1 
# for the second sample
# Reducing the age1 data by 236 observations for "G" treatment
red4 <- head(age1.assG, -236)

# Counting the observations for the reduced data set for "G" treatment in
# age1 for which the RAR = 1
count7 <- length(which(red4$data1.RAR == 1))
count7

# Number of successes
T <- count7

# Number of observations of the reduced data set
n <- nrow(red4)

# Number of successes in the second sample with "S" treatment and RAR = 1
count8 <- length(which(age1.assS$data1.RAR == 1))
count8

# Number of observations for the second sample
nuT <- nrow(age1.assS)
nuT

# Rho hat estimate for under 40 people, receiving "G" and "S" treatments
rhohat <- (nuT - T)/(n+1-T)

# Tau^2 estimate for people under the age of 40
tausquare <- ((nuT-T)/(n+1-T))*((((n-T)/(T+1))*((n+1)/(T+1)))+
                                  (((nuT/T)-1)*((n+1)/(T+1))^2))*((T/(n+1-T))^2)

# Standard error estimates for people under 40 years of age
s2 <- tausquare/n

# 95% confidence intervals for people under 40 who receive "G" and "S"
# treatments
lower <- rhohat - qnorm(0.975)*sqrt(s2)
upper <- rhohat + qnorm(0.975)*sqrt(s2)
CI <- c(lower, upper)
CI

# Reducing the age2 data by 275 observations for "G" treatment
red5 <- head(age2.assG, -275)

# Counting the observations for the reduced data set for "G" treatment in
# age2 for which the RAR = 1
count9 <- length(which(red5$data1.RAR == 1))
count9

# Number of successes in reduced data of "G" treatment
T <- count9

# Number of observations of the reduced data set
n <- nrow(red5)

# Number of successes in the second sample with "S" treatment and RAR = 1
count10 <- length(which(age2.assS$data1.RAR == 1))
count10

# Number of observations for the second sample
nuT <- nrow(age2.assS)
nuT

# Rho hat estimate for over 40 people, receiving "G" and "S" treatments
rhohat <- (nuT - T)/(n+1-T)

# Tau^2 estimate for people over 40, receiving "G" and "S" treatments
tausquare <- ((nuT-T)/(n+1-T))*((((n-T)/(T+1))*((n+1)/(T+1)))+
                                  (((nuT/T)-1)*((n+1)/(T+1))^2))*((T/(n+1-T))^2)

# Standard error for rho hat for people over 40 years of age 
s2 <- tausquare/n

# 95% confidence intervals for people over 40 who receive "G" and "S"
# treatments
lower <- rhohat - qnorm(0.975)*sqrt(s2)
upper <- rhohat + qnorm(0.975)*sqrt(s2)
CI <- c(lower, upper)
CI
