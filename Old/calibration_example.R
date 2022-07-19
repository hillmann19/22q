# example of Tyler's calibration script (for CNB task vaildity)

# Akira Di Sandro, 4.7.22

# load packages ----
library(PerFit)

# load data ----

ADT30_1 <- read.csv("/Users/hillmann/Projects/22q/Data/forNoah_ADT30_Scored.csv")     # ADT = non-memory task
CPF_1 <- read.csv("/Users/hillmann/Projects/22q/Data/forNoah_CPF_Scored.csv")         # CPF = memory task

# write.csv(ADT30_1, "../forNoah_ADT30_Scored.csv",row.names = F)
# write.csv(CPF_1, "../forNoah_CPF_Scored.csv",row.names = F)

# validation ----

dat <- ADT30_1[,c(grep("_CORR",colnames(ADT30_1)),grep("_TTR",colnames(ADT30_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
  mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
  res[,j] <- scale(residuals(mod,na.action=na.exclude))
}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)

ADT30_1 <- data.frame(ADT30_1,sc)

x <- ADT30_1[,grepl("_CORR",colnames(ADT30_1)) | grepl("PFscores",colnames(ADT30_1))]
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
x <- x[,colnames(x) != "PFscores"]

dat <-CPF_1[,c(grep("_CORR",colnames(CPF_1)),grep("_TTR",colnames(CPF_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
  mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
  res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPF_1<- data.frame(CPF_1,sc)

x <-CPF_1[,grepl("_CORR",colnames(CPF_1)) | grepl("PFscores",colnames(CPF_1))]
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
x <- x[,colnames(x) != "PFscores"]



