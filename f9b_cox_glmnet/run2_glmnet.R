# module load gcc/7.1.0  openmpi/3.1.4
# module load R/4.1.1
# use local R in mac
# install.packages("glmnet", repos = "https://cran.us.r-project.org")

library(glmnet)
library(Matrix)


x <- read.csv("data/top31TF_RPKM.csv", sep=",",row.names=1) 
x <- as.matrix(x)
# x <- read.csv("data/top31TF_CPM.csv", sep=",",row.names=1) 
# x <- as.matrix(x)

y <- read.csv("data/clinical_data.csv", sep=",",row.names=1) 
y <- as.matrix(y)


# fit glmnet
fit <- glmnet(x, y, family = "cox")
pdf(file=paste("f2_figs/fit.pdf",sep=""))
plot(fit,xvar='lambda')
dev.off()
coef <- coef(fit,s=0.08)
write.table(as.matrix(coef),'f2_figs/coef_s0p08.csv',col.names=FALSE,sep=",")

# cv for the best lambda/s
for (cvi in 1:10){
pdf(file=paste("f2_figs/cvfit_",cvi,".pdf",sep=""))
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C",nfolds=5)
plot(cvfit)
dev.off()
}


# tf9 <- c('GATA1','RAD21','CUX1','KLF16','ADNP','ZEB2','TEAD4','STAG1','SIRT6')
# xx <- x[,tf9]
# 

