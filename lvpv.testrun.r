### TEST RUN
source("lvpv.r")

## create TestData dataset
library(evolqg)
library(MASS)

set.seed(254)
m = rnorm(8, mean=10, sd=4)
s = m*0.1 + rnorm(8, mean = 0, sd = 1)
V = s^2
VCV = RandomMatrix(8, variance=V, LKJ=FALSE)
TestData = mvrnorm(n=240, mu=m, Sigma=VCV)
colnames(TestData) = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8")

## analyze TestData using Newman's leading eigenvector method
pvTD = PV.complete(TestData, iseed=767)
PV.dendro(pvTD, print.pv = F); PV.text(pvTD) 
plot(PV.graph(pvTD))
PV.graph.highlight(pvTD, col="lightblue")
PV.summary(pvTD)

## compare to the results of cluster analysis and Sammon's mapping

library(pvclust)
hclTD = pvclust(TestData, method.hclust = "average")
plot(hclTD)
pvrect(hclTD)

sTD = sammon(as.dist(1-cor(TestData)), k = 2)
plot(sTD$points[,1], sTD$points[,2], bg="springgreen", pch=21, cex=2, xlab="Axis 1", ylab="Axis 2"); text(sTD$points[,1], sTD$points[,2], label=rownames(sTD$points), pos=4)

# the results of cluster analysis and Sammon's method differ from the results of Newman's leading eigenvector method because v5 negatively correlates with other variables

# to obtain the results of cluster analysis and Sammon's method similar to the results of Newman's leading eigenvector method, use
TestData5 = TestData
TestData5[,5] = TestData[,5]*(-1)

hclTD5 = pvclust(TestData5, method.hclust = "average")
plot(hclTD5)
pvrect(hclTD5)

sTD5 = sammon(as.dist(1-cor(TestData5)), k = 2)
plot(sTD5$points[,1], sTD5$points[,2], bg="springgreen", pch=21, cex=2, xlab="Axis 1", ylab="Axis 2"); text(sTD5$points[,1], sTD5$points[,2], label=rownames(sTD5$points), pos=4)

# to obtain the results of Newman's leading eigenvector method similar to the results of cluster analysis and Sammon's method, use
pvTD.n = PV.complete(TestData, iseed=767, square=FALSE)
PV.dendro(pvTD.n, print.pv = F); PV.text(pvTD.n)
PV.graph.highlight(pvTD.n, col="lightblue")
PV.summary(pvTD.n)

## compare to the modularity detection algorithm realized in evolqg

# use evolqg
lmod = LModularity(cor(TestData))
lmod.partition = integer(dim(TestData)[2])
for (j in 1:dim(TestData)[2]) lmod.partition[j] = which(lmod$Modularity_hypothesis[j,]==1)
names(lmod.partition) = colnames(TestData)
lmod$LModularity
lmod$Modularity_hypothesis
lmod.partition

# use lvpv, specify the same calculation algorithm as in evolqg 
pvTD.diag = PV.complete(TestData, iseed=767, q.cut=0, diag=TRUE)
PV.dendro(pvTD.diag)
PV.graph.highlight(pvTD.diag, col="lightblue")
PV.summary(pvTD.diag)

# compare
c(pvTD.diag$modularity, lmod$LModularity)
lmod.partition
pvTD.diag$membership

## use Gaussian transformation before analyzing the data
pvTD.gauss = PV.complete(TestData, iseed=767, method="kendall", gauss=TRUE)
PV.dendro(pvTD.gauss, print.pv = F); PV.text(pvTD.gauss) 
PV.graph.highlight(pvTD.gauss, col="lightblue")
PV.summary(pvTD.gauss)
