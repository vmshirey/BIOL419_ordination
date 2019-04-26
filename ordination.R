## Code for testing fourth-corner method with Cascadia butterflies ##

## Call packages ##
library(ade4)
library(vegan)

## Read in matrixes ##
R <- read.csv("EnvironmentalMatrix.csv", header = TRUE) # read environmental data
L <- read.csv("AbundanceMatrixNormalized.csv", header = TRUE) # read transect data
Q <- read.csv("SpeciesTraitMatrix.csv", header = TRUE) # read trait data

## Convert to appropriate format ##
R2 <- R[,-1]
rownames(R2) <- R[,1]
L2 <- L[,-1]
rownames(L2) <- L[,1]
Q2 <- Q[,-1]
rownames(Q2) <- Q[,1]

## Run RLQ analysis ##
coa1 <- dudi.coa(L2, scannf=FALSE, nf=2)
dudimil <- dudi.hillsmith(R2, scannf=FALSE, nf=2, row.w = coa1$lw)
duditrait <- dudi.hillsmith(Q2, scannf=FALSE, nf=2, row.w = coa1$cw)

rlq1 <- rlq(dudimil, coa1, duditrait, scannf=FALSE, nf=2) # conduct RLQ

plot(rlq1)

summary(rlq1)
randtest(rlq1)

par(mfrow=c(1,1))
s.arrow(rlq1$l1, boxes=FALSE, grid=FALSE) # environmental variables
s.arrow(rlq1$c1, boxes=FALSE, grid=FALSE) # traits
s.label(rlq1$lQ, boxes=FALSE, grid=FALSE, clabel=0.75) # species
s.label(rlq1$lR, boxes=FALSE, grid=FALSE) # sites

## Validate with NMDS ##
ord <- metaMDS(L2)
ord

plot(ord, type="n")
text(ord, display="sites", cex=0.7, col="blue")

plot(ord, type="n")
points(ord, display="sites", cex=0.5, col="black")
text(ord, display="species", cex=0.5, col="blue")

ordiplot(ord, type="n")
text(ord, display="sites", cex=0.7, col="blue")
ordihull(ord, groups = R2$Region, display="sites")
points(ord, display="species", cex=0.5, col="red")


## Anosim ##
L.dist <- vegdist(L2)
L.ano <- anosim(L.dist, R2$Region, permutations = 999)
summary(L.ano)

## Vector Fit ##
vf <- envfit(ord, R2, permutations = 999)

plot(ord, type="n")
text(ord, display="sites", cex=0.7, col="red")
plot(vf)


## 4th Corner ##
fcQ <- fourthcorner.rlq(rlq1, type="Q.axes") # test trait signficance to rlq axes
fcR <- fourthcorner.rlq(rlq1, type="R.axes") # test environment signficance to rlq axes

plot(fcQ, alpha=0.05, type="table", stat="D2")
plot(fcR, alpha=0.05, type="table", stat="D2")

plot(fcQ, alpha=0.05, type="biplot", stat="D2", col=c("black", "blue", "orange", "green"))
plot(fcR, alpha=0.05, type="biplot", stat="D2", col=c("black", "blue", "orange", "green"))