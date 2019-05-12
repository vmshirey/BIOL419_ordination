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

## Drop transects with no abundance data ##
dropNames <- row.names(L2[rowSums(L2[]) == 0, ]) # get row names from 0 transect data
R2 <- R2[!row.names(R2) %in% dropNames, ]
L2 <- L2[!row.names(L2) %in% dropNames, ]

#############################################################################
## RLQ                                                                     ##
#############################################################################

coa1 <- dudi.coa(L2, scannf=FALSE, nf=2)
dudimil <- dudi.hillsmith(R2, scannf=FALSE, nf=2, row.w = coa1$lw)
duditrait <- dudi.hillsmith(Q2, scannf=FALSE, nf=2, row.w = coa1$cw)

rlq1 <- rlq(dudimil, coa1, duditrait, scannf=FALSE, nf=2) # conduct RLQ

plot(rlq1)

summary(rlq1)
rlq_rand <- randtest(rlq1)

## Final Plot for RLQ: Elevation Class ##
xlim <- range(rlq1$lR[,1])
ylim <- range(rlq1$lR[,2])
colVec <- c("black", "red")

plot.new()
plot.window(xlim=c(-5,5), ylim=c(-6,2), asp=0.5)

abline(h = 0, lty="dotted")
abline(v = 0, lty="dotted")

with(R2, points(rlq1$lR, col=colVec[ElevationClass], pch=21, bg=colVec[ElevationClass]), cex=0.5)
with(R2, legend("topright", legend=levels(ElevationClass), bty="n", col=colVec, 
                pch=21, pt.bg = colVec))
legend("bottomleft", legend=c("Elevation", "Broad Oligophagy", "Multiple Diapause"),
       bty="n", col=c("black", "red", "blue"), lty=1)

arrows(0,0,rlq1$li[3,1], rlq1$li[3,2])
arrows(0,0,rlq1$co[4,1], rlq1$co[4,2], col="red")
arrows(0,0,rlq1$co[11,1], rlq1$co[11,2], col="blue")
text(-5,6, "Trait/Env. Factors Not Signficant")

axis(side = 1)
axis(side = 2)
title(xlab = "S1", ylab="S2", main="RLQ for Cascadian Butterfly Transects By Elevation Class")
box()

## Final Plot for RLQ: Habitat ##
xlim <- range(rlq1$lR[,1])
ylim <- range(rlq1$lR[,2])
colVec <- c("black", "blue", "red")

plot.new()
plot.window(xlim=c(-5,5), ylim=c(-6,6), asp=0.5)

abline(h = 0, lty="dotted")
abline(v = 0, lty="dotted")

with(R2, points(rlq1$lR, col=colVec[Habitat], pch=21, bg=colVec[Habitat]), cex=0.5)
with(R2, legend("topright", legend=levels(Habitat), bty="n", col=colVec, 
                pch=21, pt.bg = colVec))
legend("bottomleft", legend=c("Elevation", "Broad Oligophagy", "Multiple Diapause"),
       bty="n", col=c("black", "red", "blue"), lty=1)

arrows(0,0,rlq1$li[3,1], rlq1$li[3,2])
arrows(0,0,rlq1$co[4,1], rlq1$co[4,2], col="red")
arrows(0,0,rlq1$co[11,1], rlq1$co[11,2], col="blue")
text(-5,6, "Trait/Env. Factors Not Signficant")

axis(side = 1)
axis(side = 2)
title(xlab = "S1", ylab="S2", main="RLQ for Cascadian Butterfly Transects By Habitat")
box()

##################################################################################
## NMDS                                                                         ##
##################################################################################

ord <- metaMDS(L2)
ord

## ANOSIM: Elevation Class  ##
L.dist <- vegdist(L2)
L.ano.elev <- anosim(L.dist, R2$ElevationClass, permutations = 999)
summary(L.ano.elev)

## ANOSIM: Habitat  ##
L.dist <- vegdist(L2)
L.ano.hab <- anosim(L.dist, R2$Habitat, permutations = 999)
summary(L.ano.hab)

## Vector-fit: Elevation ##
vf <- envfit(ord, R2[,c(2,4:5)])
vf

par(mfrow=c(1,1))

## Final Plot for NMDS: Elevation Class ##
scrs <- scores(ord, display=c("sites", "species"), scaling=scl)
colVec <- c("black", "red")
xlim <- range(scrs[,1])
ylim <- range(scrs[,2])

plot.new()
plot.window(xlim=xlim, ylim=ylim, asp=1)

abline(h = 0, lty="dotted")
abline(v = 0, lty="dotted")

with(R2, points(scrs, col=colVec[ElevationClass], pch=21, bg=colVec[ElevationClass]))

with(R2, legend("topright", legend=levels(ElevationClass), bty="n", col=colVec, 
                     pch=21, pt.bg = colVec))

with(R2, ordiellipse(ord, ElevationClass, col=colVec))
plot(vf, p.max=0.1)
text(2,-1.5, paste("ANOSIM p-value: ", L.ano.elev$signif))

axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab="NMDS 2", main="NMDS Plot for Cascadian Butterfly Transects By Elevation Class")
box()

## Final Plot for NMDS: Habitat ##
scrs <- scores(ord, display=c("sites", "species"), scaling=scl)
colVec <- c("black", "blue", "red")
xlim <- range(scrs[,1])
ylim <- range(scrs[,2])

plot.new()
plot.window(xlim=xlim, ylim=ylim, asp=1)

abline(h = 0, lty="dotted")
abline(v = 0, lty="dotted")

with(R2, points(scrs, col=colVec[Habitat], pch=21, bg=colVec[Habitat]))

with(R2, legend("topright", legend=levels(Habitat), bty="n", col=colVec, 
                     pch=21, pt.bg = colVec))

with(R2, ordiellipse(ord, Habitat, col=colVec))
plot(vf, p.max=0.1)
text(2,-1.5, paste("ANOSIM p-value: ", L.ano.hab$signif))

axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab="NMDS 2", main="NMDS Plot for Cascadian Butterfly Transects By Habitat")
box()

################################################################################
## 4th Corner                                                                 ##
################################################################################

fcQ <- fourthcorner.rlq(rlq1, type="Q.axes") # test trait signficance to rlq axes
fcR <- fourthcorner.rlq(rlq1, type="R.axes") # test environment signficance to rlq axes

plot(fcQ, alpha=0.05, type="table", stat="D2")
plot(fcR, alpha=0.05, type="table", stat="D2")

plot(fcQ, alpha=0.05, type="biplot", stat="D2", col=c("black", "blue", "orange", "green"))
plot(fcR, alpha=0.05, type="biplot", stat="D2", col=c("black", "blue", "orange", "green"))