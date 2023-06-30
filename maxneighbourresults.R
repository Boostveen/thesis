setwd("~/Master IS/Master Thesis/Onderzoeks fase (huidig)/R/R 22-05-2023 Max neighbour effect/pretty r")

library(dplyr)# for transforming data
library(car)  # for VIF test
library(lmtest) # for BP test
library(spdep) #for creating links between the regions
library(sfdep)# for further spaial datashit
library(rgdal) #for importing shapedata
library(RColorBrewer) # Kleurtjes voor model
library(spatialreg) #spatial regressions
library(stargazer) # tabelletjes
library(vtable) # sumtable
library(corrplot)  # Correlation
library(matrixStats) # matrix stats for max lagged values

data <- read.delim("~/Master IS/Master Thesis/Onderzoeks fase (huidig)/R/R 22-05-2023 Max neighbour effect/data4v14maxlag.csv", sep =',')
sumtable(data)



#shapefile import ####
file_path <- "~/Master IS/Master Thesis/Onderzoeks fase (huidig)/R/Latest data versions (keep updating)/shapefile good" 
nuts2_shapefile <- readOGR(dsn = file_path , layer = "shapefile")
#create nb file and nb list
nb <- poly2nb(nuts2_shapefile) # maak buurman bestand van de shapefile
# connecting DK01(151) & SE22 (254)
nb[[151]] <- sort(as.integer(c(nb[[151]], 254)))
nb[[254]] <- sort(as.integer(c(151, nb[[254]])))
# connecting DK02(159) & DK03 (174)
nb[[159]] <- sort(as.integer(c(nb[[159]], 174)))
nb[[174]] <- sort(as.integer(c(159, nb[[174]])))
lw <- nb2listw(nb) # neighbour weight matrix maken
nuts2_shapefile@data <- inner_join(nuts2_shapefile@data, data, by = "NUTS_ID")

spat.data1 <- nuts2_shapefile

# demand Run separate regression models for each category
cv <- lm(Crunch17_21 ~ internal_structure, data = spat.data1)
dem <- lm(Crunch17_21 ~ internal_structure  + Demand , data = spat.data1)
demmax <- lm(Crunch17_21 ~ internal_structure  + Demand + lagmaxDemand, data = spat.data1)
stargazer(cv, dem, demmax, type = "text", out = "Demandmax.htm")
vif(dem)
vif(demmax)
bptest(demmax)
bptest(dem)

moran.test(spat.data1$Demand, lw, 500)
moran.mc(spat.data1$Demand, lw, 500)
moran.mc(dem$residuals, lw, 500)
moran.plot(dem$residuals,lw)
moran.plot(spat.data1$Crunch17_21, lw)
plot(spat.data1$Crunch17_21, spat.data1$lagmaxDemand )
plot(demmax$residuals, spat.data1$Crunch17_21)

# Residual Dem model
residualsdem <- lm(formula = internal_structure ~ Demand, data = spat.data1@data)
cvdem <- lm(Crunch17_21 ~ residualsdem$residuals, data = spat.data1@data)
dem2 <- lm(Crunch17_21 ~ residualsdem$residuals  + Demand , data = spat.data1)
demmax2 <- lm(Crunch17_21 ~ residualsdem$residuals  + Demand + lagmaxDemand, data = spat.data1)
stargazer(cvdem, dem2, demmax2, type = "text", out = "Demandmaxresmodel.htm")
bptest(cv)

# Talent
tal <- lm(Crunch17_21 ~ internal_structure  + Talent , data = spat.data1)
talmax <- lm(Crunch17_21 ~ internal_structure  + Talent + lagmaxTalent, data = spat.data1)
stargazer(cv, tal, talmax, type = "text", out = "Talentmax.htm")
vif(tal)
vif(talmax)
bptest(talmax)
bptest(tal)
moran.test(spat.data1$Talent, lw, 500)

# intermediates
int <- lm(Crunch17_21 ~ internal_structure  + Intermediate , data = spat.data1)
intmax <- lm(Crunch17_21 ~ internal_structure  + Intermediate + lagmaxIntermediate, data = spat.data1)
stargazer(cv, int, intmax, type = "text", out = "Intermediatemax.htm")
vif(int)
vif(intmax)
bptest(intmax)
bptest(int)
moran.test(spat.data1$Intermediate, lw, 500)
moran.test(int$residuals, lw, 500)

# knowledge
kno <- lm(Crunch17_21 ~ internal_structure  + Knowledge2 , data = spat.data1)
knomax <- lm(Crunch17_21 ~ internal_structure  + Knowledge2 + lagmaxKnowledge2, data = spat.data1)
stargazer(cv, kno, knomax, type = "text", out = "Knowledge2max.htm")
vif(kno)
vif(knomax)
bptest(knomax)
bptest(kno)
moran.test(spat.data1$Knowledge2, lw, 500)

# leadership
lea <- lm(Crunch17_21 ~ internal_structure  + Leadership , data = spat.data1)
leamax <- lm(Crunch17_21 ~ internal_structure  + Leadership + lagmaxLeadership, data = spat.data1)
stargazer(cv, lea, leamax, type = "text", out = "Leadershipmax.htm")
vif(leamax)
bptest(leamax)
bptest(lea)
moran.test(spat.data1$Leadership, lw, 500)

# finance
fin <- lm(Crunch17_21 ~ internal_structure  + Finance , data = spat.data1)
finmax <- lm(Crunch17_21 ~ internal_structure  + Finance + lagmaxFinance, data = spat.data1)
stargazer(cv, fin, finmax, type = "text", out = "Financemax.htm")
vif(finmax)
bptest(finmax)
bptest(fin)
moran.test(spat.data1$Finance, lw, 500)

# networks
net <- lm(Crunch17_21 ~ internal_structure  + Networks2 , data = spat.data1)
netmax <- lm(Crunch17_21 ~ internal_structure  + Networks2 + lagmaxNetworks2, data = spat.data1)
stargazer(cv, net, netmax, type = "text", out = "Networks2max.htm")
vif(net)
vif(netmax)
bptest(netmax)
moran.test(spat.data1$Networks2, lw, 500)
bptest(net)

# overall tests
moran.test(spat.data1$Crunch17_21, lw, 500)

sumtable(data)

# show descriptives maxes
data2 <- data[,c(15:21)]
sumtable(data2)

#view spatial relation
cor(data$lagmaxDemand, data$Demand)
cor(data$Demand, data$Crunch17_21)

moran.plot(nuts2_shapefile@data$Crunch17_21, lw)

#correlation table creation 

coreverything <- data[,c(2,6:12,14:21)]
cordatamatrix <- cor(coreverything)
cordatamatrixp <- cor.mtest(coreverything, conf.level = 0.95)

cordatamatrix2 <- round(cordatamatrix, 2)


corrplot(round(cordatamatrix,2), p.mat = cordatamatrixp$p, method = 'color', diag = F,
         type = 'lower',tl.cex=0.5, number.cex=0.5, insig='blank',
         addCoef.col ='black')

cordatamatrix2 <- data.frame(cordatamatrix2)
write.csv(cordatamatrix2, file = "cordatamatrix3.csv")

#data sumtable
data1 <- data[,c(2,6:12,14)]
sumtable(data1)


# Create flinke output tabel lol.
all <- lm(Crunch17_21 ~ internal_structure  + Networks2 + lagmaxNetworks2 + Finance, data = spat.data1)
summary(all)
stargazer(cv, all, all, all, all, all, all, all,  type = "text", out ="all.htm")
