library(V.PhyloMaker)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)
library(phytools)
library(picante)
library(smatr)

#Read in Fan et al (2017) plant rooting database
infile = "C:/Phylo_GW_Study_120920/Plant_RD_Trees_List.csv"
RD = read.csv(infile,header=TRUE)
RD$FullName = paste(RD$Genus,RD$Species)
RD$MatchName = paste(RD$Genus,"_",RD$Species,sep="")

#Clear isotope column
RD$Isotope[RD$Isotope > 1] = 1

#Compute Root Depth / Water Table Depth
RD$RDWT_Fan = RD$RD_m/RD$WT_Fan_m
RD$RDWT_Bierke = RD$RD_m/RD$WT_Bierkie_m
RD = RD[is.finite(RD$RDWT_Fan),]

#Validate plant names against "The Plant List"
r1 <- TPL(RD$FullName, corr = TRUE)
TPL_Names = data.frame(species = r1$Taxon, genus = r1$Genus, family = r1$Family)
PlantNames = TPL_Names[which(TPL_Names$family != ""),]
PlantNames$matchname = gsub(" ", "_", PlantNames$species)

# Run V.Phyolomaker
result <- phylo.maker(PlantNames, scenarios=c("S1","S2","S3"))
result_analysis = result$scenario.3
#result_analysis$node.label <- NULL

#PICS/PLS requires one value per species
for (i in 1:length(result_analysis$tip.label))
{
  speciesname = result_analysis$tip.label[i]
  family = as.character(PlantNames[PlantNames$matchname == speciesname,3])
  Records = RD[paste(RD$Genus,"_",RD$Species,sep="") == speciesname,]
  RD = median(Records$RD_m)
  RDWT = median(Records$RDWT_Fan)
  AG = max(Records$Clade)
  WT = median(Records$WT_Fan_m)
  if (maxRD < 2.5) {result_analysis$trait.state_RD[i] = 0}
  if (maxRD >= 2.5) {result_analysis$trait.state_RD[i] = 1}
  result_analysis$family[i] = family[1]
  result_analysis$trait.RD[i] = RD
  result_analysis$trait.WT[i] = WT
  result_analysis$trait.RDWT[i] = RDWT
  result_analysis$trait.AG[i] = AG
  result_analysis$matchname[i] = paste(Records$Genus[1],"_",Records$Species[1],sep="")
}

#Tree to use for PICs analysis
tree<-multi2di(result$scenario.1)

#PICs
pic.RD = pic(result_analysis$trait.RD, tree)
pic.WT = pic(result_analysis$trait.WT, tree)
pic.RDWT = pic(result_analysis$trait.RDWT, tree)

PICS <- data.frame(pic.RD,pic.WT,pic.RDWT,result_analysis$trait.AG[1:159])
colnames(PICS) <- c("RD", "WT","RDWT","Clade")

#SMA regression of ln(RD) and ln(WT)
#Hypothesis test that slope of Gymnos = Angios
fit_SMA <- with(RD, slope.com(log(RD_m), log(WT_Fan_m), Clade, intercept=TRUE))
fit_SMA

#SMA regression of PIC (equivelant to PLS regression)
#Hypothesis test that slope of Gymnos = Angios
fit_PICS <- with(PICS, slope.com(RD, WT, Clade, intercept=TRUE))
fit_PICS

