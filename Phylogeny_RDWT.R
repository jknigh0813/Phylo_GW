library(V.PhyloMaker)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)
library(phytools)
library(picante)

#Read in Fan et al (2017) plant rooting database
infile = "E:/Global_Rooting/Plant_RD_TreesShrubs_List.csv"
#infile = "E:/Global_Rooting/Plant_RD_All_List.csv"
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

# Run Phyolomaker
result <- phylo.maker(PlantNames, scenarios=c("S1","S2","S3"))
result_analysis = result$scenario.3
#result_analysis$node.label <- NULL

#Assign RDWT and Iso as traits
for (i in 1:length(result_analysis$tip.label))
  {
    speciesname = result_analysis$tip.label[i]
    family = as.character(PlantNames[PlantNames$matchname == speciesname,3])
    Records = RD[paste(RD$Genus,"_",RD$Species,sep="") == speciesname,]
    rootdepth = median(Records$RD_m)
    maxRDWT_Fan = median(Records$RDWT_Fan)
    maxRDWT_Bierke = median(Records$RDWT_Bierke)
    if (maxRDWT_Fan < 1) {result_analysis$trait.state_RD_Fan[i] = 0}
    if (maxRDWT_Fan == 1) {result_analysis$trait.state_RD_Fan[i] = 1}
    if (maxRDWT_Fan > 1) {result_analysis$trait.state_RD_Fan[i] = 1}
    if (maxRDWT_Bierke < 1) {result_analysis$trait.state_RD_Bierke[i] = 0}
    if (maxRDWT_Bierke == 1) {result_analysis$trait.state_RD_Bierke[i] = 1}
    if (maxRDWT_Bierke > 1) {result_analysis$trait.state_RD_Bierke[i] = 1}
    if (rootdepth < 2) {result_analysis$trait.state_RD[i] = 0}
    if (rootdepth >= 2) {result_analysis$trait.state_RD[i] = 1}
    maxIso = round(max(Records$Isotope))    
    result_analysis$family[i] = family[1]
    result_analysis$trait.RD[i] = rootdepth
    result_analysis$trait.RDWT_Fan[i] = maxRDWT_Fan
    result_analysis$trait.RDWT_Bierke[i] = maxRDWT_Bierke
    result_analysis$trait.state_Iso[i] = maxIso
    result_analysis$matchname[i] = paste(Records$Genus[1],"_",Records$Species[1],sep="")
  }


#Reformat traits
Traits = data.frame(result_analysis$family, result_analysis$trait.state_RD_Bierke, result_analysis$trait.state_RD_Fan, result_analysis$trait.state_RD)
row.names(Traits) <- result_analysis$tip.label
colnames(Traits) <- c("Family","RDWT_Bierke","RDWT_Fan","RD")
write.csv(Traits,"RDWT_Fan.csv")

#K and Lambda tests for trait occurence
L1 = phylosig(result_analysis, result_analysis$trait.RDWT_Bierke, method="lambda",test=TRUE)
L2 = phylosig(result_analysis, result_analysis$trait.RDWT_Fan, method="lambda",test=TRUE)
L3 = phylosig(result_analysis, result_analysis$trait.RD, method="lambda",test=TRUE)

K1 = phylosig(result_analysis, result_analysis$trait.RDWT_Bierke, method="K",test=TRUE)
K2 = phylosig(result_analysis, result_analysis$trait.RDWT_Fan, method="K",test=TRUE)
K3 = phylosig(result_analysis, result_analysis$trait.RD, method="K",test=TRUE)


#Phylogeny plot
tiff("E:/Global_Rooting/Figures/Phylogeny_RDWT_TreesShrubs_031820.tiff", units="in", width=5, height=5, res=600)
trait.plot(result_analysis, dat=Traits, class=result_analysis$family, cols = list(RDWT_Bierke = c("lightblue1","blue4"),RDWT_Fan = c("pink", "red")),cex.lab=0.4,cex.legend=0.5)
text(x=0, y=150, paste("n = ",length(result_analysis$tip.label)," species",sep=""),col ="black", cex=0.7)
text(x=0, y=120, paste("??p = ",round(L2$lambda,3), " (p = ",round(L2$P,3),")",sep=""),col ="red", cex=0.7)
text(x=0, y=90, paste("Kb = ",round(K2$K,3), " (p = ",round(K2$P,3),")",sep=""),col ="red", cex=0.7)
text(x=0, y=60, paste("??p = ",round(L1$lambda,3), " (p = ",round(L1$P,3),")",sep=""),col ="blue", cex=0.7)
text(x=0, y=30, paste("Kb = ",round(K1$K,3), ", (p = ",round(K1$P,3),")",sep=""),col ="blue", cex=0.7)
dev.off()

#Export traits for barplots
infile = "E:/Global_Rooting/Plant_RD_TreesShrubs_List.csv"
RD = read.csv(infile,header=TRUE)
RD$FullName = paste(RD$Genus,RD$Species)
RD$MatchName = paste(RD$Genus,"_",RD$Species,sep="")
RD$RDWT_Fan = RD$RD_m/RD$WT_Fan_m
RD$RDWT_Bierke = RD$RD_m/RD$WT_Bierkie_m

Export_RD = data.frame(RD$MatchName, RD$RD_m,RD$RDWT_Fan,RD$RDWT_Bierke,"A",stringsAsFactors=FALSE)
colnames(Export_RD) <- c("FullName","RD_m","RDWT_F","RDWT_B","Family")
Export_RD[is.na(Export_RD)] <- -999

for (i in 1:length(Export_RD[,1]))
{
  family = as.character(PlantNames[PlantNames$matchname == Export_RD$FullName[i],3])
  Export_RD$Family[i] = family[1]
}
Export_RD <- Export_RD[order(Export_RD$Family),]

write.csv(Export_RD,"RootTraits_Family.csv")
