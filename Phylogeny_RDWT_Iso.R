library(V.PhyloMaker)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)
library(phytools)

#Read in Fan et al (2017) plant rooting database
infile = "E:/Global_Rooting/Plant_RD_TreesShrubs_List.csv"
RD = read.csv(infile,header=TRUE)
RD$FullName = paste(RD$Genus,RD$Species)

#Clear isotope column
RD$Isotope[RD$Isotope > 1] = 1

#Compute Root Depth / Water Table Depth
RD$RDWT_Fan = RD$RD_m/RD$WT_Fan_m
RD$RDWT_Bierke = RD$RD_m/RD$WT_Bierkie_m

#Remove very shallow systems
#RD = RD[RD$RD > 0.2,]
#RD = RD[RD$WT > 0.2,]
RD = RD[is.finite(RD$WT_Fan_m),]

#Validate plant names against "The Plant List"
r1 <- TPL(RD$FullName, corr = TRUE)
TPL_Names = data.frame(species = r1$Taxon, genus = r1$Genus, family = r1$Family)
PlantNames = TPL_Names[which(TPL_Names$family != ""),]

# Run Phyolomaker
result <- phylo.maker(PlantNames, scenarios=c("S1","S2","S3"))
result_analysis = result$scenario.3
result_analysis$node.label <- NULL

#Assign RDWT and Iso as traits
for (i in 1:length(result_analysis$tip.label))
  {
    speciesname = result_analysis$tip.label[i]
    Records = RD[paste(RD$Genus,"_",RD$Species,sep="") == speciesname,]
    maxRDWT_Fan = median(Records$RDWT_Fan)
    maxRDWT_Bierke = median(Records$RDWT_Bierke)
    #Classify RD/WT as <0, =1, or >1
    if (maxRDWT_Fan < 1) {result_analysis$trait.state_RD_Fan[i] = 0}
    if (maxRDWT_Fan == 1) {result_analysis$trait.state_RD_Fan[i] = 1}
    if (maxRDWT_Fan > 1) {result_analysis$trait.state_RD_Fan[i] = 2}
    if (maxRDWT_Bierke < 1) {result_analysis$trait.state_RD_Bierke[i] = 0}
    if (maxRDWT_Bierke == 1) {result_analysis$trait.state_RD_Bierke[i] = 1}
    if (maxRDWT_Bierke > 1) {result_analysis$trait.state_RD_Bierke[i] = 2}
    maxIso = round(max(Records$Isotope))    
    result_analysis$trait.RDWT_Fan[i] = maxRDWT_Fan
    result_analysis$trait.RDWT_Bierke[i] = maxRDWT_Bierke
    result_analysis$trait.state_Iso[i] = maxIso
    result_analysis$matchname[i] = paste(Records$Genus[1],"_",Records$Species[1],sep="")
  }

#D-test for trait occurrence
Traits_Dtest = data.frame(result_analysis$matchname, result_analysis$trait.state_Iso)
colnames(Traits_Dtest) <- c("matchname","Iso")
Trees = comparative.data(result_analysis, Traits_Dtest, matchname)
redPhyloD <- phylo.d(Trees, binvar=Iso)
print(redPhyloD)

#Reformat traits for plotting
Traits = data.frame(result_analysis$trait.state_Iso)
row.names(Traits) <- result_analysis$tip.label
colnames(Traits) <- c("Iso")

#Phylogeny plot
tiff("C:/Phylo_GW_Study_120920/Figures/Phylogeny_ISO_Honly_excluded_031820.tiff", units="in", width=5, height=5, res=600)
trait.plot(result_analysis, dat=Traits, class=result_analysis$family, cols = list(Iso = c("lightblue1","blue4")),cex.lab=0.4)
text(x=0, y=60, paste("n = ",length(result_analysis$tip.label)," species",sep=""),col ="black", cex=0.8)
text(x=0, y=30, paste("D = ",round(redPhyloD$DEstimate,3), " (p =",round(redPhyloD$Pval1,3),")",sep=""),col ="blue", cex=0.8)
dev.off()

