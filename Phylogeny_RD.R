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

#Validate plant names against "The Plant List"
r1 <- TPL(RD$FullName, corr = TRUE)
TPL_Names = data.frame(species = r1$Taxon, genus = r1$Genus, family = r1$Family)
PlantNames = TPL_Names[which(TPL_Names$family != ""),]
PlantNames$matchname = gsub(" ", "_", PlantNames$species)

# Run Phyolomaker
result <- phylo.maker(PlantNames, scenarios=c("S1","S2","S3"))
result_analysis = result$scenario.3
result_analysis$node.label <- NULL

#Assign max RD as a trait
for (i in 1:length(result_analysis$tip.label))
  {
    speciesname = result_analysis$tip.label[i]
    family = as.character(PlantNames[PlantNames$matchname == speciesname,3])
    Records = RD[paste(RD$Genus,"_",RD$Species,sep="") == speciesname,]
    maxRD = max(Records$RD_m)
    if (maxRD < 2.5) {result_analysis$trait.state_RD[i] = 0}
    if (maxRD >= 2.5) {result_analysis$trait.state_RD[i] = 1}
    result_analysis$family[i] = family[1]
    result_analysis$trait.RD[i] = maxRD
    result_analysis$matchname[i] = paste(Records$Genus[1],"_",Records$Species[1],sep="")
  }

#Reformat traits
Traits = data.frame(result_analysis$trait.state_RD, result_analysis$trait.RD)
row.names(Traits) <- result_analysis$tip.label
colnames(Traits) <- c("RD")

#Run K and lambda tests for phylogenetic signal
K1 = phylosig(result_analysis, result_analysis$trait.RD, method="K",test=TRUE,nsim=1000, se=NULL, start=NULL,
         control=list())
L1 = phylosig(result_analysis, result_analysis$trait.RD, method="lambda",test=TRUE,nsim=1000, se=NULL, start=NULL,
         control=list())

#Phylogeny plot
tiff("E:/Global_Rooting/Figures/Phylogeny_RD_TreesShrubs_031820.tiff", units="in", width=5, height=5, res=600)
trait.plot(result_analysis, dat=Traits, class=result_analysis$family, cols = list(RD = c("pink","red")),cex.lab=0.4)
text(x=0, y=90, paste("n = ",length(result_analysis$tip.label)," species",sep=""),col ="black", cex=0.8)
text(x=0, y=60, paste("??p = ",round(L1$lambda,3), " (p = ",round(L1$P,3),")",sep=""),col ="red", cex=0.8)
text(x=0, y=30, paste("Kb = ",round(K1$K,3), " (p = ",round(K1$P,3),")",sep=""),col ="red", cex=0.8)
dev.off()


