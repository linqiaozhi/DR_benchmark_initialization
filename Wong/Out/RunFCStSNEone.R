#Evan Newell 2015, Ming 2016, Etienne Becht 2016, Evan 2017,2018

## PARAMETERS TO EDIT:
#####  to make sure python path is correct: here are some examples...
#Sys.setenv(PATH = paste("/Users/evan/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH = paste("/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH = paste("/Users/NewellEW/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
#Sys.setenv(PATH=paste("/Users/cadutertre/anaconda3/bin",Sys.getenv("PATH"),sep=":"))

#Choose functions to run
RunAlgorithms = T   # change to F if you just want to make meaning plots or heatplots using output files (no need to rerun tSNE etc.)
RunMeaningPlot = T
RunFileHeatplots =F
    #IF TRUE  ... In progress
    DoClusters=F # makes heatmaps to describe phenograph clusters - median marker intensity and frequencies (lin and log) within samples
    clustParam="Phenograph"  # e.g., "Phenograph" or "FlowSOM"
Run3DPlot = F
  

#Main Parameters
sourceFcsFolder = "WongTraffickGatedTNK" #Folder with samples (even 1 sample ok)
outputSuffix = "out-All"
useCSV = F   # T if using .csv files instead of fcs files - will make csv output als
CyTOF = T  # T for CyTOF logicle transform settings F for fluorescence flow
TransformOutput = F  # F Puts derived data on linear scale 1-10k - best for viewing in FJ10 otherwise T is setup for FJ9. Talk to Evan about implications for one-sense
ceil = 1000000
#number of events to take per sample (unless sample has less)
FNnames="trafficknames2.csv" #Parameters to analyze in csv file
OutputSuffix = outputSuffix # Adds this to all new directories and file names
DotSNE = F
tSNEperpelxity = 30 #default = 30; increased perplexity => increased spread
DoOneSENSE = F  #needs a modified names.csv file with column names for each category
Dophenograph = F #clusters cells using Rphenograpy. 
kValue = 20  #default = 30; k is low => more clusters and inversely
DoFlowSOM = F 
MaxClusters = 30
DoIsomap = F
DoDiffMap = F
DoPCA = F
Do3DtSNE = F #T for running 3D tSNE
DoUMAP =T
Do3DUMAP = T
DoOneSUMAP = F # umap version of One-SENSE 
  

## meaningplot paramters:
MeaningTopPercentile = 1
MeaningBotPercentile = 0
Xaxis = "UMAP1"   ## examples: tSNE1 or UMAP1  Add * to color by an additional round of tSNE or UMAP
Yaxis = "UMAP2"  ## examples: tSNE2 or UMAP1
DoOneForEach = F #T => will do meaning plots for each individual file, F => for the concat file
prefix = paste0(sourceFcsFolder,OutputSuffix)
palette=c("black","blue","green","yellow","red")
color.scale.type="relative"  # choose "relative" or "absolute"
resolution=150
cex=0.2
pch=16

# FileHeatplot (medians)
fileHPoutputsuffix = OutputSuffix
HPpalette=c("black","blue","lightblue","green","yellow","darkorange","darkred")



#3D plot options

Color3Dby = "InFile" ## eg "FlowSOM", "Phenograph", use "InFile" to color by sample file
parN1 = "TD_UMAP1"  #e.g., TD_UMAP1 or 3DtSNE1
parN2 = "TD_UMAP2" #e.g., TD_UMAP2 or 3DtSNE2
parN3 = "TD_UMAP3" #e.g., TD_UMAP3 or 3DtSNE3
TDOutputSuffix = "3D"
Ncolors = 0  # set to zero to have this automatically set to maximum number
library(RColorBrewer) # code for making random colors, good for clustering
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
paletteFor3D <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ColorList =c("red", "green", "blue") # can choose colors if numbers match otherwise it will interpolate
  # ColorList = c("colorX", "colorY", ...) put n color names if n clusters in order to manually define the colors of the clusters
  # Color palette in R: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf)





## No edit below here

source('FCStSNEone.R')
sourceFcsForMeaningPlot <- paste0(sourceFcsFolder,"_",outputSuffix)
if (RunAlgorithms) 
  FCStSNEone(LoaderPATH = sourceFcsFolder, 
           useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output als
           CyTOF = CyTOF, # T for CyTOF logicle transform settings F for fluorescence flow
           TransformOutput = TransformOutput, # F Puts derived data on linear scale 1-10k - best for viewing in FJ10 otherwise T is setup for FJ9. Talk to Evan about implications for one-sense
           ceil = ceil, #number of events to take per sample (unless sample has less)
           FNnames=FNnames,#Parameters to analyze in csv file
           OutputSuffix = outputSuffix, # Adds this to all new directories and file names
           DotSNE = DotSNE,
           tSNEperpelxity = tSNEperpelxity, #default = 30; increased perplexity => increased spread
           DoOneSENSE = DoOneSENSE, #needs a modified names.csv file with column names for each category
           Dophenograph = Dophenograph, #clusters cells using Rphenograpy
           kValue = kValue, #default = 30; k is low => more clusters and inversely
           DoFlowSOM = DoFlowSOM,
           MaxClusters = MaxClusters,
           DoIsomap = DoIsomap,
           DoDiffMap = DoDiffMap,
           DoPCA = DoPCA, 
           Do3DtSNE = Do3DtSNE, #T for running 3D tSNE
           DoUMAP = DoUMAP, #still in prep- will add many parameters for this
           Do3DUMAP = Do3DUMAP,
           DoOneSUMAP = DoOneSUMAP) 


## Still in prep:
if (RunMeaningPlot)
   meaningPlot(LoaderPATH =sourceFcsForMeaningPlot,
            useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output als
            CyTOF = CyTOF,
            FNnames = FNnames,
            ceil = ceil,
            TransformOutput = TransformOutput,
            MeaningTopPercentile = .99,
            MeaningBotPercentile = .01,
            PC1 = Xaxis,
            PC2 = Yaxis,
            DoOneForEach = DoOneForEach,
            prefix = prefix,
            palette=palette,
            color.scale.type=color.scale.type,
            resolution=resolution,
            cex=cex,
            pch=pch)

if(RunFileHeatplots)
  fileHeatplot(LoaderPATH =sourceFcsForMeaningPlot,
                          useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output als
                          CyTOF = CyTOF,
                          ceil = ceil,
                          TransformOutput = TransformOutput,
                          FNnames = FNnames,
                          OutputSuffix = fileHPoutputsuffix,
                          palette=HPpalette,
                         DoClusters=DoClusters,
                         clustParam=clustParam)

if(Run3DPlot)
  ThreeDPlot (LoaderPATH =sourceFcsForMeaningPlot,
              CyTOF = CyTOF,
              FNames = FNames,
              ceil = ceil, 
              OutputSuffix = TDOutputSuffix,
              parN1 = parN1,
              parN2 = parN2, 
              parN3= parN3, 
              Ncolors = Ncolors,
              colparN = Color3Dby,
              ColorList = ColorList)  # can use "InFile" to color by sample file

##### Useful script for making names.csv file:
# if (!require(flowCore)) { 
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("flowCore")
# } 
# 
# inFileN <- "test.fcs"        #FCS file to make names file for 
# outNamesCsv <- "names.csv"   #Output csv file name 
# 
# FF <- read.FCS(inFileN)
# colNames <- FF@parameters$desc
# empties <- which(is.na(colNames) | colNames== " ")
# colNames[empties] <- FF@parameters$name[empties]
# write.csv(colNames, outNamesCsv, row.names = F)