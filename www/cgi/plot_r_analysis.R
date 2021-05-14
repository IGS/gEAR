#!/usr/bin/env Rscript
#
# THE PRE-PROCESSED ELEMENTS THAT HAVE ALREADY BEEN PRODUCED ARE IN:
# /local/projects-t2/idea/projectoR_output
# -rw-r--r--. 1 ccolantuoni igs   1229888 Sep 17 23:30 BrainSpanDFC_Hclust.RData
# -rw-r--r--. 1 ccolantuoni igs   3789272 Aug  9 14:36 BrainSpanDFC_Kmeans.RData
# -rw-r--r--. 1 ccolantuoni igs  13810223 Aug  9 14:16 BrainSpanDFC_pca.RData
# -rw-r--r--. 1 ccolantuoni igs   1341017 Aug  9 14:10 POU3F2oe_Hclust.RData
# -rw-r--r--. 1 ccolantuoni igs   1542291 Aug  9 13:47 POU3F2oe_Kmeans.RData
# -rw-r--r--. 1 ccolantuoni igs   8602446 Aug  9 13:44 POU3F2oe_pca.RData
#
# SO FOR EACH OF OUR 2 EXAMPLE DATA SETS, HCLUST, KMEANS, AND PCA HAVE BEEN CALCULATED.
#
# DIMENSIONS FOR Y-VALUE PLOTTING ARE BURIED IN EACH OF THESE OBJECTS DIFFERENTLY - I guide you thru this below for each pre-processing type:
###########################################################################

args <- commandArgs(TRUE)

# dataset label (and file name base)
lb01 = args[1]

#data.path = "/usr/local/projects/gEAR/carlo/"
data.path = args[2]

# analysis_method
ana = args[3]

image.base = args[4]
cluster.count = args[5]

load(file=paste(data.path, lb01, "_data.RData", sep=""))
load(file=paste(data.path, lb01, "_dataROWanno.RData", sep=""))
load(file=paste(data.path, lb01, "_dataCOLanno.RData", sep=""))

##############################################
# PCA file load and parameter handling
##############################################
if (ana == "pca") {
    # PCA contians as many dimensions as there are samples in the data, so the user could request between 1 and the number of samples.
    load(file=paste(data.path, lb01, "_pca.RData", sep=""))

    # each column in pca$x is the Y-values for PCA plotting - so we would just plot those requested by the user
    # str(pca$x)
    
    KK = cluster.count  # users asked for N PCs
}

##############################################
# Kmeans 
##############################################
if (ana == "kmeans") {
    # Kmeans is run for pre-set (we could make this user specified during data import) numbers of clusters.
    load(file=paste(data.path, lb01, "_Kmeans.RData", sep=""))
    
    # in this example, sets of 5,10,15,20,25,30, and 35 clusters were defined:
    names(KmeansList)
    
    # [1] "K5"  "K10" "K15" "K20" "K25" "K30" "K35"
    KK = cluster.count  # users asked for 5 clusters
    
    # this will produce an error if the user requests a KK that has not been precomputed
    indx=match(KK, gsub("K", "", names(KmeansList)))

    # str(KmeansList[[indx]])# the "$cluster" element of this object contains gene cluster assignments
    cutp=KmeansList[[indx]]$cluster

    # cutp is now a vector of integers defining which of the 5 clusters each gene of this data sets is in.
    # code for calculating average clusters is below - see "Average Clusters"
    # code for plotting the average clusters is below - see "Plotting Code"
}

##############################################
# Hclust 
##############################################
if (ana == "hclust") {
    # A USER CAN ASK FOR AS MANY DIMENSIONS AS THEY WANT OUT OF THE Hclust OBJECTS:
    # SAY THE USERS ASKS FOR 5 DIMENSIONS/CLUSTERS FROM THE Hclust OF 'BrainSpanDFC' - IN R WE WOULD USE THIS CODE TO PLOT THE 5 CLUSTERS THE USER IS ASKING FOR:
    load(file=paste(data.path, lb01, "_Hclust.RData", sep=""))
    
    # the pre-processed 'Hclust' files contain only 1 object: 'clust1' that defines the clustering for all #s
    KK = cluster.count  # users asked for 5 clusters
    cutp=cutree(clust1,k=KK)  # extract 5 clusters
    
    # cutp is now a vector of integers defining which of the 5 clusters each gene of this data sets is in. This code will plot the 5 clusters:
    # code for calculating average clusters is below - see "Average Clusters"
    # code for plotting the average clusters is below - see "Plotting Code"
}

##############################################
# "Average Clusters" # Used for Hclust AND Kmeans, NOT pca
##############################################
if (ana == "kmeans" | ana == "hclust") {
    # code for looking at average cluster values and how tight the clusters are
    cls = sort(unique(cutp))
    cMNs = matrix(ncol = dim(data)[2],
                  nrow = length(cls))
    
    meanRRs=vector(length = length(cls))
    
    for(i in cls) {
        if (sum(cutp == i) > 1) {
            exprsMTXc=data[cutp==i,]
            exprsMTXcMN=colMeans(exprsMTXc)
            cMNs[i,]=exprsMTXcMN
            meanRRs[i]=mean(apply(exprsMTXc,1,cor,y=exprsMTXcMN))
            
        } else if (sum(cutp == i) == 1) {
            cMNs[i, ] = data[cutp == i,]
            meanRRs[i]=1
            
        } else if (sum(cutp == i) == 0) {
            print("cluster error !")
        }
    }
    
    # these warnings are OK:
    # 1: In FUN(newX[, i], ...) : the standard deviation is zero
    # 2: In FUN(newX[, i], ...) : the standard deviation is zero
}

##############################################
# "Plotting Code"
##############################################
# for pca we have "pca$x" for Yvalues
if (ana == "pca") {
    YY = pca$x
    Ylb = "PC "
    Ylb2 = pcvars
    xtrlb = "% variance"

# For both Kmeans and Hclust which have produced "cutp" objects above that together
# with the expression data give us what we need to plot - done above in "Average Clusters".
} else if (ana == "kmeans") {
    YY = t(cMNs)
    Ylb = "Kmeans "
    Ylb2 = round(meanRRs,2)
    xtrlb = "= within cluster r"
    
} else if (ana == "hclust") {
    YY = t(cMNs)
    Ylb = "Hclust "
    Ylb2 = round(meanRRs,2)
    xtrlb = "= within cluster r"
}

# the "data" object used below contains the expression data and was loaded above already
# PLOT
for(i in 1:KK) {
    image.path = paste(image.base, i, ".png", sep="")
    
    png(width = 780,
        height = 480,
        file = image.path)
    
    plot(x = dataCOLanno[,X0],
         y = YY[,i],
         col = dataCOLanno[,clr0],
         pch = dataCOLanno[,pch0],
         cex = 2,
         main = paste(lb0, ": ", ana," (from ", YexprsLB0, ")", sep=""),
         xlab = xlb0,
         ylab = paste(Ylb, i, " :  ", Ylb2[i], xtrlb, sep=" "))
    
    if (ylo0) {
        lines(x=dataCOLanno[,X0],
              y=loess(YY[,i]~dataCOLanno[,X0],spn=spn0)$fit)
    }
    
    abline(v   = abl0,
           lty = 2)
    
    text(x = abl0,
         y = par()$yaxp[1],labels=ablb0)
    
    dev.off()
}

# objects loaded from "_dataCOLanno" file, and used in above plotting code:
	# lb0 - data set name/identifier
	# dataCOLanno - column/sample annotation
	# colID - column name in "dataCOLanno" where column names (sample names) of data are
	# xlb0 - X-axis label
	# ylo0 - T/F - do you want to plot a local mean of Y-values across X?
	# spn0 - smoothing parameter for local mean of Y-values across X, if requested in ylo0
	# abl0 - vertical dotted lines to be displayed at these X-values
	# ablb0 - labels for abl0 dotted lines
	# X0 - column in dataCOLanno where X-values should come from
	# clr0 - column in dataCOLanno where plotting colors should come from
	# pch0 - column in dataCOLanno where plotting symbols should come from
	# YexprsLB0 - Data type label
# objects loaded from "_dataROWanno" file, and used in above plotting code:
	# rowID - column name in "dataROWanno" where row names (gene/locus names) of data are
	# ROWplotID -  column name in "dataROWanno" where gene IDs

# i set all these defualts when i created the _dataCOLanno and _dataROWanno files along with the primary _data file
# users should set these values when they import data for the 1st time
