#!/usr/bin/Rscript

library("optparse")
 
option_list = list(
    make_option(c("-ps", "--projection-source"), type="character", default=NULL, 
                help="projection source", metavar="character"),
    make_option(c("-sop", "--set-of-patterns"), type="character", default=NULL, 
                help="set of patterns", metavar="character")
    make_option(c("-sd", "--source-dataset-id"), type="character", default=NULL, 
                help="Source dataset ID", metavar="character")
    make_option(c("-sd", "--target-dataset-id"), type="character", default=NULL, 
                help="Target dataset ID", metavar="character")
    
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pattern_base_dir = '/usr/local/projects/gEAR/projectr/HuttCtxDevoLMDhs_to_ARKctxDevo3Ksc/'

input_pattern_file = paste(pattern_base_dir, opt$set-of-patterns)
target_h5_file = paste(pattern_base_dir, opt$target-dataset-id,'.h5ad')

# read target data from h5ad
library('rhdf5')
obsdata = h5read(target_h5_file, 'obs')
h5closeAll()

# If loading matrix from TAB file instead of H5AD
input_mtx = paste(pattern_base_dir, 'HuttCtxDevoLMDhs_DataMTX.tab')
my.data = as.matrix(read.delim(as.is=TRUE,file=input_mtx))

# This is, for example, a PCA of the Huttner dataset.  It could be stored within
# 'uns' and then applied to others.  For example, we could generate a menu showing that
# PCA and NMF patterns are available from this dataset to be projected onto others
CTXload = as.matrix(read.delim(as.is=TRUE, file=input_pattern_file))

# ENSGids need to be rownames, currently they are the first column
rnmsData=my.data[,1]
my.data = my.data[,-1]
rnmsLoad=CTXload[,1]
CTXload =  CTXload[,-1]
sum( rnmsData %in% rnmsLoad )

my.data=matrix(data=as.numeric(my.data),nrow=dim(my.data)[1],ncol=dim(my.data)[2])
CTXload=matrix(data=as.numeric(CTXload),nrow=dim(CTXload)[1],ncol=dim(CTXload)[2])

rownames( my.data ) = rnmsData
rownames( CTXload ) = rnmsLoad

# Do the indexes match?
sum( rownames( my.data ) %in% rownames( CTXload ) )

library(projectR)
my.proj = projectR(data=my.data, loadings = CTXload, full = FALSE)

# Get the obs data of the target of the projection, switch to h5
colmeta_path = paste(pattern_base_dir, 'HuttCtxDevoLMDhs_COLmeta.tab')
colmeta = read.delim(as.is=TRUE, file=colmeta_path)

# Here we will actually rely on the NeMO displays to define the curations
colmeta$AgeGWj=jitter(colmeta$AgeGW,amount=0.25)

# lets look at pattern 1 projection
rowN=1

png(filename='/tmp/test.png')

# product will be a matrix like an expression matrix but the row indexes are PC#
# or pattern# in the case of NMF
plot(
    x=colmeta$AgeGWj,
    y=my.proj[rowN,],
    xaxt="n",
    main="Huttner Hs Fetal LMD",
    cex=3,
    col=colmeta$Color,
    pch=19,
    xlab="Human Embryonic Age (GW)",
    ylab=paste("Projection of pattern #",rowN,sep="")
)
legend(x="topright",legend=c("CP","oSVZ","iSVZ","VZ"),text.col=c(3,4,2,1),bty="n")
axis(side=1,at=c(13:16),labels=c("GW13","GW14","GW15","GW16"))

