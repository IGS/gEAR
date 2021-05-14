#!/usr/bin/env Rscript

#############################
# 
#############################

args <- commandArgs(TRUE)
#gene0 = "E2F1"
gene0 = args[1]

#lb0="BrainSpanDFC"
lb0 = args[2]

#data_path = "/usr/local/projects/gEAR/carlo/"
data_path = args[3]

#image_path=paste(data_path, "GenePlot_", lb0, "_", gene0, ".png", sep="")
image_path = args[4]

load(file=paste(data_path,lb0,"_data.RData",sep=""))
load(file=paste(data_path,lb0,"_dataCOLanno.RData",sep=""))
load(file=paste(data_path,lb0,"_dataROWanno.RData",sep=""))

# colID;xlb0;ylo0;spn0;abl0;ablb0;rowID;X0;clr0;pch0;YexprsLB0;ROWplotID

dataROWanno.ROWindx=match(gene0,dataROWanno[,ROWplotID])
data.ROWindx=match(dataROWanno[dataROWanno.ROWindx,rowID],rownames(data))

png( width=780,
     height=480,
     file=image_path)

main_label = paste(lb0,": ",gene0,sep="")
y_label = paste(gene0," Expression: ",YexprsLB0,sep="")

plot( x=dataCOLanno[,X0],
      y=data[data.ROWindx,],
      col=dataCOLanno[,clr0],
      pch=dataCOLanno[,pch0],
      cex=2,
      main=main_label,
      xlab=xlb0,
      ylab=y_label )

if (ylo0) {
   lines( x=dataCOLanno[,X0],
          y=loess(data[data.ROWindx,]~dataCOLanno[,X0],spn=spn0)$fit
        )
}

abline(v=abl0,lty=2);text(x=abl0,y=par()$yaxp[1],labels=ablb0)

dev.off()
