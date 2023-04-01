library(QoRTs)

qortsdir <- "/athena/angsd/scratch/yus4008/project/dataset/qc_analysis/qorts/"
resdir <- paste0(qortsdir,"results/")
dcdtxt <- paste0(qortsdir,"decoder.txt")

decoder.data <- read.table(dcdtxt,header=T,stringsAsFactors=F)
res <- read.qc.results.data(resdir,decoder=decoder.data,calc.DESeq2=T,calc.edgeR=T)
makeMultiPlot.colorByGroup(res,outfile.dir=qortsdir,plot.device.name="pdf")