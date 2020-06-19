if(!require(optparse)){
    install.packages("optparse", repos='https://ftp.fau.de/cran/')
  library(optparse)
}

if(!require(mixtools)){
    install.packages("mixtools", repos='https://ftp.fau.de/cran/')
   library(mixtools)
}

 


option_list<- list(
make_option(c("-i","--input"),type="character",dest="input",default=NULL,help="input file name"),
make_option(c("-o","--output"),type="character",dest="output",default="statOut",help="output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))


x<-read.table(opt$input,header=T)

x$log <- log(x$MeanBTReadCount)
K = normalmixEM(x$log)
hist(K$x, freq=FALSE, ylim=c(0, 0.5),xlab = "median of all bootstrap dna read counts",main = "distribution of read counts")
xaxis = seq(0, 12, 0.1)
fp = 1 - pnorm(x$log, K$mu[1], K$sigma[1])
tp = pnorm(x$log, K$mu[2], K$sigma[2])
x$FP <- round(fp,digits = 2)
x$TP <- round(tp,digits=2)

write.table(x,file=opt$output,sep="\t",row.names=FALSE,quote=FALSE)

