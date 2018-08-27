

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

plot_function<-function(num=num, center=center, UTR5=FALSE, UTR3=FALSE, fulllength=FALSE,le, sRNA=FALSE, le_sRNA=le_sRNA,thres=0.99){
x<-read.csv("intarna_websrv_table.csv", sep=";", header=FALSE)

##moved code from below

pars<-read.csv("intarna_pars.txt", header=FALSE)
up<-as.numeric(as.character(pars[2,1]))
down<-as.numeric(as.character(pars[3,1]))

# selects top50
dens_num<-nrow(x)                      # old with pval thresh ||      length(which(x[,1]<=thres))
num<-nrow(x)
if(dens_num<num){
	num<-dens_num
}

s<-6
e<-7

if(sRNA==TRUE){
	s<-13
	e<-14
}

le<-1

#mi<-min(target_table[[1]][,1],na.rm = TRUE)
#ma<-max(target_table[[1]][,2],na.rm = TRUE)

mi<-1
if(sRNA==TRUE){
	x[1,4]<-le_sRNA
}
ma<-as.numeric(x[1,4])

dis<-7
collist<-c("steelblue2","orangered1","olivedrab3","darkgoldenrod","darkorchid3","dodgerblue3","chocolate3","chartreuse3")


	
if(sRNA==TRUE){
	plot(0,0, xlim=c(mi,ma*1.25), ylim=c(0,((num*dis+num-1)*2)*1.01), type="n", ylab="", main="Predicted interaction regions sRNA", xlab="",yaxt="n",bty="n",xaxt="n")
	} else {
        plot(0,0, xlim=c(mi,ma*1.25), ylim=c(0,((num*dis+num-1)*2)*1.01), type="n", ylab="", main="Predicted interaction regions mRNAs", xlab="",yaxt="n",bty="n",xaxt="n")
}

#Interaction plots
#print ("Interaction plots")
count<-0
count2<-1

for(i in 1:(num)){
	
	
	count<-count+1
	count2<-count2+1
	xx<-count2/8
	while(xx > 1) xx<-xx-1
	colnumb<-xx*8
	
	lines(c(1,x[1,4]),c(i+count*dis,i+count*dis),col="gray85", lwd=4)
	
        # remve comment here to make colors
	#lines(c(x[i,s],x[i,e]),c(i+count*dis,i+count*dis), col=collist[colnumb], lwd=4)
        # add comment here to add colors
	lines(c(x[i,s],x[i,e]),c(i+count*dis,i+count*dis), col="orangered1", lwd=4)
	text(ma*1.15, i+count*dis+0.5, x[i,3] ,cex=0.7 )
	lines(c(mi,ma*1.10),c(i+count*dis+dis/2+0.5,i+count*dis+dis/2+0.5))
}
	


	
	
#Density plot
#print ("Density plot")
density_pars<-c()
#print (dens_num)
for(i in 1:dens_num){
	density_pars<-c(density_pars,seq(x[i,s],x[i,e]))
	#print (i)
	#print (seq(x[i,s],x[i,e]))
}

hhh<-density(density_pars)
maxi<-localMaxima(hhh$y)
mmm<-max(hhh$y)
mmm2<-max(hhh$y)
mmm<-(num*dis+num-1)/mmm

hhh$y<-hhh$y*mmm+(num*dis+num-1)*1.01+(num*dis+num-1)*0.04

points(hhh, type="l", col="black", lwd=2)






#Ablines for local maxima in density plot
#print ("ABlines for local maxima")
for(ii in 1:length(maxi)){
	lines(c(hhh$x[maxi[ii]],hhh$x[maxi[ii]]),c((num*dis+num-1)*1.01+(num*dis+num-1)*0.04,hhh$y[maxi[ii]]),lwd=2)
	}
#parting line between density plot and interaction plots
abline(h=(num*dis+num-1)*1.01+(num*dis+num-1)*0.04)

#Axis for density plot
#print ("Axis for density plot")
if(sRNA==FALSE){
	ylen<-seq(0,mmm2+0.02, by=0.001)
	ylen1<-ylen
	ylen<-ylen*mmm+(num*dis+num-1)*1.01+(num*dis+num-1)*0.04
	axis(2,at=ylen,labels=ylen1)
	}

else{
	ylen<-seq(0,mmm2+0.02, by=0.01)
	ylen1<-ylen
	ylen<-ylen*mmm+(num*dis+num-1)*1.01+(num*dis+num-1)*0.04
	axis(2,at=ylen,labels=ylen1)
	}
	
#Annotations for interaction plot
#print ("Annotations for interaction plot")
if(UTR5==TRUE){
	#lines(c(center*1.2,center*1.2),c(-100,num+num*dis+5), col="gray50")
	lines(c(center,center),c(-100,num+num*dis+5), col="gray50")

	#se3<-c(center*1.2-up,center*1.2-25,center*1.2,center*1.2+25,center*1.2+down)
	se3<-c(center-up,center-(up/2),center,center+(down/2),center+down)

	axis(1,at=se3,labels=c(up,(up/2),"start codon",(down/2),down))
}


if(UTR3==TRUE){
	#lines(c(center*1.2,center*1.2),c(-100,num+num*dis+5), col="gray50")
	lines(c(center,center),c(-100,num+num*dis+5), col="gray50")

	#se3<-c(center*1.2-50,center*1.2-25,center*1.2,center*1.2+25,center*1.2+50)
	se3<-c(center-up,center-(up/2),center,center+(down/2),center+down)

	axis(1,at=se3,labels=c(up,up/2,"stop codon",down/2,down))
}

if(fulllength==TRUE){


	se<-seq(mi,ma,by=100)
	selabels<-seq(0,100*(length(se)-1),by=100)
	axis(1,at=se,labels=selabels)
	title(xlab="nt")
	}

if(sRNA==TRUE){
	xval<-ma
	leng<-seq(0,xval,by=50)
	leng[1]<-1
		if((xval-max(leng))>20){
		leng<-c(leng, xval)}
		else{
		leng[length(leng)]<-xval
		}
	axis(1,at=leng,labels=leng)
	title(xlab="nt")
	}
}

#plot_function(num=25, center=70, sRNA=TRUE, le_sRNA=201)

IntaRNA_plot<-function(thres=0.99, num=nrow(x),name="ncRNA" ){
require(seqinr)
sRNA<-read.fasta("binding_RNAs.fa")
le_sRNA<-length(sRNA[[1]])

na<-"sRNA_regions.pdf"
na_ps<-"sRNA_regions.ps"

pdf(file=na, paper="a4r", width=0, height=0)
plot_function(num, sRNA=TRUE,le_sRNA=le_sRNA,thres=thres) ## error here: 'from' cannot be NA, NaN or infinite
dev.off()

postscript(file=na_ps)
plot_function(num, sRNA=TRUE,le_sRNA=le_sRNA,thres=thres)  
dev.off()

pars<-read.csv("intarna_pars.txt", header=FALSE)
up<-as.numeric(as.character(pars[2,1]))
down<-as.numeric(as.character(pars[3,1]))
type<-as.character(pars[1,1])


UTR5<-FALSE
UTR3<-FALSE
fulllength<-FALSE
if(type=="CDS"){
	fulllength<-TRUE
}

if(type=="5utr"){
	UTR5<-TRUE
	center<-up
}

if(type=="3utr"){
	center<-up
	UTR3<-TRUE
}


na<-"mRNA_regions.pdf"
na_ps<-"mRNA_regions.ps"

pdf(file=na, paper="a4r", width=0, height=0)
plot_function(num, center,UTR5=UTR5,UTR3=UTR3,fulllength=fulllength,le_sRNA=le_sRNA,thres=thres)
dev.off()

postscript(file=na_ps)
plot_function(num, center,UTR5=UTR5,UTR3=UTR3,fulllength=fulllength,le_sRNA=le_sRNA,thres=thres)
dev.off()

}	
	IntaRNA_plot()
	
