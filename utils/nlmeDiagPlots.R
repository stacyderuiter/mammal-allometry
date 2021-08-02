library(nlme)

#  Routine to plot residuals against fitted, qqnorm plot of 
#  residuals, plus residuals against any column(s) in the 
#  dataset

#  Input is:
#  model:  the model you want diagnostics for
#  x:  the dataset
#  col.nos.to.plot:  the columns (as column numbers) you want
#                    to plot against the residuals

add.loess.line<-function(f,r) {
	rf.model<-loess(r~f,span=0.75)
	f.seq<-seq(min(f),max(f),(max(f)-min(f))/100)
	r.seq<-predict(rf.model,newdata=data.frame(f=f.seq))
	lines(f.seq,r.seq,col="red")
}

nlme.diag.plots<-function(model,x,col.nos.to.plot) {	
	n.plot.cols<-1+trunc((length(col.nos.to.plot)+1)/2)
	par(mfcol=c(2,n.plot.cols))
	r<-resid(model,type="normalized")  #  just to save typing
	plot(r~fitted(model),xlab="Fitted values",ylab="Normalized residuals")
	abline(h=0,lty=3)
	add.loess.line(fitted(model),r)
	qqnorm(r)
	#  if there are NA's in the dataset, r and x[,j] will
	#  have different lengths, resulting in an error.
	if (nrow(x)>length(r))  {x<-na.omit(x)}
	for (j in col.nos.to.plot) {plot(r~x[,j],xlab=colnames(x)[j])}
}

#  Example:  I've commented out the commands here, run each without #
#  str(Orthodont)
#  explain distance using age and Sex, with Subject as random factor
#  example.model<-lme(distance~age*Sex,random=~1|Subject,data=Orthodont)
#  we want to plot residuals against age, Subject, and Sex, that's
#  columns 2, 3, and 4 of Orthodont
#  nlme.diag.plots(example.model,Orthodont,2:4)

#  BUT DON'T ANSWER YET!!  BONUS FUNCTION!!
#  boxplot.color.by.numbers produces a boxplot of y~x where each box
#  is colored by factor f (MUST be a factor, else an error)

#  syntax:  boxplot.color.by.numbers(y~x|f,data=your.data.frame,legend.where=NA)
#  by default, the legend appears top center; if you want it somewhere else,
#  specify legend.where as needed:  see example below

boxplot.color.by.numbers<-function(formula,data,legend.where=NA) {
	f.string<-as.character(formula)
	y.name<-f.string[2]
	y.col<-(1:ncol(data))[colnames(data)==y.name]
	splat<-strsplit(f.string[3],"|")
	x.name<-""; i<-1
	while (splat[[1]][i]!=" ") {x.name<-paste(x.name,splat[[1]][i],sep="");i<-i+1}
	x.col<-(1:ncol(data))[colnames(data)==x.name]
	f.name<-""; i<-i+3
	while (i<=length(splat[[1]])) {f.name<-paste(f.name,splat[[1]][i],sep="");i<-i+1}
	f.col<-(1:ncol(data))[colnames(data)==f.name]
	data.sorted<-data[order(data[,x.col]),]
	first.rows<-1
	for (i in 2:nrow(data.sorted)) {
		if (data.sorted[i,x.col]!=data.sorted[i-1,x.col]) {first.rows<-c(first.rows,i)}
	}
	boxplot(data[,y.col]~data[,x.col],col=sapply(data[first.rows,f.col],color.me<-function(x) {(1:length(unique(data[,f.col])))[levels(data[,f.col])==x]}))
	if (is.na(legend.where)) {legend.where="top"}
	legend(legend.where,legend=levels(data[,f.col]),pch=15,col=1:length(levels(data[,f.col])))
}

#  Example:  I've commented out the commands here, run each without #
#  boxplot of distance by Subject, colored by Sex, from Orthodont
#  boxplot.color.by.numbers(distance~Subject|Sex,data=Orthodont)
#  legend in the wrong place, it would look better in the upper right
#  hand corner:
#  boxplot.color.by.numbers(distance~Subject|Sex,data=Orthodont,"topright")