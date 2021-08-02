#  examples of use
#  qqnorm.with.sim.bounds(airquality$Wind)

#  different rules for drawing red/blue/green lines
#  for this example, not much difference
#  qqnorm.with.sim.bounds(airquality$Wind,robust=F)

#  add Shapiro-Wilk test results, title, suppress legend:
#  qqnorm.with.sim.bounds(airquality$Wind,sw=T,main="Wind",legend=F)

qqnorm.with.sim.bounds<-function(the.data,sw=F,robust=T,main=NA,legend=T) {
	if (is.na(main)) {main<-"Normal Q-Q Plot"}
	n<-length(the.data)
	y.mx<-matrix(0,5000,n)
	#  robust--use median and mad; !robust--use mean and sd
	if (robust) {x.center<-median(the.data)} else {x.center<-mean(the.data)}
	if (robust) {s.center<-mad(the.data)} else {s.center<-sd(the.data)}
	for (i in 1:5000) {y.mx[i,]<-sort(rnorm(n,x.center,s.center))}
	lo.bds<-apply(y.mx,2,lo.bd<-function(x) {quantile(x,0.025)})
	hi.bds<-apply(y.mx,2,hi.bd<-function(x) {quantile(x,0.975)})
	lo.lo.bds<-apply(y.mx,2,lo.lo.bd<-function(x) {quantile(x,0.005)})
	hi.hi.bds<-apply(y.mx,2,hi.hi.bd<-function(x) {quantile(x,0.995)})
	meds<-apply(y.mx,2,median)
	ideal.x<-qnorm(ppoints(n))
	if (sw) {
plot(ideal.x,sort(the.data),ylim=c(min(c(lo.lo.bds,the.data)),max(c(the.data,hi.hi.bds))),main=main,xlab="Theoretical Quantiles",ylab="Sample Quantiles",sub=paste("SW p-value",signif(shapiro.test(the.data)$p.value,4)))
	}
	else {plot(ideal.x,sort(the.data),ylim=c(min(c(lo.lo.bds,the.data)),max(c(the.data,hi.hi.bds))),main=main,xlab="Theoretical Quantiles",ylab="Sample Quantiles")}
	abline(thud<-lm(meds~ideal.x),col="blue")
	lines(ideal.x,lo.bds,col="green",lty=2)
	lines(ideal.x,hi.bds,col="green",lty=2)
	lines(ideal.x,lo.lo.bds,col="red",lty=2)
	lines(ideal.x,hi.hi.bds,col="red",lty=2)
	if (legend) legend("topleft",legend=c("median","95% bounds","99% bounds"),lty=1,col=4:2)
}

