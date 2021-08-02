#  MuMIn utilities:

#  MuMIn.shrink.fit returns the predicted y-values based on the model averaged
#  coefficients with shrinkage.  It will need a copy of the model average summary
#  and the model used for averaging (commonly, the full model).

#  Example:

#  aq.no.na<-na.omit(airquality)  #  remove NA rows from the airquality dataset
#  big.model<-lm(Ozone~.,data=aq.no.na)
#  big.model.avg<-summary(model.avg(dredge(big.model)))
#  fitted.Ozone<-MuMIn.shrink.fit(big.model.avg,big.model)

MuMIn.shrink.fit<-function(avg.model,the.model) {
	the.coefs<-avg.model$coef.shrinkage
  X<-model.matrix(the.model)
  X<-X[,names(the.coefs)]
	return(X%*%the.coefs)
}

#  MuMIn.shrink.resid returns the residuals of the model using average coefficients
#  with shrinkage (in other words, residuals based on the results of MuMIn.shrink.fit).
#  It will need a copy of the model average summary and the original model and, 
#  optionally, the fitted values from MuMIn.shrink.fit

MuMIn.shrink.resid<-function(avg.model,the.model,shrink.fit=c()) {
	if (length(shrink.fit)==0) {shrink.fit<-MuMIn.shrink.fit(avg.model,the.model)}
	y<-the.model$model[,1]
	return(y-shrink.fit)
}

MuMIn.resid.df<-function(avg.model) {
      return(attr(avg.model,"nobs")-sum(avg.model$summary$df*avg.model$summary$Weight))
}
