
#' Overdispersed ProxECAT
#'
#' Runs the overdispersed Proxy External Controls Association Test (ProxECAT) method given the number of rare alleles in each gene region and an overdispersion parameter
#'
#' @param x1 vector of the number of functional rare alleles in cases for each gene region
#' @param x2 vector of the number of synonymous rare alleles (or other proxy) in cases for each gene region
#' @param x3 vector of the number of functional rare alleles in controls for each gene region
#' @param x4 vector of the number of synonymous rare alleles (or other proxy) in controls for each gene region
#' @param size overdispersion parameter. The smaller size, the more over dispersed. The larger size the less overdispersed.  Size of >1000 is approximately Poisson.
#'
#' @return input data, likelihood ratio test statistic and p-value for overdispersed ProxECAT. LRT is asymptotically chi-squared with 1 df.
#'
#' @examples
#' ##runs overdispersed proxy-ECAT for 5 gene regions
#' x1<-c(12,10,5,103,89)
#' x2<-c(6,5,3,10,3)
#' x3<-c(20,20,8,5,8)
#' x4<-c(9,11,8,7,3)
#' temp<-proxecat_over(x1=x1, x2=x2, x3=x3, x4=x4, size=100)
#' temp
#'
#'@importFrom stats dnbinom dpois
#'
#' @export
proxecat_over<-function(x1, x2,x3,x4, size=1000){
  if(min(c(x1, x2, x3, x4))<0) {warning(paste("row(s)", paste(which(apply(cbind(x1, x2, x3, x4), 1, min)<0), collapse=","), "contain at least one negative value"))}
  
  den<-x1+x2+x3+x4
lambda1<-(x1**2+x1*x2+x1*x3+x2*x3)/den
lambda2<-(x2**2+x1*x2+x1*x4+x4*x2)/den
lambda3<-(x3**2+x1*x3+x3*x4+x1*x4)/den
lambda4<-(x4**2+x3*x4+x2*x3+x4*x2)/den

L_const<-dnbinom(round(x1),size=size,mu=lambda1)*dnbinom(round(x2),size=size,mu=lambda2)*dnbinom(round(x3),size=size,mu=lambda3)*dnbinom(round(x4),size=size,mu=lambda4)
L_unconst<-dnbinom(round(x1),size=size,mu=x1)*dnbinom(round(x2),size=size,mu=x2)*dnbinom(round(x3),size=size,mu=x3)*dnbinom(round(x4),size=size,mu=x4) 
test.statistic <- -2*log(L_const/L_unconst)
p.value <- pchisq(test.statistic, 1, lower.tail=F)
dat<-data.frame(x1=x1, x2=x2, x3=x3, x4=x4, size=size)

out<-list(dat = dat, test.statistic = test.statistic, p.value = p.value)
return(out)
}#end of function