# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function(x) {
  x^2
}

form_dat <- function(dat_mrMed){
	if(!"id.X"%in%colnames(dat_mrMed)){id.X <- rep("X",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,id.X)}
	if(!"X"%in%colnames(dat_mrMed)){X <- rep("X",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,X)}
	if(!"id.M"%in%colnames(dat_mrMed)){id.M <- rep("M",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,id.M)}
	if(!"M"%in%colnames(dat_mrMed)){M <- rep("M",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,M)}
	if(!"id.Y"%in%colnames(dat_mrMed)){id.Y <- rep("Y",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,id.Y)}
	if(!"Y"%in%colnames(dat_mrMed)){Y <- rep("Y",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,Y)}
	return(dat_mrMed)
}

#CI_mtd: 0:no  1:Delta Method 2:Bootstrap 3:ALL
Mtd1 <- function(u_delta,se_delta,u_tau,se_tau,cov_deltatau,dat_mrMed,gamma=0.05){ 

	#point estimates
	IE <- u_tau-u_delta
	DE <- u_delta
	rho <- IE/u_tau
	z=qnorm(1-gamma/2)

	TE_CI <- c(lower=u_tau-z*se_tau,upper=u_tau+z*se_tau)	
	DE_se=c(NA);IE_se=c(NA);rho_se=c(NA)
	IE_CI_Delta=c(NA,NA);DE_CI_Delta=c(NA,NA);rho_CI_Delta=c(NA,NA)
	#IE_CI_boot=c(NA,NA);DE_CI_boot=c(NA,NA);rho_CI_boot=c(NA,NA)

	#===CI: Delta method		
	#IE_var <- max(se_tau^2 + se_delta^2 - 2*cov_deltatau,(se_tau-se_delta)^2)
	IE_var <- se_tau^2 + se_delta^2 - 2*cov_deltatau
	IE_se <- sqrt(IE_var)
	IE_CI_Delta <- c(lower=IE-z*IE_se,upper=IE+z*IE_se)

	DE_var <- se_delta^2
	DE_se <- se_delta	
	DE_CI_Delta <- c(lower=DE-z*DE_se,upper=DE+z*DE_se)

	rho_var <- max(DE_var/u_tau^2 + DE^2*se_tau^2/u_tau^4 - 2*DE*cov_deltatau/u_tau^3,0)
	rho_se <- sqrt(rho_var)
	#rho_CI_Delta <- c(lower=rho-z*rho_se,upper=rho+z*rho_se)
	rho_CI_Delta <- c(rho-z*rho_se,rho+z*rho_se)
	names(rho_CI_Delta) <- c(paste0(100*(gamma/2),"%"),paste0(100*(1-gamma/2),"%"))
	

	return(list(TE=u_tau,TE_se=se_tau,TE_CI=TE_CI,DE=DE,DE_se=DE_se,DE_CI=DE_CI_Delta,IE=IE,IE_se=IE_se,IE_CI=IE_CI_Delta,rho=rho,rho_se=rho_se,rho_CI=rho_CI_Delta))
}


M1_IVW <- function(dat_mrMed,gamma=0.05){


	dat_mrMed <- form_dat(dat_mrMed)

	#estimate TE using Gx rather than Gx_plum
	indx_Gxy <- !(dat_mrMed$Gx==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <-mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_X <- res$coefficients[1,1]
	se_X <- res$coefficients[1,2]/min(1,res$sigma)
	pval_X <- 2*(1-pnorm(abs(b_X/se_X)))

	cov_deltatau <- 0

	return(Mtd1(b_X,se_X,res_3$b,res_3$se,cov_deltatau[1],dat_mrMed,gamma))
}
