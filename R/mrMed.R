

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
	rho_CI_Delta <- c(lower=rho-z*rho_se,upper=rho+z*rho_se)
	#rho_CI_Delta <- c(rho-z*rho_se,rho+z*rho_se)
	#names(rho_CI_Delta) <- c(paste0(100*(gamma/2),"%"),paste0(100*(1-gamma/2),"%"))
	

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
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_X <- res$coefficients[1,1]
	se_X <- res$coefficients[1,2]/min(1,res$sigma)
	pval_X <- 2*(1-pnorm(abs(b_X/se_X)))

	cov_deltatau <- 0

	return(Mtd1(b_X,se_X,res_3$b,res_3$se,cov_deltatau[1],dat_mrMed,gamma))
}

M1_IVWa <- function(dat_mrMed,gamma=0.05){


	dat_mrMed <- form_dat(dat_mrMed)

	#estimate TE
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_X <- res$coefficients[1,1]
	se_X <- res$coefficients[1,2]/min(1,res$sigma)
	pval_X <- 2*(1-pnorm(abs(b_X/se_X)))

	#indx_mr <- (is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|ctg_Gk$Gx_plum!=1)
	Gamma_mr <- replace(dat_mrMed$beta.X,which(indx_Gxy==0),0)
	W <- diag(replace(1/dat_mrMed$se.Y^2,which(indx_mvmr==0),0))
	Gamma_mvmr <- cbind(replace(dat_mrMed$beta.X,which(indx_mvmr==0),0),replace(dat_mrMed$beta.M,which(indx_mvmr==0),0))

	sigma_uvmr <- summary(lm(dat_mrMed$beta.Y[indx_Gxy]~0+dat_mrMed$beta.X[indx_Gxy],weights=1/dat_mrMed$se.Y[indx_Gxy]^2))$sigma
	sigma_mvmr <- res$sigma

	cov_deltatau <- solve(t(Gamma_mvmr)%*%W%*%Gamma_mvmr)%*%t(Gamma_mvmr)%*%W%*%Gamma_mr/sum(t(Gamma_mr)%*%W%*%Gamma_mr)
	#cov_deltatau <- cov_deltatau*max(1,sqrt(sigma_uvmr*sigma_mvmr))^2
	#cov_deltatau <- cov_deltatau*max(1,sqrt(sigma_uvmr*sigma_mvmr))^2
	cov_deltatau <- cov_deltatau*max(1,min(sigma_uvmr,sigma_mvmr))^2

	return(Mtd1(b_X,se_X,res_3$b,res_3$se,cov_deltatau[1],dat_mrMed,gamma))
}

M1_Egger <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	dat_mrMed$beta.Y <- sign(dat_mrMed$beta.X)*dat_mrMed$beta.Y
	dat_mrMed$beta.M <- sign(dat_mrMed$beta.X)*dat_mrMed$beta.M
	dat_mrMed$beta.X <- sign(dat_mrMed$beta.X)*dat_mrMed$beta.X

	#estimate TE
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <-TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_egger_regression"))

	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_X <- res$coefficients[2,1]
	se_X <- res$coefficients[2,2]/min(1,res$sigma)
	pval_X <- 2*(1-pnorm(abs(b_X/se_X)))

	#indx_mr <- (is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|ctg_Gk$Gx_plum!=1)
	Gamma_mr <- cbind(replace(rep(1,dim(dat_mrMed)[1]),which(indx_Gxy==0),0), replace(dat_mrMed$beta.X,which(indx_Gxy==0),0))
	W <- diag(replace(1/dat_mrMed$se.Y^2,which(indx_mvmr==0),0))
	Gamma_mvmr <- cbind(replace(rep(1,dim(dat_mrMed)[1]),which(indx_mvmr==0),0),replace(dat_mrMed$beta.X,which(indx_mvmr==0),0),replace(dat_mrMed$beta.M,which(indx_mvmr==0),0))

	sigma_uvmr <- summary(lm(dat_mrMed$beta.Y[indx_Gxy]~dat_mrMed$beta.X[indx_Gxy],weights=1/dat_mrMed$se.Y[indx_Gxy]^2))$sigma
	sigma_mvmr <- res$sigma	

	cov_deltatau <- solve(t(Gamma_mvmr)%*%W%*%Gamma_mvmr)%*%t(Gamma_mvmr)%*%W%*%Gamma_mr%*%solve(t(Gamma_mr)%*%W%*%Gamma_mr)
	cov_deltatau <- cov_deltatau*max(1,min(sigma_uvmr,sigma_mvmr))^2

	return(Mtd1(b_X,se_X,res_3$b,res_3$se,cov_deltatau[2,2],dat_mrMed,gamma))
}


M1_Median <- function(dat_mrMed,gamma=0.05,Nboot=1000){

	dat_mrMed <- form_dat(dat_mrMed)

	#estimate TE
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_weighted_median"))
	u_tau <- res_3$b
	se_tau <- res_3$se
	z=qnorm(1-gamma/2)
	TE_CI <- c(lower=u_tau-z*se_tau,upper=u_tau+z*se_tau)	
	
	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	u_delta <- quantreg::rq(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2)$coefficients[[1]]
	IE <- u_tau - u_delta
	rho <- IE/u_tau

	#===CI: Bootstrap

	DE_se=c(NA);IE_se=c(NA);rho_se=c(NA)
	#IE_CI_Delta=c(NA,NA);DE_CI_Delta=c(NA,NA);rho_CI_Delta=c(NA,NA)
	IE_CI_boot=c(NA,NA);DE_CI_boot=c(NA,NA);rho_CI_boot=c(NA,NA)

	res_boot <- data.frame(IE_i=rep(0,Nboot),DE_i=rep(0,Nboot),rho_i=rep(0,Nboot))
	dat_mrMed_i <- dat_mrMed
		for(i in 1:Nboot){
			indx_X <- !is.na(dat_mrMed_i$beta.X)
			indx_M <- !is.na(dat_mrMed_i$beta.M)
			indx_Y <- !is.na(dat_mrMed_i$beta.Y)
			dat_mrMed_i$beta.X[indx_X] <- as.vector(rmvnorm(1,dat_mrMed$beta.X[indx_X],diag(dat_mrMed$se.X[indx_X]^2)))
			dat_mrMed_i$beta.M[indx_M] <- as.vector(rmvnorm(1,dat_mrMed$beta.M[indx_M],diag(dat_mrMed$se.M[indx_M]^2)))
			dat_mrMed_i$beta.Y[indx_Y] <- as.vector(rmvnorm(1,dat_mrMed$beta.Y[indx_Y],diag(dat_mrMed$se.Y[indx_Y]^2)))


			b_iv <-dat_mrMed_i$beta.Y[indx_Gxy]/dat_mrMed_i$beta.X[indx_Gxy]
			VBj <- (dat_mrMed_i$se.Y[indx_Gxy]^2)/(dat_mrMed_i$beta.X[indx_Gxy]^2) + (dat_mrMed_i$beta.Y[indx_Gxy]^2)*(dat_mrMed_i$se.X[indx_Gxy]^2)/(dat_mrMed_i$beta.X[indx_Gxy]^4)
			tau_mvmr  <- weighted_median(b_iv, 1 / VBj)
			#tau_mvmr <- mr(cbind(dat_XY_i,mr_keep),method_list=c("mr_weighted_median"))$b
			#tau_mvmr <- rq(gwasY_i[indx_Gx]~0+gwasX_i[indx_Gx],weights=1/gwasY_rmna$se.outcome[indx_Gx]^2)$coefficients[[1]]

			delta_mvmr <- rq(dat_mrMed_i$beta.Y[indx_mvmr]~0+dat_mrMed_i$beta.X[indx_mvmr]+dat_mrMed_i$beta.M[indx_mvmr],weights=1/dat_mrMed_i$se.Y[indx_mvmr]^2)$coefficients[[1]]

			res_boot[i,1] <- tau_mvmr-delta_mvmr
			res_boot[i,2] <- delta_mvmr
			res_boot[i,3] <- 1-delta_mvmr/tau_mvmr
		}
	DE_se=c(NA);IE_se=c(NA);rho_se=c(NA)
	IE_se <- sd(res_boot[,c("IE_i")])
	DE_se <- sd(res_boot[,c("DE_i")])
	rho_se <- sd(res_boot[,c("rho_i")])
	IE_CI_boot=quantile(res_boot[,c("IE_i")],c(gamma/2,1-gamma/2))
	DE_CI_boot=quantile(res_boot[,c("DE_i")],c(gamma/2,1-gamma/2))
	rho_CI_boot=quantile(res_boot[,c("rho_i")],c(gamma/2,1-gamma/2))
	names(IE_CI_boot) <- c("lower","upper")
	names(DE_CI_boot) <- c("lower","upper")
	names(rho_CI_boot) <- c("lower","upper")
	
	return(list(TE=u_tau,TE_se=se_tau,TE_CI=TE_CI,DE=u_delta,DE_se=DE_se,DE_CI=DE_CI_boot,IE=IE,IE_se=IE_se,IE_CI=IE_CI_boot,rho=rho,rho_se=rho_se,rho_CI=rho_CI_boot))
}

Mtd2 <- function(u_alpha,se_alpha,u_beta,se_beta,u_tau,se_tau,gamma=0.05){ 
	IE <- u_alpha*u_beta
	DE <- u_tau-u_alpha*u_beta
	rho <- IE/u_tau
	z=qnorm(1-gamma/2)

	TE_CI <- c(lower=u_tau-z*se_tau,upper=u_tau+z*se_tau)
	DE_se=c(NA);IE_se=c(NA);rho_se=c(NA)
	IE_CI=c(NA,NA);DE_CI=c(NA,NA);rho_CI=c(NA,NA)
		
	#===CI: Delta Method
		
	IE_var <- u_alpha^2*se_beta^2 + u_beta^2*se_alpha^2
	IE_se <- sqrt(IE_var)
	IE_CI_Delta <- c(lower=IE-z*IE_se,upper=IE+z*IE_se)

	DE_var <- se_tau^2 + IE_var
	DE_se <- sqrt(DE_var)
	DE_CI_Delta <- c(lower=DE-z*DE_se,upper=DE+z*DE_se)

	rho_var <- IE_var/u_tau^2 + IE^2*se_tau^2/u_tau^4
	rho_se <- sqrt(rho_var)	
	rho_CI_Delta <- c(lower=rho-z*rho_se,upper=rho+z*rho_se)
	#rho_CI_Delta <- c(rho-z*rho_se,rho+z*rho_se)
	#names(rho_CI_Delta) <- c(paste0(100*(gamma/2),"%"),paste0(100*(1-gamma/2),"%"))


	#====
	return(list(TE=u_tau,TE_se=se_tau,TE_CI=TE_CI,DE=DE,DE_se=DE_se,DE_CI=DE_CI_Delta,IE=IE,IE_se=IE_se,IE_CI=IE_CI_Delta,rho=rho,rho_se=rho_se,rho_CI=rho_CI_Delta))
}

M2_IVW <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)
 
	indx_Gxm <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	indx_Gmy <- !(dat_mrMed$Gm_plum==0|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y))
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))


	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]

	dat_MY <- dat_mrMed[indx_Gmy,c("SNP","M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	#estimate alpha 
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_ivw"))

	#estimate beta 
	names(dat_MY) <- gsub("M","exposure",names(dat_MY))
	names(dat_MY) <- gsub("Y","outcome",names(dat_MY))
	mr_keep <- rep(TRUE,dim(dat_MY)[1])
	res_2 <- TwoSampleMR::mr(cbind(dat_MY,mr_keep),method_list=c("mr_ivw"))

	#estimate TE
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	return(Mtd2(res_1$b,res_1$se,res_2$b,res_2$se,res_3$b,res_3$se,gamma=0.05))
} 

M2_Egger <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	indx_Gxm <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	indx_Gmy <- !(dat_mrMed$Gm_plum==0|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y))
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))


	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]

	dat_MY <- dat_mrMed[indx_Gmy,c("SNP","M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	#estimate alpha 
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_egger_regression"))

	#estimate beta 
	names(dat_MY) <- gsub("M","exposure",names(dat_MY))
	names(dat_MY) <- gsub("Y","outcome",names(dat_MY))
	mr_keep <- rep(TRUE,dim(dat_MY)[1])
	res_2 <- TwoSampleMR::mr(cbind(dat_MY,mr_keep),method_list=c("mr_egger_regression"))

	#estimate TE
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_egger_regression"))

	return(Mtd2(res_1$b,res_1$se,res_2$b,res_2$se,res_3$b,res_3$se,gamma=0.05))
} 


M2_Median <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	indx_Gxm <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	indx_Gmy <- !(dat_mrMed$Gm_plum==0|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y))
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))


	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]

	dat_MY <- dat_mrMed[indx_Gmy,c("SNP","M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	#estimate alpha 
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_weighted_median"))

	#estimate beta 
	names(dat_MY) <- gsub("M","exposure",names(dat_MY))
	names(dat_MY) <- gsub("Y","outcome",names(dat_MY))
	mr_keep <- rep(TRUE,dim(dat_MY)[1])
	res_2 <- TwoSampleMR::mr(cbind(dat_MY,mr_keep),method_list=c("mr_weighted_median"))

	#estimate TE
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_weighted_median"))

	return(Mtd2(res_1$b,res_1$se,res_2$b,res_2$se,res_3$b,res_3$se))
}
