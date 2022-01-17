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
