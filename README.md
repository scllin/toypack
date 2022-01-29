# toypack
devtools::install_github("scllin/toypack")

library(toypack)

data(WHR_SMK_CAD)

M1_IVW_0(WHR_SMK_CAD)

M1_IVW(WHR_SMK_CAD)

M1_Egger(WHR_SMK_CAD)

M1_Median(WHR_SMK_CAD,Nboot=100)

M2_IVW(WHR_SMK_CAD)

M2_Egger(WHR_SMK_CAD)

M2_Median(WHR_SMK_CAD)
