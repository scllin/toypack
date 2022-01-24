# toypack
devtools::install_github("scllin/toypack")

library(toypack)

data(WHR_SML_CAD)

M1_IVW_0(WHR_SML_CAD)

M1_IVW(WHR_SML_CAD)

M1_Egger(WHR_SML_CAD)

M1_Median(WHR_SML_CAD,Nboot=100)

M2_IVW(WHR_SML_CAD)

M2_Egger(WHR_SML_CAD)

M2_Median(WHR_SML_CAD)
