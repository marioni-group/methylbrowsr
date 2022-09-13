# Define M2Beta function
M2Beta <-
function (M) 
{
    return((2^M)/(2^M + 1))
}


# Load required data to be used throughout lifetime of app
samps <- readRDS("samps.rds")
# samps <- readRDS("Z:/Daniel/MethylBrowsR/GS20k_Targets.rds")
samps$Sex <- factor(samps$sex, labels = c("Female", "Male"))

# mvals = readRDS("Z:/Daniel/MethylBrowsR/apoe_elovl2_dat.rds")
# beta_vals = M2Beta(mvals)
anno = readRDS("EPIC_AnnotationObject_df_fileind.rds")



