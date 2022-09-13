setwd("/Cluster_Filespace/Marioni_Group/Daniel/MethylBrowsR")

# dir.create("datafiles")
dat = readRDS("wave3_mvals.rds")

# Read Anno
anno = readRDS("EPIC_AnnotationObject_df.rds")


# Subset to DNAm sites
anno = anno[rownames(dat),]
# Order by chr/pos
anno = anno[order(anno$chr, anno$pos), ]
dat = dat[rownames(anno),]

ind = round(seq(1, nrow(dat), length=1000))
for(i in ind[-length(ind)]){
	file = which(ind == i)
	start = ind[file]
	if(file==length(ind[-length(ind)])){
		end = ind[length(ind)]
		}
		else {
			end = ind[file+1]-1
		}
	out = dat[start:end,]
	saveRDS(out, file=paste0("datafiles/dat_", file, ".rds"))
}

anno$file_ind = NA
for(i in ind){
	d = length(ind)
	d1 = which(ind==i)
	if(i < ind[d]){
		anno$file_ind[i:(ind[d1+1]-1)] = d1
	} else {
		anno$file_ind[i:nrow(anno)] = length(ind)-1
	}
}

saveRDS(anno, file="EPIC_AnnotationObject_df_fileind.rds")