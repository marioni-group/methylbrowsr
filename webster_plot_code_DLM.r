# # # # # # # # # #
# REGIONAL PLOTS  #
# # # # # # # # # #
#(very long) code for making nice regional plots - expects a beta matrix called 'beta' and m value matrix called mval. You will also need to change file path to Illumina manifest, used as input for 'probe.features' object

#prepare manifest object ***edit with your file path***

# Edited by Daniel McCartney for MethylBrowsR
probe.features <- anno
probe.features$chr = gsub("chr", "", probe.features$chr)



probe.features$gene_name <- unlist(lapply(probe.features$UCSC_RefGene_Name, function(x) {
  y <- strsplit(x, ";")[[1]][1]
}))

probe.features$genomic_feature <- unlist(lapply(probe.features$UCSC_RefGene_Group, function(x) {
  y <- strsplit(x, ";")[[1]][1]
}))

probe.features$genomic_region <- unlist(lapply(probe.features$genomic_feature, function(x) {
  if(is.na(x)) {
    return("IGR")
  }
  else if(x=="5'UTR" | x=="TSS200" | x=="TSS1500" | x=="1stExon") {
    return("Promoter")
  }
  else if(x=="3'UTR" | x=="Body" | x=="ExonBnd") {
    return("Body")
  }
}))

probe.features$TargetID <- probe.features$Name
rownames(probe.features) <- probe.features$TargetID


# #Only include probe features which are present in methylation data
# Not applicable - already subsetted
# probe.features.all <- probe.features
# probe.features <- subset(probe.features, rownames(probe.features) %in% rownames(mval))

#Add feature column
probe.features$feature1 <- unlist(lapply(probe.features$UCSC_RefGene_Group, function(x) {
  y <- strsplit(x, ";")[[1]][1]
  ifelse(is.na(y), "IGR", y)
}))

#EPIC
probe.features$feature1 <- ifelse(is.na(probe.features$feature), "IGR", probe.features$feature)
unique(probe.features$feature1)

probe.features$cgi <- probe.features$Relation_to_Island #changed from 'Relation_to_Island'
probe.features$cgi <- ifelse(probe.features$Relation_to_Island=="Island", "island", probe.features$cgi)
probe.features$cgi <- ifelse(probe.features$Relation_to_Island=="N_Shore", "shore", probe.features$cgi)
probe.features$cgi <- ifelse(probe.features$Relation_to_Island=="S_Shore", "shore", probe.features$cgi)
probe.features$cgi <- ifelse(probe.features$Relation_to_Island=="S_Shelf", "shelf", probe.features$cgi)
probe.features$cgi <- ifelse(probe.features$Relation_to_Island=="N_Shelf", "shelf", probe.features$cgi)
probe.features$cgi <- ifelse(probe.features$Relation_to_Island=="OpenSea", "opensea", probe.features$cgi)
head(probe.features$cgi)
unique(probe.features$cgi)

probe.features$gene <- unlist(lapply(probe.features$UCSC_RefGene_Name, function(x) {
  y <- strsplit(x, ";")[[1]][1]
}))



lwidth <- 0.5
mtitle <- 1.5
cexaxis <- 1.5
cexlab <- 1.75
cexmain <- 2

info <- samps  #replace pd_100 with your phenotype data
# group.index <- as.character(info$sex)  #Groups of interest to be plotted against each other, expected to contain 0s and 1s to indicate group
# group.index <- replace(group.index, group.index=="M", 0) #made vector binary: Female=1, Male=0
# group.index <- replace(group.index, group.index=="F", 1) #made vector binary: Female=1, Male=0
# #group.index <- replace(group.index, group.index=="metastases", 2) #made vector binary: blinded group 1=1, blinded group 2=0

# all.equal(colnames(beta), rownames(info))
# all.equal(colnames(mval), rownames(info))

# pheno.0 <- "0"  #"male"
# pheno.1 <- "1"  #"female"
#pheno.2 <- "2"  #"metastases"

col.0 <- "darkblue" # "#A1FB8E"  #col.im colour for male
# col.1 <- "#B681F7" #col.co colour for female
#col.2 <- "#8EBEFA" #col.co colour for group metastases

lty.0 <- 1
# lty.1 <- 1
#lty.2 <- 1


#load packages
# library("limma")
# library("missMethyl")
# library("DMRcate")
# library("data.table")
# library("stringr")

#Code for regional plots below adapted from code received by Dirk Paul

#Assign color to features
labelPalette.feature <- c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR","IGR", "ExonBnd")
names(labelPalette.feature) <- c("#c2a5cf", "#9970ab", "#762a83", "#c51b7d", "#b2182b", "#b35806", "#C4C8A9", "#9e2a7e")

labelPalette.island <- c("opensea", "island", "shore", "shelf")
names(labelPalette.island) <- c("white","black","darkgrey","lightgrey")

#Data information
#info <- annot #Sample annotation
 # dist <- 5000 #Range around the given site to be plotted (Added dist to function (DLM)



drawRegionalPlot <- function(probe_id, chr, start, end, title, sig.cpgs, dist=5000, # img_path, img_name, 
							dataset) {


  #probe_id: probes to plot
  #chr: only necessary when plotting regions and giving start and end
  #sig_cpgs: complete list of significant cpgs to mark them in the plot
  #start: genomic coordinate to start, can be NA, then dist will be used
  #end: genomic coordinate to end, can be NA, then dist will be used
  #title: title of the plot
  #img_path: path to save the images
  #dataset: either a beta value or m-value matrix of methylation values

  plot_mval <- TRUE;
  if(min(dataset, na.rm=T) >= 0 & max(dataset, na.rm=T) <= 1) { #If all values between 0 and 1 we have beta values
    plot_mval <- FALSE
  }

  cpg <- probe_id
  # cat("Creating plot for", img_name, "...\n")

	  if (is.na(chr)) {
		target.chr <- probe.features$chr[which(rownames(probe.features)==cpg)]
	  }
	  else {
		target.chr <- chr
	  }

  target.pos <- probe.features$pos[which(rownames(probe.features)==cpg)]

  #Start must be smaller than end
  if(!is.na(start) & !is.na(end) & start > end) {
    xtmp <- start
    start <- end
    end <- xtmp
  }

  #If start or end are not given explicitely use all probes within given distance dist
  if(is.na(start) | is.na(end)) {
    tmp <- subset(probe.features, chr==target.chr & pos>target.pos-dist & pos<target.pos+dist)
  } else {
    tmp <- subset(probe.features, chr==target.chr & pos>=start & pos<=end)
  }

  #If there are too many probes within the distance or between start and end reduce list until only <= 25 probes are left to be plotted
  # disttmp <- dist
  # while(length(rownames(tmp)) > 12) {
  #   disttmp <- disttmp-100
  #   tmp <- subset(probe.features, chr==target.chr & pos>target.pos-(disttmp) & pos<target.pos+disttmp)
  # }

  #Get all info for the plot
  gene.probes <- data.frame(
    "targetID"=rownames(tmp),
    "chr"=tmp$chr,
    "pos"=tmp$pos,
    "feature1"=tmp$feature1,
    "island"=tmp$cgi)
  gene.probes <- gene.probes[order(gene.probes$pos),]
  if(is.na(title)) {
    title <- probe.features$gene[which(rownames(probe.features)==cpg)]
  }

  #Get data
  feature.col <- names(labelPalette.feature)[match(gene.probes$feature1, labelPalette.feature)]
  island.col <- names(labelPalette.island)[match(gene.probes$island, labelPalette.island)]
  data.0 <- dataset[match(gene.probes$targetID, rownames(dataset)), ] 
  # removed group
  # data.0 <- dataset[match(gene.probes$targetID, rownames(dataset)), grep(0, group.index)]  #changed grep(pheno.0 to grep(0 #edited by Amy Oct2021
  # data.1 <- dataset[match(gene.probes$targetID, rownames(dataset)), grep(1, group.index)]  #changed grep(pheno.1 to grep(1 #edited by Amy Oct2021
  # data.2 <- dataset[match(gene.probes$targetID, rownames(dataset)), grep(2, group.index)]  #changed grep(pheno.1 to grep(2  #edited by Amy Oct2021

  if(!is.null(dim(data.0))) { #Only if there is more than 1 probe within the defined distance
    gene.info <- list(
      "chr"=gene.probes$chr,
      "pos"=gene.probes$pos,
      "feature"=gene.probes$feature,
      "island"=gene.probes$island,
      "feature.col"=feature.col,
      "island.col"=island.col,
      "data.0"=data.0) #,
      # "data.1"=data.1)

    # Plot data
#    png(paste(img_path, "/", img_name, ".png", sep=""), width=max(ifelse(length(gene.info$pos)<30, length(gene.info$pos), 30),2.75), height=7, units="in", res=300, pointsize=6)
    layout(matrix(1:2, ncol=1), respect=FALSE, widths=c(1,1), heights=c(2.75,1))
     par(mar=c(-0.25,7.2,4.5,7.2)+0.25, mgp=c(5,1,0), lwd=lwidth, cex.main=cexmain*1.5, cex.lab=cexlab*1.5, cex.axis=cexaxis*1.5) # par: bottom, left, top and right margins
    if (plot_mval==FALSE) { #Plot beta values
      # plot(1:length(gene.info$pos), rowMeans(gene.info$data.0), col=col.0, type="l", lty=lty.0, lwd=0.5, xaxt="n", yaxt="n", xlab="", ylim=c(0,1), ylab="Beta value", main="", xaxs="i", xlim=c(0.5,length(gene.info$pos)+0.5))
	  boxplot(t(gene.info$data.0), col=col.0, type="l", lty=lty.0, lwd=0.5, xaxt="n", yaxt="n", xlab="", ylim=c(0,1), ylab="Beta value", main="", xaxs="i", xlim=c(0.5,length(gene.info$pos)+0.5))
    }
    else { #Plot M values
	  ymin = min(gene.info$data.0)- 1
	  ymax = max(gene.info$data.0) + 1 # DLM: Changed m-val scale to min/max -/+1
      # plot(1:length(gene.info$pos), rowMeans(gene.info$data.0), col=col.0, type="l", lty=lty.0, lwd=0.5, xaxt="n", yaxt="n", xlab="", ylim=c(ymin,ymax), ylab="M-value", main="", xaxs="i", xlim=c(0.5,length(gene.info$pos)+0.5))
	boxplot(t(gene.info$data.0), col=col.0, type="l", lty=lty.0, lwd=0.5, xaxt="n", yaxt="n", xlab="", ylim=c(ymin,ymax), ylab="M-value", main="", xaxs="i", xlim=c(0.5,length(gene.info$pos)+0.5))
    }

    title(title, line=mtitle+1)
    # for(i in 1:length(gene.info$pos)){
      # points(rep(i, ncol(gene.info$data.0)), gene.info$data.0[i,], pch=16, col=col.0, cex=1.5)
    # }
    # points(1:length(gene.info$pos), rowMeans(gene.info$data.0), pch=16, col=col.0, cex=3.5)

    # spacr <- (-0.1)
    # lines(1:length(gene.info$pos)+spacr, rowMeans(gene.info$data.1), col=col.1, lty=lty.1, lwd=0.5)
    # for(i in 1:length(gene.info$pos)){
      # points(rep(i+spacr, ncol(gene.info$data.1)), gene.info$data.1[i,], pch=16, col=col.1, cex=1.5)
    # }
    # points(1:length(gene.info$pos)+spacr, rowMeans(gene.info$data.1), pch=16, col=col.1, cex=3.5)

    axis(side=2, lwd=0, lwd.ticks=lwidth)

    # ## Mark significant probes
    # if(plot_mval==FALSE) { #plot beta values
      # points(which(rownames(gene.info$data.0) %in% sig.cpgs==TRUE), rep(-0.025, length(which(rownames(gene.info$data.0) %in% sig.cpgs==TRUE))), pch=17, col=col.0, cex=3)
    # }
    # else { #plot m-values
      # points(which(rownames(gene.info$data.0) %in% sig.cpgs==TRUE), rep(-11.5, length(which(rownames(gene.info$data.0) %in% sig.cpgs==TRUE))), pch=17, col=col.0, cex=3)
    # }

    # ## Legend
    # if (plot_mval==FALSE) {
      # legend(ifelse(rowMeans(gene.info$data.0)[1] > 0.5, "bottomleft", "topleft"), c(pheno.0, pheno.1), col=c(col.0, col.1), pch=16, lty=c("solid","solid"), bty="n", cex=1.5)
    # } else {
      # legend(ifelse(rowMeans(gene.info$data.0)[1] > 0, "bottomleft", "topleft"), c(pheno.0, pheno.1), col=c(col.0, col.1), pch=16, lty=c("solid","solid"), bty="n", cex=1.5)
    # }

    ## Features Box
    par(mar=c(13,7.2,-0.25,7.2)+0.25, mgp=c(5,1,0))
    plot(1:length(gene.info$pos), axes=FALSE, xlab="", ylab="", ylim=c(0,0.1), xlim=c(0.5,length(gene.info$pos)+0.5), xaxs="i")
    rect(seq(1:length(gene.info$pos))-0.5, 0, seq(1:length(gene.info$pos))+0.5, 0.05, col=gene.info$feature.col, border=gene.info$feature.col)
    rect(seq(1:length(gene.info$pos))-0.5, 0.05, seq(1:length(gene.info$pos))+0.5, 0.1, col=gene.info$island.col, border=gene.info$island.col)
    axis(1, at=1:length(gene.info$pos), labels=rownames(gene.info$data.0), las=2, cex.axis=0.9) # , cex=0.9, lwd=lwidth)

	# Add legend for genomic features
	legend(x = "topright",
		   inset = c(1, 0.3),
		   legend = labelPalette.feature, 
		   fill = names(labelPalette.feature),
		   xpd = TRUE, title="Genomic Feature", bty="n"
	)

	# Add second legend
	legend(x = "topleft",
		   inset = c(1, 0.3),
		   legend = labelPalette.island, 
		   fill = names(labelPalette.island),
		   xpd = TRUE, title="CpG Island", bty="n"
	)

    
    # ## Labels (adapted to work when the same label occurrs more than once in the plot)
    # featurelabels <- rle(as.character(gene.info$feature))
    # lastpos <- 0
    # pos <- 0
    # for (i in 1:length(featurelabels$lengths)) {
      # pos <- (lastpos+(featurelabels$lengths[i]/2))+0.5
      # if(!is.na(featurelabels$values[i])) text(pos, 0.022, featurelabels$values[i], col="white", cex=2)
      # lastpos <- lastpos+featurelabels$lengths[i]
    # }

    # featurelabels <- rle(as.character(gene.info$island))
    # lastpos <- 0
    # pos <- 0
    # for (i in 1:length(featurelabels$lengths)) {
      # pos <- (lastpos+(featurelabels$lengths[i]/2))+0.5
      # if(featurelabels$values[i]=="island") text(pos, 0.072, "CGI", col="white", cex=2)
      # lastpos <- lastpos+featurelabels$lengths[i]
    # }
#    dev.off()
  } else { #Create empty plot if there is only 1 probe within the distance
#    png(paste(img_path, "/", img_name, ".png", sep=""), width=5, height=2, units="in", res=300)
    plot(1,0,col="#ffffff00", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    mtext("Too few individuals/probes. Adjust genomic range or sample parameters")
#    dev.off()
  }
}



# #try plotting ABCG1 regional plot
# #cg06500161 - ABCG1#In methyldetectR 1st (most influential) for HDL cholesterol and 2nd for BMI.#In Braun 2017 paper below, assoc with HDL-C
# rownames(probe.features) <- probe.features$Name
# probe.features["cg06500161",]   #matches up with above


# drawRegionalPlot(probe_id = "cg05593176",
                # chr = 22,
                # start = 16000000, #gene start 43619799 [cpg is at 43656587]
                # end= 18000000,  #gene end 43717354
                # title = "ABCG1 gene",
                # sig.cpgs= NA,
                # img_path="./",
                # img_name="test",
                # dataset=beta_vals)
# dev.off()
