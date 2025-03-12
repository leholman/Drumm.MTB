####### Data cleaning for Drumm #####

## some functions 

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

CountTable <- function(in.taxonomy,in.data,output="Count",some.unassigned=T){
  if(length(in.taxonomy)!=length(in.data[,1])){stop("Dataframe and corresponding taxonomy are not the same length")}
  in.taxonomy[is.na(in.taxonomy)] <- ""
  out.dat <- as.data.frame(matrix(ncol=length(in.data[1,]),nrow=length(unique(in.taxonomy))))
  rownames(out.dat) <- sort(unique(in.taxonomy))
  colnames(out.dat) <- colnames(in.data)    
  out.dat.abundance <- out.dat
  for (sample in 1:length(in.data[1,])){
    out.dat[,sample] <- table(in.taxonomy[in.data[,sample]>0])[match(sort(unique(in.taxonomy)),names(table(in.taxonomy[in.data[,sample]>0])))]
    out.dat.abundance[,sample] <- aggregate(in.data[,sample], by=list(Category=in.taxonomy), FUN=sum)[,2]
  }
  out.dat[is.na(out.dat)] <- 0
  if(some.unassigned==T){rownames(out.dat)[1] <- "Unassigned"}
  if(output=="Count"){return(out.dat)}else if(
    output=="Abundance"){return(out.dat.abundance)}
}

minAbundance <- function(inputtable = NA, minAbun = 0.01) {
  others <- rep(0, ncol(inputtable))
  
  for (col in 1:ncol(inputtable)) {
    threshold = sum(inputtable[, col]) * minAbun
    below_threshold_indices = inputtable[, col] < threshold
    others[col] = sum(inputtable[below_threshold_indices, col])
    inputtable[below_threshold_indices, col] = 0
  }
  
  inputtable <- rbind(inputtable, others)
  rownames(inputtable)[nrow(inputtable)] = "Others"
  
  # Remove rows that sum to zero
  inputtable <- inputtable[rowSums(inputtable) != 0, ]
  
  return(inputtable)
}
getPalette = colorRampPalette(brewer.pal(9, "Set3"))


library(RColorBrewer)

### read in metadata 

metadata <- readxl::read_excel("LV til kilde_nix pille.xlsx")
metadata$`Nyt navn uden formel`

metadata$sample <- sapply(strsplit(metadata$`Nyt navn uden formel`, "\\."), `[`, 2)
metadata$exp <- sapply(strsplit(metadata$`Nyt navn uden formel`, "\\."), `[`, 1)


## First V9
eukv9.dat <- read.csv("rawdata/eukv9.raw.names.csv.gz")
eukv9.asv <- seqinr::read.fasta("rawdata/OTUS/EUKv9.DADA2.ASVs.fasta",as.string = TRUE)
eukv9.pr2 <- read.csv("taxonomy/eukv9.tax.PR2.csv")

unique(metadata$exp)

metadata$LV[metadata$sample=="C"]

controls <- eukv9.dat[,match(metadata$LV[metadata$sample=="C"],colnames(eukv9.dat))]
colnames(controls)

waterAVR <- rowMeans(eukv9.dat[,match(metadata$LV[metadata$sample=="W"],colnames(eukv9.dat))])
sedAVR <- rowMeans(eukv9.dat[,match(metadata$LV[metadata$sample=="S"],colnames(eukv9.dat))])


## division
eukv9.c <-   as.matrix(minAbundance(CountTable(as.character(eukv9.pr2$Division),cbind(controls,sedAVR,waterAVR),output = "Abundance"),minAbun=0.01))
row.names(eukv9.c)[1] <- "Unknown"

pdf("figures/controls.div.pdf",width=11,height=6)
par(mar=c(8.5, 6.1, 1.1, 6.1),xpd=TRUE)
barplot(eukv9.c[,dim(eukv9.c)[2]:1],las=2,cex.names=0.6,col=getPalette(dim(eukv9.c)[1]),ylab="",border = NA,
        names=rev(c(metadata$`Nyt navn uden formel`[match(colnames(controls),metadata$LV)],"sediment_mean","water_mean")))
mtext("Read Abundance",side = 2,line=5)
legend(48.5,1200000,rev(rownames(eukv9.c)),col=getPalette(dim(eukv9.c)[1]),cex=0.7,pch=15,pt.cex = 2,bty = "n", xpd = TRUE)
dev.off()


## class
eukv9.c <-   as.matrix(minAbundance(CountTable(as.character(eukv9.pr2$Class),cbind(controls,sedAVR,waterAVR),output = "Abundance"),minAbun=0.01))
row.names(eukv9.c)[1] <- "Unknown"

pdf("figures/controls.cls.pdf",width=11,height=6)
par(mar=c(8.5, 6.1, 1.1, 6.1),xpd=TRUE)
barplot(eukv9.c[,dim(eukv9.c)[2]:1],las=2,cex.names=0.6,col=getPalette(dim(eukv9.c)[1]),ylab="",border = NA,
        names=rev(c(metadata$`Nyt navn uden formel`[match(colnames(controls),metadata$LV)],"sediment_mean","water_mean")))
mtext("Read Abundance",side = 2,line=5)
legend(48.5,1400000,rev(rownames(eukv9.c)),col=getPalette(dim(eukv9.c)[1]),cex=0.4,pch=15,pt.cex = 2,bty = "n", xpd = TRUE)
dev.off()

## division prop
eukv9.c <-   as.matrix(minAbundance(CountTable(as.character(eukv9.pr2$Division),prop.table(as.matrix(cbind(controls,sedAVR,waterAVR)),margin=2),output = "Abundance"),minAbun=0.01))
row.names(eukv9.c)[1] <- "Unknown"

pdf("figures/controls.div.prop.pdf",width=11,height=6)
par(mar=c(8.5, 6.1, 1.1, 6.1),xpd=TRUE)
barplot(eukv9.c[,dim(eukv9.c)[2]:1],las=2,cex.names=0.6,col=getPalette(dim(eukv9.c)[1]),ylab="",border = NA,
        names=rev(c(metadata$`Nyt navn uden formel`[match(colnames(controls),metadata$LV)],"sediment_mean","water_mean")))
mtext("Read Abundance",side = 2,line=5)
legend(48.5,1,rev(rownames(eukv9.c)),col=getPalette(dim(eukv9.c)[1]),cex=0.7,pch=15,pt.cex = 2,bty = "n", xpd = TRUE)
dev.off()

## class prop
eukv9.c <-   as.matrix(minAbundance(CountTable(as.character(eukv9.pr2$Class),prop.table(as.matrix(cbind(controls,sedAVR,waterAVR)),margin=2),output = "Abundance"),minAbun=0.01))
row.names(eukv9.c)[1] <- "Unknown"

pdf("figures/controls.class.prop.pdf",width=11,height=6)
par(mar=c(8.5, 6.1, 1.1, 6.1),xpd=TRUE)
barplot(eukv9.c[,dim(eukv9.c)[2]:1],las=2,cex.names=0.6,col=getPalette(dim(eukv9.c)[1]),ylab="",border = NA,
        names=rev(c(metadata$`Nyt navn uden formel`[match(colnames(controls),metadata$LV)],"sediment_mean","water_mean")))
mtext("Read Abundance",side = 2,line=5)
legend(48.5,1,rev(rownames(eukv9.c)),col=getPalette(dim(eukv9.c)[1]),cex=0.4,pch=15,pt.cex = 2,bty = "n", xpd = TRUE)
dev.off()




c(metadata$`Nyt navn uden formel`[match(colnames(controls),metadata$LV)],"sediment_mean","water_mean")
  
match(colnames(controls),metadata$LV)
colnames(controls)





