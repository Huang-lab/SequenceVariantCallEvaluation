##### examine_batch_ethnicity_effects.R #####
# Kuan-lin Huang @ MSSM

# implementation of LMG to find independent contribution of regressors
# library("relaimpo") # couldn't get this to work, seems to require all quantitative variables
library("hier.part")

system("mkdir out")

# take input arguments
args = commandArgs(trailingOnly=TRUE)
# read in files
if (length(args)<1) {
  stop("At least two argument must be supplied [input file] [output tag name].n", call.=FALSE)
}
input_file = args[1]
out_name = args[2]

# read in clinical data file
clin_f = "PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)

# read in and process meta data file
meta_f = "tcga_meta_data.normals_only.Cases_N10389.txt"
meta = read.table(header=F, sep="\t", file=meta_f, fill=T)
colnames(meta) = c("bcr_patient_barcode","ID","cancer","X","assay","center","date","platform","Y","reagent","capture")
meta$capture_brief = gsub("\\|.*","",meta$capture)
meta$year = gsub(".*/.*/([0-9]+)","\\1",meta$date)
meta = meta[order(meta$year,decreasing = T),]# keep only the later year
meta_uniq = meta[!duplicated(meta$bcr_patient_barcode),]
meta_uniq$analyte = gsub("TCGA-..-....-...-..(.)-.*","\\1",meta_uniq$ID)
meta_uniq$analyte[meta_uniq$analyte=="X"] = "D"

# read in and process istat variant count file
istat_all = read.table(header=F, sep="\t", file=gzfile(input_file), fill=T)
colnames(istat_all) = c("ID","NALT","NMIN","NHET","NVAR","RATE","SING","TITV","PASS","PASS_S","QUAL","DP","geneName","fileName")
istat_all$geneName = gsub("istat/","",istat_all$geneName)
istat_all$bcr_patient_barcode = substr(istat_all$ID,1,12)

##### individual level stats #####

for (gene in unique(istat_all$geneName)){
  tryCatch({
  # gene = unique(istat_all$geneName)[1] # trouble-shoot
  istat = istat_all[istat_all$geneName==gene,]
  istat_clin = merge(istat,clin[,c("bcr_patient_barcode","consensus_call")],by="bcr_patient_barcode")
  istat_clin_meta = merge(istat_clin, meta_uniq[,-which(colnames(meta_uniq) %in% c("bcr_sample_barcode","X","Y","assay"))], by=c("ID"))
  
  # note: can consider using capture_brief instead to avoid over-fitting/attribution
  # LMG test
  istat_clin_eur = istat_clin_meta[istat_clin_meta$consensus_call=="eur",]
  # pdf("variation_BRCA1.pdf",width=10)
  variation_explained = hier.part(istat_clin_meta$NVAR, istat_clin_meta[,c("capture","platform","reagent","year","center","analyte","consensus_call")], family = "gaussian", gof = "Rsqu")$I.perc
  # dev.off()
  variation_explained_eur = hier.part(istat_clin_eur$NVAR, istat_clin_eur[,c("capture","platform","reagent","year","center","analyte")], family = "gaussian", gof = "Rsqu")$I.perc
  
  # # regression test
  # fit = lm(data = istat_clin, NVAR ~ consensus_call + type + Center + Analyte) # the order of importance found using LMG
  # summary(fit)
  # results = data.frame(anova(fit))
  # 
  # fit = lm(data = istat_clin, NVAR ~ consensus_call + type + Center + Analyte) # the order of importance found using LMG
  # summary(fit)
  # results = data.frame(anova(fit))
  # 
  # fit = lm(data = istat_clin_eur, NVAR ~ type + Center + Analyte) # the order of importance found using LMG
  # summary(fit)
  # results = data.frame(anova(fit))
  
  # concatenate gene name into output 
  variation_explained$gene = gene
  variation_explained_eur$gene = gene
  
  write.table(variation_explained, quote=F,col.names = F,file = paste("out/all_gene_istat_full_exon_",out_name,".tsv",sep=""),
              append = TRUE, sep = "\t")
  write.table(variation_explained_eur, quote=F,col.names = F,file = paste("out/all_gene_istat_eur_full_exon_",out_name,".tsv",sep=""),
              append = TRUE, sep = "\t")
  }, error=function(e){cat("ERROR for gene:",gene, "\n")})
}

