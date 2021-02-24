# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#TCGA use case 1

# biomaRt and rms packages each take several minutes to install
# to improve the load time of RStudio, they are optionally included here
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
#BiocManager::install("rms")

## R setup to connect to the Athena JDBC and example queries
options(java.parameters = "-Xmx16g")
stringsAsFactors = FALSE

#set up libraries needed
library(rJava)
library(RJDBC)
library(AWR.Athena)
library(biomaRt)
library(DESeq2)
library(RColorBrewer)

##connect to athena through JDBC
##Use the region you are working in, eg 'us-east-1', and the output S3 bucket here

# Connect to Athena & Output bucket. Note that the bucket must be in the same region.
#con <- dbConnect(AWR.Athena::Athena(), region='<REGION>', S3OutputLocation='s3://<BUCKETNAME>', Schema='default')
con <- dbConnect(AWR.Athena::Athena(), region='us-west-2', S3OutputLocation='s3://newtcgabucket', Schema='default')

# Use Case: Cluster Breast Cancer samples by gene expression and identify patterns in survival
# GDC- Users are encouraged to normalize raw read count values if a subset of genes is investigated.
brca_expression<-dbGetQuery(con,"SELECT * from tcgatables.tcga_expression")

#brca_pt<-dbGetQuery(con,"SELECT bcr_patient_barcode, gender from tcgatables.tcga_clinical_patient_brca")
brca_pt<-dbGetQuery(con,"SELECT * from tcgatables.tcga_clinical_patient_brca")

# Tables are stored in a gene by sample matrix
dim(brca_expression)

# Prep the data - first column contains ensemble geneIDs and last is the projectID
row.names(brca_expression) <- brca_expression[,1]
brca_expression<- brca_expression[,-c(1, ncol(brca_expression))]
colnames(brca_expression) <- toupper(colnames(brca_expression))
brca_expression <- data.matrix(brca_expression)

row.names(brca_pt) <- brca_pt$bcr_patient_barcode
brca_pt <- brca_pt[-c(1,2),]
# first 2 rows of clinical data are header/metadata

## Normalized Data - transpose for sample (row) x gene (cols)
brca_log2 <- apply(brca_expression, 1, function(row){log2(row+1)})
brca_vsd <- t(vst(brca_expression, blind=FALSE))
brca_norm <- brca_vsd

# Cluster by Geneset
# We will use the PAM50 gene expression based subtype predictor for subtypes luminal A", "luminal B", "HER2-enriched", "and basal-like
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2667820/

PAM50_symbols <- c("UBE2C", "PTTG1", "MYBL2", "BIRC5", "CCNB1", "TYMS", "MELK", "CEP55", "KNTC2", "UBE2T", "RRM2", "CDC6", "ANLN", "ORC6L", "KIF2C", "EXO1", "CDCA1", "CENPF", "CCNE1", "MKI67", "CDC20", "MMP11", "GRB7", "ERBB2", "TMEM45B", "BAG1", "PGR", "MAPT", "NAT1", "GPR160", "FOXA1", "BLVRA", "CXXC5", "ESR1", "SLC39A6", "KRT17", "KRT5", "SFRP1", "BCL2", "KRT14", "MLPH", "MDM2,FGFR4", "MYC", "MIA", "FOXC1", "ACTR3B", "PHGDH", "CDH3", "EGFR")
PAM50_ensembl <- c("ENSG00000133627", "ENSG00000011426", "ENSG00000107262", "ENSG00000171791", "ENSG00000089685", "ENSG00000106605", "ENSG00000134057", "ENSG00000105173", "ENSG00000117399", "ENSG00000094804", "ENSG00000062038", "ENSG00000117724", "ENSG00000138180", "ENSG00000171604", "ENSG00000146648", "ENSG00000141736", "ENSG00000091831", "ENSG00000174371", "ENSG00000160867", "ENSG00000129514", "ENSG00000054598", "ENSG00000173890", "ENSG00000141738", "ENSG00000142945", "ENSG00000186847", "ENSG00000128422", "ENSG00000186081", "ENSG00000277956", "ENSG00000276155", "ENSG00000186868", "ENSG00000135679", "ENSG00000165304", "ENSG00000261857", "ENSG00000148773", "ENSG00000115648", "ENSG00000275365", "ENSG00000099953", "ENSG00000101057", "ENSG00000136997", "ENSG00000171428", "ENSG00000082175", "ENSG00000092621", "ENSG00000164611", "ENSG00000171848", "ENSG00000104332", "ENSG00000141424", "ENSG00000151715", "ENSG00000176890", "ENSG00000175063", "ENSG00000077152")

# [optional] Map symbols - using the biomaRt package & ensembl dataset
# https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
PAM50_ensembl <- getBM(attributes='ensembl_gene_id', 
      filters = 'hgnc_symbol', 
      values = PAM50_symbols, 
      mart = ensembl)


# Note that not all the symbols are measured and need to be excluded
setdiff(PAM50_ensembl, colnames(brca_norm))
geneset <- intersect(PAM50_ensembl,colnames(brca_norm) )

# PCA also requires 

# Run Principal Component Analysis (PCA) for clustering
# The prcomp method uses SVD decomposition as opposed to eigen in princomp
# to cluster the samples by gene values, we must transform the matrix values

pca<- princomp(brca_norm[,geneset ])

# The scree plot shows that the variation is contained in the first to components (axes)
plot(pca, type='lines')

#category = "her2_status_by_ihc"  #not enriched
category = "er_status_by_ihc"  #top / bottom
#category =  "pr_status_by_ihc"    #top / bottom
fields <- unique(brca_pt[, category])
colors <- brewer.pal(length(fields), name="Accent")

col = sapply(match(brca_pt[substr(rownames(pca$scores), 0, 12), category], fields), function(i) colors[i] )
# sample IDs have additional suffix, pt ID is within first 12 characters
# map samples to pt information for given category & assign color

plot(pca$scores[,c(2,1)], pch=16, col=col, main=paste("BRCA PAM50 expression",category, sep="\n"), xlab="PC1", ylab="PC2")
legend("topleft", legend=fields, fill=colors , bty ="n", ncol=4, cex=0.5)

## See Survival Difference in ER, PR, and rest of population
brca_pt$death_days_to[brca_pt$death_days_to == "[Not Applicable]"] <- NA
brca_pt$death_days_to <- as.numeric(brca_pt$death_days_to)
brca_pt$last_contact_days_to[brca_pt$last_contact_days_to == "[Not Applicable]"] <- NA
brca_pt$last_contact_days_to <- as.numeric(brca_pt$last_contact_days_to)

status = brca_pt$vital_status
# normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death).
status[status=="Alive"] = 1; status[status=="Dead"] = 2; 

brca_surv <- data.frame(
            pt = brca_pt$bcr_patient_barcode, 
            time = brca_pt$death_days_to, 
            status = as.numeric(status),
            classify = brca_pt[, category])

# Alive & no event timing, assign last contact date
brca_surv$time <- as.numeric(apply(brca_surv, 1,
 function(event){
    if (is.na(event["time"]) & event["status"]==1){
        event["time"] = brca_pt[event["pt"], "last_contact_days_to"]
    }
    event["time"]
}))

library(survival)
brca_survfit_all <- survfit(Surv(time, status) ~ 1, data = brca_surv)
brca_survfit_cat <- survfit(Surv(time, status) ~ classify, data = brca_surv)

#plot(brca_survfit_all) #show overall survival
#plot(brca_survfit_cat) # and by category

ggsurvplot(
    fit = survfit(Surv(time, status) ~ classify, data = brca_surv), 
    xlab = "Days", 
    ylab = "Overall survival probability", 
    title=paste("Survival by Biomarker:", category, sep="\n"), palette=brewer.pal(length(fields), name="Accent"))
