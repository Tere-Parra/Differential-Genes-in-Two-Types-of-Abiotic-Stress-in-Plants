#librerias
library(dplyr)
library(DESeq2)
library(data.table)
library(tidyr)
library(gridExtra)
library(tidyverse)
library(tidyr)



setwd("C:/Users/Teresita Parra/Pictures/Abiotic_Stress/DEG")
#################################################################
##################   METADATOS ########################
######################################################3
# Realizar la tabla con con los metadatos y condiciones, para deseq2, es 
#necesario tener una columna con las muestras y otra con la condicion

sampledata <- data.frame(
  samples=c("WS.Car1_Read_Count", "WS.Car2_Read_Count",       
            "WS.Car3_Read_Count", "WS.Cleo1_Read_Count" ,      
            "WS.Cleo2_Read_Count" , "WS.Cleo3_Read_Count",
            "WS.HL.Car1.bis_Read_Count",
            "WSmasHL.Car2_Read_Count", 
            "WSmasHL.Car3_Read_Count", 
            "WSmasHL.Cleo1_Read_Count", 
            "WSmasHL.Cleo2_Read_Count",
            "WSmasHL.Cleo3_Read_Count"
            ),
  
  Condition=c( "Water Stress", "Water Stress", "Water Stress","Water Stress",
               "Water Stress","Water Stress", "Water Stress and high irradiance",
               "Water Stress and high irradiance","Water Stress and high irradiance",
               "Water Stress and high irradiance","Water Stress and high irradiance",
               "Water Stress and high irradiance")
  
)



    
   
   
head(sampledata)

#################################################################
##################   conteos ########################
######################################################

data <- read.table("counts.txt", header=TRUE, sep="\t")

#Delete columns we will not use
colnames(data)

data <- filter( data, Type=="protein_coding")
#12 658 coding genes


data <- select(data, -Feature_GID, -Feature_TID, -Type, -Gene_Symbol,
               -Gene_Synonym, -Protein_ID, -Product)

data <- data[, 1:49]

counts <- data %>% 
         gather(key = 'samples', value = 'counts', -Entrez_Gene_ID) %>% 
          merge(sampledata,  by="samples") %>%
             select(1,2,3) %>%
          mutate(samples = gsub('\\_', '.', samples)) 


  
counts <- counts %>%
         spread(key = 'samples', value = 'counts') 
         

counts <- na.omit(counts)

counts <- counts %>%
        column_to_rownames(var="Entrez_Gene_ID")



#apply changes
sampledata <- sampledata %>%
  mutate(samples = gsub('\\_', '.', samples)) 

save(sampledata, counts, file="Entrada_Deseq2.RData")
############################################################################

###############    ENTRADA A DESEQ2 ###########################

###################################################
load("Entrada_Deseq2.RData")
# making the rownames and column names identical
all(rownames(counts) %in% colnames(sampledata))
all(rownames(counts) == colnames(sampledata))


# create dds
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampledata,
                              design = ~ Condition) 

#Eliminar 0 counts

dds10 <- dds[rowSums(counts(dds)) >= 10 ]
nrow(dds10) # 10,613 genes

##########################################################################
##################### b) Estimar factores estabilizadores #########

dds_wt <- estimateSizeFactors(dds10)

##########################################################################
##################### c) Estimar factores estabilizadores #########

rld<-rlog(dds_wt, blind = TRUE)
head(assay(rld))


##################################################################
################# Heatmap y evaluacion de expresion genica #######

rld_mat_wt <- assay(rld)
vsd_cor_wt <- cor(rld_mat_wt)

#guardar matriz
write.table(vsd_cor_wt, file = "rld_cor_abiotic.tsv", sep = "\\t", row.names = TRUE, col.names = TRUE)


#graficar heatmap
library(pheatmap)
# Preparar los datos
data_for_heatmap <- as.matrix(vsd_cor_wt)

# Convert tissueType a un vector de caracter
annotation_row <- as.character(sampledata$Condition)

# Añadir espacios entre palabras 
annotation_row_with_spaces <- paste(" ", annotation_row, " ")

# Graficar el heatmap usando la librearia pheatmap 
pheatmap(data_for_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main="Abiotic Stress",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = FALSE,
         show_colnames = TRUE,
         row_names_side = "left",
         annotation_colors = "black",
         annotation_names_row = FALSE,
         labels_row = annotation_row_with_spaces,
         fontsize_row = 8,     
         fontsize_col = 12,    
         angle_col = "45")      

##################################################################
##############  graficar pca y boxplot ###########################################

boxplot(counts, outline=FALSE, main="Antes de la estabilización", xaxt="n")
boxplot(rld_mat_wt, outline=FALSE, main="Después de la estabilización", xaxt="n")

#PCA
plotPCA(rld, intgroup = "Condition")



save(dds,dds10, dds_wt, rld, file="rlds_Abiotic_Stress.RData")

####################################################################
##################### PASO 4 Analisis de Expresion Diferencial
load("rlds_Abiotic_Stress.RData")

# a) Utilizar la funcion DESEQ
dds <- DESeq(dds10)
head(dds)
## estimating size factors (normalizacion)
## estimating dispersions  (σ2)
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing

#b) Establecer las condiciones para el DEG
#- Brev Flor antesis vs  Long Flor Antesis 

res <- results(dds, contrast=c("Condition","Water Stress",
                               "Water Stress and high irradiance"))
res

summary(res)

#up 410 genes
#down 458 genes 

#crear una tabla con los genes up and down regulated
resSig <- subset(res, res$padj < 0.05 )


#Un total de upregulated genes de 1194 y down_regulated genes 72

#c) Generación de graficas 
#- Volcano plot simple
#- MAT


#ordenamos antes res de menor a mayor
res<-res[order(res$padj),]

res_table<- as.data.frame(res)

resig_table<- as.data.frame(resSig)
#guardar res
write.csv( as.data.frame(res), file="resultsDEG_AbioticStress.csv" )

#guardar regulates genes
write.csv( as.data.frame(resSig), file="Down_Up_AbioticStress.csv" )

save(res, dds,resSig, file="DEG_AbioticStress.RData")


sum( res$padj < 0.05, na.rm=TRUE )
#695  genes 


############################################################################
############################### gráficas ####################

##############################################################
############### volcano plot ##################
#reset 
par(mfrow=c(1,1))

# volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Abiotic Stress", xlim=c(-9,9)))

# azul si padj<0.05, rojo si log2FC>2 y padj<0.05)
with(subset(res, padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<0.05 & abs(log2FoldChange)>2), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


######## Heatmap ###
library(dplyr)
library(RColorBrewer)
library(pheatmap)

#1 ordenamos el DEG por el valor padj
res_ordered <- res[order(res$padj), ]
resSig <- subset(res_ordered, res$padj < 0.05 ) 


#2. Extraemos los valores normalizados
norm.data <- (counts(dds, normalized=TRUE))
norm_sig <- norm.data[rownames(resSig),]
head(norm_sig)

write.table(norm.data, "norm_AbioticStress.csv")
# 3. añadimos las anotaciones
annotation_row <- as.character(sampledata$Condition)

# Añadir espacios entre palabras 
annotation_row_with_spaces <- paste(" ", annotation_row, " ")
head(annotation_row_with_spaces )
### Color del heatmap
heat.colors <- brewer.pal(9, "YlOrRd")

#anotaciones
annotation <- data.frame(sampletype=sampledata[,'Condition'], 
                         row.names=rownames(sampledata))

# colnames(mat) <- str_sub(colnames(mat), 1, -3)
 rownames(annotation) <- colnames(norm.data)
### corremos el heatmap
pheatmap(norm_sig, main="Abiotic Stress", color = heat.colors, cluster_rows = T,
         show_rownames=F,labels_row= annotation_row_with_spaces, border_color=NA, scale="row",
         annotation= annotation, fontsize_row = 8, fontsize_col = 12, angle_col = "45")


save(norm.data, file="NormData.RData")
