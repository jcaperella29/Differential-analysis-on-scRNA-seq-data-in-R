library(Seurat)
library(SeuratData)
InstallData("ifnb")
ifnb <-LoadData("ifnb")

#Normalize the data
ifnb <- NormalizeData(ifnb)

# Find DE features between CD16 Mono and CD1 Mono
Idents(ifnb) <- "seurat_annotations"
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono")

# view results
head(monocyte.de.markers)
monocyte.de.markers_hits <- monocyte.de.markers[monocyte.de.markers$p_val_adj<=.05,] 

hits<- row.names(monocyte.de.markers_hits)
write.csv(monocyte.de.markers_hits,"C:/Users/ccape/OneDrive/Attachments/Desktop/SCRNA/monocyte_results.csv")

write.table(hits,"C:/Users/ccape/OneDrive/Attachments/Desktop/SCRNA/monocyte_hits.txt")
# Find differentially expressed features between CD14+ Monocytes and all other cells, only
# search for positive markers
monocyte.de.markers_postive <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = NULL, only.pos = TRUE)
# view results
head(monocyte.de.markers_postive)
monocyte.de.markers_postive_hits<- monocyte.de.markers_postive[monocyte.de.markers_up_reg$p_val_adj<=.05,] 

postive_hits<- row.names(monocyte.de.markers_postive)
write.csv(monocyte.de.markers_postive_hits,"C:/Users/ccape/OneDrive/Attachments/Desktop/SCRNA/monocyte_results_positve.csv")

write.table(postive_hits,"C:/Users/ccape/OneDrive/Attachments/Desktop/SCRNA/monocyte_positive_markers.txt")



#Perform DE analysis within the same cell type across conditions

ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
head(mono.de, n = 10)


#Perform DE analysis after pseudobulking

# load the inferred sample IDs of each cell
ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)

# rename the cell IDs by substituting the '-' into '.'
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)

# only keep the cells with high-confidence sample ID
info <- info[grep(pattern = "SNG", x = info$BEST), ]

# remove cells with duplicated IDs in both ctrl and stim groups
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]

# now add the sample IDs to ifnb 
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
ifnb <- AddMetaData(ifnb, metadata = info)

# remove cells without donor IDs
ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
ifnb <- subset(ifnb, subset = donor_id != "unknown")

# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))

pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
Idents(pseudo_ifnb) <- "celltype.stim"
#bulk DE analysis 
bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "CD14 Mono_STIM", 
                            ident.2 = "CD14 Mono_CTRL",
                            test.use = "DESeq2")
head(bulk.mono.de, n = 15)


# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.mono.de) <- paste0(names(bulk.mono.de), ".bulk")
bulk.mono.de$gene <- rownames(bulk.mono.de)

names(mono.de) <- paste0(names(mono.de), ".sc")
mono.de$gene <- rownames(mono.de)

merge_dat <- merge(mono.de, bulk.mono.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                 merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]

print(paste0('# Common: ',length(common)))


print(paste0('# Only in single-cell: ',length(only_sc)))

print(paste0('# Only in bulk: ',length(only_bulk)))

#create a new column to annotate sample-condition-celltype in the single-cell dataset
ifnb$donor_id.stim <- paste0(ifnb$stim, "-", ifnb$donor_id)

# generate violin plot 
Idents(ifnb) <- "celltype.stim"
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')])


VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim") 


VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1) 

print(merge_dat[merge_dat$gene%in%c('SRGN','HLA-DRA'),c('gene','p_val.sc','p_val.bulk')])

VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim") 

VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1) 



write.table(common,"C:/Users/ccape/OneDrive/Attachments/Desktop/SCRNA/same_cell_type_hits")

#orginal code came from https://satijalab.org/seurat/articles/de_vignette. Lines were added to help with reporting.
