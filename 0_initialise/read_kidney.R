# biopsy.seurat.Robj is at:
## https://data.humancellatlas.org/explore/projects/027c51c6-0719-469f-a7f5-640fe57cbece/expression-matrices

load("~/Documents/main_files/GCP/explain-dpm/uploads/script_DRAFT_cell2pixel/data/main_gene_transcripts/biopsy.seurat.Robj",verbose = T)
biopsy <- Seurat::UpdateSeuratObject(biopsy)
counts <- biopsy@assays$RNA@counts
colnames(counts) <- biopsy@active.ident
write.csv(counts,"~/Documents/main_files/GCP/explain-dpm/uploads/script_DRAFT_cell2pixel/data/main_gene_transcripts/wu_counts.csv",quote=F)
