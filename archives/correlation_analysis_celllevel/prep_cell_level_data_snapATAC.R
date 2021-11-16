library(tidyverse)
library(SnapATAC)
library(viridis)
library(GenomicRanges)
library(rtracklayer)
library(tictoc)
# library(leiden)
# library(foreach)
# library(doParallel)

## testing snapATAC pipeline with one sample

snap_path <- "/cndd2/fangming/projects/scf_enhancers/data/atac_raw/"
# barcode_path <- "/cndd2/fangming/projects/neuron_primer/snap_files/"

selected_sample <- list.files(snap_path) %>% str_subset("\\.snap$") %>% str_replace("\\.snap$", "")
snap_files <- paste0(snap_path, selected_sample, ".snap")

# barcodes_files <- paste0(barcode_path, selected_sample, ".csv")

x.sp <- createSnap(
  file = snap_files[1],
  sample = selected_sample[1],
  num.cores = 18
)

# x.sp.ls <- lapply(seq_along(snap_files), function(i){
#   createSnap(
#     file = snap_files[i],
#     sample = selected_sample[i],
#     num.cores = 20
#   )
# })
# names(x.sp.ls) <- selected_sample
# 
# ## barcode selection
# # select high-quality barcodes based on two criteria:
# # 1. number of unique fragments
# # 2. fragments in promoter ratio
# barcodes.ls <- lapply(seq_along(snap_files), function(i){
#   barcodes <- read.csv(
#     barcodes_files[i],
#     head = TRUE
#   )
#   
#   # remove first row with "NO_BARCODE"
#   barcodes <- barcodes[2:nrow(barcodes),]
#   # barcodes$logUMI <- log(barcodes$passed_filters + 1, 10)
#   # barcodes$promoter_ratio <- (barcodes$promoter_region_fragments + 1) / (barcodes$passed_filters + 1)
#   barcodes
# })
# 
# qc_plots <- lapply(seq_along(snap_files), function(i){
#   p1 <- ggplot(barcodes.ls[[i]], aes(x = logUMI, y = promoter_ratio)) +
#     geom_point(size = 0.1, color = "grey") +
#     theme_classic() +
#     ggtitle(selected_sample[i]) +
#     ylim(0, 1) +
#     xlim(0, 6) +
#     labs(x = "log10(UMI)", y = "promoter ratio")
#   
#   p1
# })
# # qc_plots
# 
# x.sp.ls
# 
# # use filtering criteria in snapATAC's tutorial
# # cutoff.logUMI.low = c(3.5, 3.5);
# # cutoff.logUMI.high = c(5, 5);
# # cutoff.FRIP.low = c(0.4, 0.4);
# # cutoff.FRIP.high = c(0.8, 0.8);
# 
# # use 10x solution to filter barcode multiplets
# # barcodes.sel = barcodes[which(barcodes$is__cell_barcode == 1), ]
# 
# # use both
# # tmp <- barcodes %>% as.tibble()
# # tmp %>% mutate(UMI = log(passed_filters+1, 10), promoter_ratio = (promoter_region_fragments + 1)/(passed_filters + 1)) %>% mutate(group = case_when(is__cell_barcode == 1 & UMI >=3 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6 ~ "yy", is__cell_barcode ==0 & (UMI >= 3 & UMI <=5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6) ~ "ny", is__cell_barcode ==1 & !(UMI >= 3 & UMI <= 5 &promoter_ratio >=0.15 & promoter_ratio <= 0.6) ~ "yn", TRUE ~ "nn"))-> tmp2
# # tmp2 %>% ggplot(aes(UMI, promoter_ratio)) + geom_point(size = 0.1, aes(color = group)) + theme_classic()
# 
# # rownames(barcodes.sel) <- barcodes.sel$barcode
# # x.sp <- x.sp[which(x.sp@barcode %in% barcodes.sel$barcode), ]
# # x.sp@metaData <- barcodes.sel[x.sp@barcode, ]
# # x.sp
# 
# barcodes.ls <- lapply(seq(snap_files), function(i){
#   barcodes <- barcodes.ls[[i]]
#   # use snapATAC recommended filtering criteria
#   # idx <- which(
#   #   barcodes$logUMI >= cutoff.logUMI.low[i] & 
#   #     barcodes$logUMI <= cutoff.logUMI.high[i] & 
#   #     barcodes$promoter_ratio >= cutoff.FRIP.low[i] &
#   #     barcodes$promoter_ratio <= cutoff.FRIP.high[i]
#   # )
#   
#   # use 10x solution to filter barcode multiplets
#   idx <- which(barcodes$is__cell_barcode ==1)
#   barcodes[idx,]
# })

# x.sp.ls <- lapply(seq(snap_files), function(i){
#   barcodes <- barcodes.ls[[i]]
#   x.sp <- x.sp.ls[[i]]
#   barcode.shared <- intersect(x.sp@barcode, barcodes$barcode)
#   x.sp <- x.sp[match(barcode.shared, x.sp@barcode),]
#   barcodes <- barcodes[match(barcode.shared, barcodes$barcode),]
#   x.sp@metaData <- barcodes
#   x.sp
# })
# names(x.sp.ls) <- selected_sample
# x.sp.ls
# 
# # combine the snap objects
# x.sp <- Reduce(snapRbind, x.sp.ls)
# x.sp@metaData["sample"] <- x.sp@sample
# x.sp

# save the snap object
saveRDS(x.sp, "/cndd2/fangming/projects/scf_enhancers/data/raw_atac/all_samples_merged_snap_files.rds")

# ## run dimensionality reduction with diffusion maps
# # to overcome the computational cost, use the Nystrom landmark diffusion maps algorithm
# landmarks_size <- 30000
# 
# tic("sampling landmark")
# # density-based sampling
# row_covs_dens <- density(
#   x = x.sp@metaData[, "logUMI"],
#   bw = 'nrd',
#   adjust = 1
# )
# sampling_prob <- 1 / (approx(x = row_covs_dens$x, y = row_covs_dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps)
# set.seed(666)
# idx_landmark_ds <- sort(sample(x = seq(nrow(x.sp)), size = landmarks_size, prob = sampling_prob))
# 
# # split the snap object
# x.sp <- x.sp[idx_landmark_ds, ]
# x.query.sp <- x.sp[-idx_landmark_ds, ]
# landmark.used.sample <- names(which(table(x.sp@metaData$sample)>1))
# query.used.sample <- names(which(table(x.query.sp@metaData$sample)>1))
# #used.sample <- intersect(landmark.used.sample, query.used.sample)
# x.sp <- x.sp[which(x.sp@metaData$sample %in% landmark.used.sample), ]
# x.query.sp <- x.query.sp[which(x.query.sp@metaData$sample %in% query.used.sample), ]
# toc()

## add cell-by-bin matrix
tic("addBmatToSnap")
x.sp <- addBmatToSnap(x.sp, do.par = TRUE, num.cores = 18, bin.size = 5000);
## Matrix binarization
# convert the cell-by-bin count matrix to binary matrix
# remove the top 0.1% values in the count matrix and then convert the remaining non-zero values to 1 
x.sp <- makeBinary(x.sp, mat = "bmat", outlier.filter = 0.001);
toc()

## Bin filtering
# filter out bins overlapping with the ENCODE blacklist
# version 2 from: https://github.com/Boyle-Lab/Blacklist/tree/master/lists
tic("filterBins")
black_list <- read_tsv("/cndd/junhao/genomes/encode_blacklist/mm10-blacklist.v2.bed.gz",
                       col_names = c("chrom", "start", "end", "type"))
black_list_gr <- GRanges(black_list$chrom, IRanges(black_list$start, black_list$end))
idy <- queryHits(findOverlaps(x.sp@feature, black_list_gr))
if(length(idy) > 0){
  x.sp <- x.sp[, -idy, mat = "bmat"]
}

# remove unwanted chromosomes
chr_exclude <- seqlevels(x.sp@feature)[grep("_|chrM|chrEBV", seqlevels(x.sp@feature))]
idy <- grep(paste(chr_exclude, collapse = "|"), x.sp@feature)
if(length(idy) > 0){
  x.sp <- x.sp[, -idy, mat = "bmat"]
}

# remove bins with top 5% bin coverage that overlap with invariant features such as promoters of the house keeping genes
bin_cov <- log10(Matrix::colSums(x.sp@bmat) + 1)
bin_cutoff <- quantile(bin_cov[bin_cov > 0], 0.95)
idy <- which(bin_cov <= bin_cutoff & bin_cov > 0)
x.sp <- x.sp[, idy, mat="bmat"];

pdf("all_samples_random_lankmark_5k_bin_coverage.pdf")
hist(
  bin_cov[bin_cov > 0],
  xlab = "log10(bin cov)",
  main = "log10(Bin Cov)",
  col = "lightblue",
  xlim = c(0,5)
)
dev.off()

# optional but highly recommended: remove any cells of bin coverage less than 1000
idx = which(Matrix::rowSums(x.sp@bmat) > 1000)
x.sp = x.sp[idx,]

toc()


# run diffusion maps on the landmark cells
tic("runDiffusionMaps")
x.sp <- runDiffusionMaps(
  obj = x.sp,
  input.mat = "bmat",
  num.eigs = 50
)
# x.sp@metaData$landmark <- 1
toc()

## determine significant components
# ad hoc, select the number of dimensions in which the scatter plot starts looking like a blob

# pdf("all_samples_random_lankmark_5k_bin_plot_dimResuce.pdf")
# plotDimReductElbow(
#   obj = x.sp,
#   point.size = 1.5,
#   point.shape = 19,
#   point.color = "red",
#   point.alpha = 1,
#   pdf.file.name = NULL,
#   pdf.height = 7,
#   pdf.width = 7,
#   labs.title = "PCA Elbow plot",
#   labs.subtitle = NULL
# )
# 
# plotDimReductPW(
#   obj = x.sp, 
#   eigs.dims = 1:50,
#   point.size = 0.3,
#   point.color = "grey",
#   point.shape = 19,
#   point.alpha = 0.6,
#   # down.sample = 5000,
#   pdf.file.name = NULL,
#   pdf.height = 7,
#   pdf.width = 7 
# )
# dev.off()
# 
# saveRDS(x.sp, "snap_object_all_samples_landmark_subset_post_diffusionMap.rds")
# saveRDS(x.query.sp, "snap_object_all_samples_query_subset.rds")

# # project query cells to landmarks by batch/sample
# tic("extension")
# 
# query_sample_list <- unique(x.query.sp@sample)
# num_query_files <- length(query_sample_list)
# 
# run_projection <- function(query_sample, num_core = 12, bin_size = 5000){
#   
#   message(paste0("Projecting cells from sample: ", query_sample))
#   
#   message("Subseting...")
#   x.query.sub <-
#     x.query.sp[which(x.query.sp@sample == as.character(query_sample)),]
#   
#   message("Adding Bmat...")
#   x.query.sub <- addBmatToSnap(
#     x.query.sub,
#     do.par = TRUE,
#     num.cores = num_core,
#     bin.size = bin_size
#   )
#   
#   message("Making binary...")
#   x.query.sub <- makeBinary(x.query.sub, mat = "bmat", outlier.filter = 0.001)
#   
#   message("Filtering...")
#   idy <- unique(queryHits(findOverlaps(
#     x.query.sub@feature, x.sp@feature
#   )))
#   
#   x.query.sub <- x.query.sub[, idy, mat = "bmat"]
#   
#   idx <- which(Matrix::rowSums(x.query.sub@bmat) > 1000)
#   x.query.sub <- x.query.sub[idx,]
# 
#   message("Projecting...")
#   x.query.sub <- runDiffusionMapsExtension(obj1 = x.sp,
#                                            obj2 = x.query.sub,
#                                            input.mat = "bmat")
#   
#   x.query.sub@metaData$landmark <- 0
#   
#   message("Removing Bmat...")
#   x.query.sub <- rmBmatFromSnap(x.query.sub)
#   return(x.query.sub)
# }
# 
# x.query.ls <- map(query_sample_list, run_projection)
# toc()
# 
# tic("Merge snap files")
# # combine landmark and query cells
# x.sp <- rmBmatFromSnap(x.sp)
# x.query.sp <- Reduce(snapRbind, x.query.ls)
# x.sp <- snapRbind(x.sp, x.query.sp)
# x.sp <- x.sp[order(x.sp@metaData[, "sample"])]
# toc()
# 
# saveRDS(x.sp, "snap_object_all_samples_post_diffusionMap.rds")
# # x.sp <- readRDS("snap_object_all_samples_post_diffusionMap.rds")
# write_tsv(x.sp@metaData %>% as_tibble(), "snap_all_samples_post_diffusionMap_metadata.txt")
# 
# ## determine significant components, continued...
# selected_dims <- 30
# 
# ## graph-based clustering
# tic("runKNN")
# x.sp <- runKNN(
#   obj = x.sp,
#   eigs.dims = 1:selected_dims,
#   k = 15
# )
# toc()
# 
# # use louvain by igraph for now, try installing leiden later...
# tic("Run Louvain clustering")
# x.sp <- runCluster(
#   obj = x.sp,
#   tmp.folder = tempdir(),
#   louvain.lib = "R-igraph",
#   seed.use = 666,
#   resolution = 1
# )
# x.sp@metaData$cluster <- x.sp@cluster
# toc()
# 
# # use leiden clustering
# # tic("Run leiden clustering")
# # x.sp <- runCluster(
# #   obj = x.sp,
# #   tmp.folder = tempdir(),
# #   louvain.lib = "leiden",
# #   seed.use = 666,
# #   resolution = 1
# # )
# # x.sp@metaData$cluster <- x.sp@cluster
# # toc()

## visualization

x.sp <- runViz(
  obj = x.sp,
  tmp.folder = tempdir(),
  dims = 2,
  eigs.dims = 1:50,
  # method = "Rtsne",
  method = "umap",
  seed.use = 666
)

# vis_method <- "tsne"
vis_method <- "umap"


write_tsv(x.sp@barcode %>% as_tibble(), "/cndd2/fangming/projects/neuron_primer/snap3_barcode.txt")
write_tsv(x.sp@sample %>% as_tibble(), "/cndd2/fangming/projects/neuron_primer/snap3_sample.txt")
write_tsv(x.sp@umap %>% as_tibble(), "/cndd2/fangming/projects/neuron_primer/snap3_umap.txt")

plotViz(obj=x.sp, 
        method="umap",
        )



# 
# # add disease states and brain region info into metadata
# 
# metaData_df <- x.sp@metaData %>% 
#   as.tibble() %>% 
#   mutate(disease = case_when(str_detect(sample, "ALS") ~ "ALS",
#                              str_detect(sample, "FT[DP]") ~ "FTD",
#                              str_detect(sample, "control") ~ "control"),
#          brain_region = case_when(str_detect(sample, "[mM]CX") ~ "MCX",
#                                   str_detect(sample, "FCX") ~ "FCX"))
# x.sp@metaData$disease <- metaData_df$disease
# x.sp@metaData$brain_region <- metaData_df$brain_region
# 
# x.sp@metaData$umap_1 <- x.sp@umap[, "umap-1"]
# x.sp@metaData$umap_2 <- x.sp@umap[, "umap-2"]
# 
# # saveRDS(x.sp, "snap_object_all_samples_post_umap.rds")
# # x.sp <- readRDS("snap_object_all_samples_post_umap.rds")
# write_tsv(x.sp@metaData %>% as_tibble(), "snap_all_samples_post_umap_metadata.txt")
# 
# # filter outliers on umap space?(check what are those later...)
# idx = which((x.sp@umap[, "umap-1"] < 15) &
#               (x.sp@umap[, "umap-1"] > -15) &
#               (x.sp@umap[, "umap-2"] < 6.5) &
#               (x.sp@umap[, "umap-2"] > -6.5))
# x.sp = x.sp[idx,]
# 
# pdf("all_samples_cluster_and_qc.pdf", width = 8, height = 6)
# par(mfrow = c(2, 2))
# 
# plotViz(
#   obj = x.sp,
#   method = vis_method, 
#   main = "Cluster",
#   point.color = x.sp@cluster, 
#   point.size = 0.1,
#   point.shape = 19,
#   point.alpha = 0.8,
#   text.add = TRUE,
#   text.size = 1,
#   text.color = "black",
#   text.halo.add = TRUE,
#   text.halo.color = "white",
#   text.halo.width = 0.2,
#   down.sample = 10000,
#   legend.add = FALSE
# )
# 
# plotFeatureSingle(
#   obj = x.sp,
#   feature.value = x.sp@metaData[, "logUMI"],
#   method = vis_method, 
#   main = "Read Depth",
#   point.size = 0.2,
#   point.shape = 19,
#   down.sample = 10000,
#   quantiles = c(0.01, 0.99)
# ) 
# plotFeatureSingle(
#   obj = x.sp,
#   feature.value = x.sp@metaData$promoter_ratio,
#   method = vis_method, 
#   main = "FRiP",
#   point.size = 0.2,
#   point.shape = 19,
#   down.sample = 10000,
#   quantiles = c(0.01, 0.99)
# )
# plotFeatureSingle(
#   obj = x.sp,
#   feature.value = x.sp@metaData$duplicate / x.sp@metaData$total,
#   method = vis_method, 
#   main = "Duplicate rate",
#   point.size = 0.2,
#   point.shape = 19,
#   down.sample = 10000,
#   quantiles = c(0.01, 0.99)
# )
# 
# dev.off()
# 
# pdf("all_samples_umap_by_metadata.pdf", width = 8, height = 6)
# par(mfrow = c(2, 2))
# plotViz(
#   obj = x.sp,
#   method = vis_method,
#   main = "Sample",
#   point.size = 0.2,
#   point.shape = 19,
#   point.color = x.sp@sample,
#   text.add = FALSE,
#   text.size = 1.2,
#   text.color = "black",
#   down.sample = 10000,
#   legend.add = FALSE
# )
# 
# plotViz(
#   obj = x.sp,
#   method = vis_method,
#   main = "Disease",
#   point.size = 0.2,
#   point.shape = 19,
#   point.color = x.sp@metaData$disease,
#   text.add = FALSE,
#   text.size = 1.2,
#   text.color = "black",
#   down.sample = 10000,
#   legend.add = TRUE
# )
# 
# plotViz(
#   obj = x.sp,
#   method = vis_method,
#   main = "Brain region",
#   point.size = 0.2,
#   point.shape = 19,
#   point.color = x.sp@metaData$brain_region,
#   text.add = FALSE,
#   text.size = 1.2,
#   text.color = "black",
#   down.sample = 10000,
#   legend.add = TRUE
# )
# 
# plotViz(
#   obj = x.sp,
#   method = vis_method,
#   main = "Landmark",
#   point.size = 0.2,
#   point.shape = 19,
#   point.color = x.sp@metaData[, "landmark"],
#   text.add = FALSE,
#   text.size = 1.2,
#   text.color = "black",
#   down.sample = 10000,
#   legend.add = TRUE
# )
# dev.off()
# 
# ## Gene based annotation
# 
# genes_gr <- import.bed("gencode_V28_gene.bed")
# marker_genes <- c(
#   "SNAP25", "GAD2", "APOE",
#   "C1QB", "PVALB", "VIP",
#   "SST", "LAMP5", "SLC17A7",
#   "VWF", "MBP", "OLIG1",
#   "OLIG2", "PDGFRA", "VCAN",
#   "TBX18", "TLE4", "RORB",
#   "CUX2", "ADARB2", "LHX6",
#   "SOX6", "SOX10", "LRRK1",
#   "FOXP2", "PDZRN4", "TSHZ2"
# )
# 
# genes_sel_gr <- genes_gr[which(genes_gr$name %in% marker_genes)]
# 
# # re-add the cell-by-bin matrix to the snap object
# # x.sp <- addBmatToSnap(x.sp)
# # x.sp <- createGmatFromMat(
# #   obj = x.sp,
# #   input.mat = "bmat",
# #   genes = genes_sel_gr,
# #   do.par = TRUE,
# #   num.cores = 16
# # )
# 
# x.sp <- addGmatToSnap(x.sp, do.par = TRUE, num.cores = 12)
# 
# # saveRDS(x.sp, "snap_object_all_samples_post_addCellByGeneMat.rds")
# # x.sp <- readRDS("snap_object_all_samples_post_addCellByGeneMat.rds")
# 
# # subset gmat..
# x.sp@gmat <- x.sp@gmat[, marker_genes]
# 
# # normalize the cell-by-gene matrix
# # check other ways (TPM?) to normalize later...
# x.sp <- scaleCountMatrix(
#   obj = x.sp,
#   cov = x.sp@metaData$passed_filters + 1,
#   mat = "gmat",
#   method = "RPM"
# )
# 
# # subset...
# subset_size <- 10000
# row_covs_dens <- density(
#   x = x.sp@metaData[, "logUMI"],
#   bw = 'nrd',
#   adjust = 1
# )
# sampling_prob <- 1 / (approx(x = row_covs_dens$x, y = row_covs_dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps)
# set.seed(23)
# idx_subset_ds <- sort(sample(x = seq(nrow(x.sp)), size = subset_size, prob = sampling_prob))
# 
# # split the snap object
# x.sub.sp <- x.sp[idx_subset_ds, ]
# # 
# # x.sub.sp <- scaleCountMatrix(
# #   obj = x.sub.sp,
# #   cov = x.sub.sp@metaData$passed_filters + 1,
# #   mat = "gmat",
# #   method = "RPM"
# # )
# # smooth the cell-by-gene matrix
# x.sub.sp <- runMagic(
#   obj = x.sub.sp,
#   input.mat = "gmat",
#   step.size = 3
# )
# 
# par(mfrow = c(3, 3))
# 
# # for (i in 1:9){
# # for (i in 10:18){
# for (i in 19:27){
#   plotFeatureSingle(
#     obj = x.sub.sp,
#     feature.value = x.sub.sp@gmat[, marker_genes[i]],
#     method = vis_method,
#     main = marker_genes[i],
#     point.size = 0.1,
#     point.shape = 19,
#     down.sample = 10000,
#     quantiles = c(0.05, 0.95)
#   )
# }
# 
# 
# ## bar plot by sample, region, disease
# 
# meta_data <- x.sp@metaData %>% as_tibble()
# 
# tmp <- meta_data %>% count(cluster, disease)
# 
# p <- tmp %>% ggplot(aes(cluster, n))
# p + 
#   geom_bar(aes(fill = disease), stat = "identity", position = "stack") +
#   theme_bw(base_size = 10, base_family = "Helvetica") +
#   scale_fill_manual(values = c("#BEBEBE", "#E3191C", "#FFD700")) +
#   ylab("# cells") +
#   xlab("Cluster") +
#   coord_flip() +
#   theme(panel.grid.minor = element_blank())
# ggsave("atac_cell_count_by_cluster_color_by_disease.pdf", device = cairo_pdf(),
#        width = 6, height = 4)
# 
# p + 
#   geom_bar(aes(fill = disease), stat = "identity", position = "fill") +
#   theme_bw(base_size = 10, base_family = "Helvetica") +
#   scale_fill_manual(values = c("#BEBEBE", "#E3191C", "#FFD700")) +
#   ylab("# cells") +
#   xlab("Cluster") +
#   coord_flip() +
#   theme(panel.grid.minor = element_blank())
# ggsave("atac_cell_count_by_cluster_color_by_disease_percent.pdf", device = cairo_pdf(),
#        width = 6, height = 4)
# 
# tmp <- meta_data %>% count(cluster, brain_region)
# 
# p <- tmp %>% ggplot(aes(cluster, n))
# p + 
#   geom_bar(aes(fill = brain_region), stat = "identity", position = "stack") +
#   theme_bw(base_size = 10, base_family = "Helvetica") +
#   scale_fill_manual(values = c("#BEBEBE", "#E3191C", "#FFD700")) +
#   ylab("# cells") +
#   xlab("Cluster") +
#   coord_flip() +
#   theme(panel.grid.minor = element_blank())
# ggsave("atac_cell_count_by_cluster_color_by_region.pdf", device = cairo_pdf(),
#        width = 6, height = 4)
# 
# p + 
#   geom_bar(aes(fill = brain_region), stat = "identity", position = "fill") +
#   theme_bw(base_size = 10, base_family = "Helvetica") +
#   scale_fill_manual(values = c("#BEBEBE", "#E3191C", "#FFD700")) +
#   ylab("# cells") +
#   xlab("Cluster") +
#   coord_flip() +
#   theme(panel.grid.minor = element_blank())
# ggsave("atac_cell_count_by_cluster_color_by_region_percent.pdf", device = cairo_pdf(),
#        width = 6, height = 4)
# 
# sample_levels <- meta_data$sample %>% as.factor() %>% levels()
# sample_levels_order <- c(
#   sample_levels %>% str_subset("ALS"),
#   sample_levels %>% str_subset("control"),
#   sample_levels %>% str_subset("FT[DP]")
# )
# 
# tmp <- meta_data %>%
#   mutate(sample = factor(sample, levels = sample_levels_order)) %>%
#   count(cluster, sample)
# 
# p <- tmp %>% ggplot(aes(cluster, n))
# p + 
#   geom_bar(aes(fill = sample), stat = "identity", position = "stack") +
#   theme_bw(base_size = 10, base_family = "Helvetica") +
#   # scale_fill_manual(values = c("#BEBEBE", "#E3191C", "#FFD700")) +
#   ylab("# cells") +
#   xlab("Cluster") +
#   coord_flip() +
#   theme(panel.grid.minor = element_blank())
# 
# ggsave("atac_cell_count_by_cluster_color_by_sample.pdf", device = cairo_pdf(),
#        width = 8, height = 5)
# 
# 
# 