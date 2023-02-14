library(tidyverse)
library(dendextend)
library(viridis)
library(ggpubr)
library(patchwork)

read_arraytsv <- function(fname) {
  read_tsv(fname, col_select = c(11), show_col_types = FALSE) %>%
    #mutate(chip_name = colnames(.)[2] %>% str_remove('\\..*')) %>%
    rename_with(~ str_remove(., '\\..*'), c(1) ) %>%
    mutate(across(everything(), ~ ifelse(. == 'NC', NA, str_count(., 'A'))))
}

sampletable <- read_tsv('sample_table.txt') %>%
	mutate(sample_id = paste0(Chip_No, '_', Chip_Pos),
				 Celltype = str_remove(Celltype, 's$'))

all.files <- list.files('data', 'processed-data.tsv', full.names = T, recursive = T)

use.files <- all.files[str_detect(all.files, paste(sampletable$sample_id, collapse = '|'))]

full.snp.data <- lapply(use.files, read_arraytsv) %>% bind_cols()


SNP.matrix <- full.snp.data %>%
  # remove all SNP probes with any NA values (wont work properly)
	# TODO ?? and with the same value across all samples (uninformative)
  filter(if_all(everything(), ~!is.na(.)))
colnames(SNP.matrix) <- sampletable[ match(colnames(SNP.matrix), sampletable$sample_id),]$SampleName

#large table takes considerably longer -> don't do this for general report
d <- dist(SNP.matrix %>% t(), method = 'manhattan')
hc <- hclust(d)

#plot(hc)
dd <- as.dendrogram(hc) 

dend.format.df <- sampletable %>%
  mutate(SampleGroup = factor(SampleGroup, levels = sort(unique(SampleGroup)))) %>%
  arrange(SampleGroup) %>%
  mutate(shape = ifelse(str_detect(Celltype, 'iPSC'), 15, 17),
         shape = ifelse(str_detect(Celltype, 'Fibroblast'), 16, shape),
         #shape = ifelse(Celltype == 'ESC', 18, shape),
  			 Celltype = ifelse(!Celltype %in% c('iPSC', 'Fibroblast'), 'Other', Celltype),
         col = viridis_pal(option='H')(length(unique(SampleGroup)))[match(SampleGroup, unique(SampleGroup))],
  			 col = ifelse(is.na(SampleGroup), 'grey90', col),
         raw_na = sapply(sample_id, function(n) ifelse(n %in% colnames(full.snp.data), full.snp.data[,n] %>% is.na() %>% sum(), NA))
  )
dend.format.df <- dend.format.df[match(labels(dd), dend.format.df$SampleName),] 

dend <- dd %>%
  set('labels_col', dend.format.df$col) %>%
  set('leaves_col', dend.format.df$col) %>%
  set('leaves_pch', dend.format.df$shape) 
#plot(dend)

dend.format.df <- dend.format.df %>%
  arrange(SampleGroup) 

# Needed to map proper names to legends
col_map = unique(dend.format.df$col)
names(col_map) <- unique(dend.format.df$SampleGroup)
shape_map = 15:17 #unique(dend.format.df$shape)
names(shape_map) <- c('iPSC', 'Fibroblast', 'Other') #unique(dend.format.df$Celltype))

gg1 <- dend %>%
  raise.dendrogram (max(d)/50) %>% 
  set('labels_cex', .8) %>%
  set('leaves_cex', 2) %>%
  set('branches_lwd', .5) %>%
  as.ggdend() %>%
  ggplot(offset_labels = -max(d)/25, horiz = F) + 
  # horiz ylim(max(get_branches_heights(dend)), -500000) + 
  # vert
  ylim(-max(d)/2, NA) + 
  theme(legend.position = "bottom")

# make legend
gg2 <- dend.format.df %>%
  arrange(SampleGroup) %>%
  ggplot(aes(x = SampleName, y= 1, col = SampleGroup, shape = Celltype)) + 
  geom_point() + 
  scale_color_manual(values = col_map, guide = guide_legend(direction = 'horizontal',title.position = 'top', ncol = 10, byrow=T))  + 
  scale_shape_manual(values = shape_map, guide = guide_legend(direction = 'horizontal', title.position = 'top')) +
  theme(legend.box = "vertical", legend.text = element_text(size = 8), legend.title = element_text(size = 10))

gg <- gg1 + as_ggplot(get_legend(gg2)) + plot_layout(ncol = 1, heights = c(4, 1))

ggsave('clustered_samples.png', gg, height = 10, width = 12)

# #####
# 
# cl_samplesheet <- read_tsv('clustering_samplesheet.tsv') %>%
#   dplyr::rename(chip_no = `Chip/Sentrix Barcode (L&B)`,
#                 chip_pos = `SentrixPosition (L&B)`,
#                 sample_name = `Sample_ID (CORE)`,
#                 cellline_base = `Corresponding primary cells or iPSC line`,
#                 cell_type = `Cell type`
#   ) %>%
#   mutate(chip_name = paste0(chip_no, '_', chip_pos),
#          sample_id = paste(sample_name, `DNA ID/ Barcode (CORE)`))
# 
# 
# SNP.matrix <- full.snp.data %>%
#   select(one_of(cl_samplesheet$chip_name)) %>%
#   # remove all SNP probes with any NA values (wont work properly) and with the same value across all samples (uninformative)
#   filter(if_all(everything(), ~!is.na(.)))
# colnames(SNP.matrix) <- cl_samplesheet[ match(colnames(SNP.matrix), cl_samplesheet$chip_name),]$sample_id
# 
# 
# d2 <- dist(SNP.matrix %>% t())
# hc2 <- hclust(d2)
# dd2 <- as.dendrogram(hc2) 
# 
# dend.format.df <- cl_samplesheet[ match(labels(dd2), cl_samplesheet$sample_id),] %>%
#   mutate(cell_type = str_remove(cell_type, 's$'),
#          shape = ifelse(cell_type == 'iPSC', 15, 17),
#          shape = ifelse(cell_type == 'Fibroblast', 16, shape),
#          col = viridis_pal(option='H')(length(unique(cellline_base)))[match(cellline_base, unique(cellline_base))]
#          )
# 
# col_map = unique(dend.format.df$col)
# names(col_map) <- unique(dend.format.df$cellline_base)
# shape_map = unique(dend.format.df$shape)
# names(shape_map) <- unique(dend.format.df$cell_type)
# 
# gg1 <- dd2 %>%
#   raise.dendrogram (max(d2)/50) %>% 
#   set('labels_col', dend.format.df$col) %>%
#   set('leaves_col', dend.format.df$col) %>%
#   set('leaves_pch', dend.format.df$shape) %>%
#   set('labels_cex', 1) %>%
#   set('leaves_cex', 3) %>%
#   set('branches_lwd', .75) %>%
#   as.ggdend() %>%
#   ggplot(offset_labels = -max(d2)/25, horiz = F ) + 
#   ylim(-max(d2)/1.5, NA) + 
#   #scale_col_manual(values = col_map, name = 'Base cellline') + 
#   theme(legend.position = "bottom")
# # make legend
# gg2 <- ggplot(dend.format.df, aes(x = sample_name, y= 1, col = ReferenceSample, shape = Celltype)) + 
#   geom_point() + 
#   scale_color_manual(values = col_map, guide = guide_legend(direction = 'horizontal',title.position = 'top'))  + 
#   scale_shape_manual(values = shape_map, guide = guide_legend(direction = 'horizontal',title.position = 'top')) +
#   theme(legend.box = "horizontal")
# 
# gg <- gg1 + as_ggplot(get_legend(gg2)) + plot_layout(ncol = 1, heights = c(7, 1))
# 
# # D0497, D0495
# # --> I switched name & corresponding line
# 
# # 11395OC -> 11935OC

