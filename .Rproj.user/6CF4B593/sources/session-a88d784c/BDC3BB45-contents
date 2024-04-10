# 0. 加载R包和函数
require(magrittr)
require(microeco)
require(ggplot2)
require(ggpubr)
source('libs/utils.r')
source('libs/fun-red-mp.r')

# 1. 读取文件
## 创建目录
outdir_01 <- "results/01.read-files"
if (!outdir_01 %>% dir.exists) dir.create(outdir_01, recursive = T)

## 读取ASV表
otu_table <- read.table("rawdata/asv-table-resampled-8000.tsv", header = T, sep = '\t')
rownames(otu_table) <- otu_table$OTUID
otu_table <- otu_table[,-1]
otu_table[1:3, 1:3] %>% print
otu_table %>% dim # 8916 * 30

## 读取物种分类信息
tax_table <- data.table::fread("rawdata/asv-taxonomy-silva-138.tsv", sep = '\t') %>% as.data.frame
rownames(tax_table) <- tax_table$`Feature ID`
tax_table <- tax_table[rownames(otu_table),]
tax_table %>% dim # 8916 * 3
tax_rst <- lapply(X = tax_table %>% nrow %>% seq, function(x){
    tax_res <- tax_table$Taxon[x] %>% strsplit(., split = "; ", fixed = T) %>% unlist
    if (tax_res %>% length < 7) tax_res <- c(tax_res, rep(NA, 7 - length(tax_res)))
    tax_res %<>% as.data.frame %>% t
    colnames(tax_res) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    rownames(tax_res) <- rownames(tax_table)[x]
    tax_res[1,] <- lapply(tax_res[1,], function(col) gsub("^D_\\d+__", "", col)) %>% unlist
    tax_res
}) %>% do.call('rbind', .)
tax_table <- cbind(tax_table, tax_rst)
tax_table <- tax_table[,-c(1, 2, 3)]
tax_table[1:3, 1:3] %>% print
tax_table %>% dim # 8916 * 7

## 读取KEGG拷贝数
kegg_copies <- data.table::fread("rawdata/ko-copies-for-each-asv.tsv.gz", 
                                 sep = '\t', header = T, nThread = 5) %>% as.data.frame()
rownames(kegg_copies) <- kegg_copies$OTUID; kegg_copies <- kegg_copies[,-1]
kegg_copies <- kegg_copies[rownames(otu_table),]
kegg_copies[1:4, 1:4] %>% print
kegg_copies %>% dim() # 8916 * 7309

## 读取分组表
gro_table <- xlsx::read.xlsx("rawdata/metadata-table.xlsx", sheetIndex = 1, as.data.frame = T)
rownames(gro_table) <- gro_table$ID
gro_table$Elevation <- as.character(gro_table$Elevation)
gro_table <- gro_table[colnames(otu_table),]
gro_table %>% head
gro_table %>% dim # 30 * 4

## 读取土壤参数
met_table <- xlsx::read.xlsx("rawdata/metadata-table.xlsx", sheetIndex = 2, as.data.frame = T)
rownames(met_table) <- met_table$ID; met_table <- met_table[,-1]
met_order <- c("pH", "MC", "ST", "CD", "TOC", "TN", "NH", "NO", "TP", "AP")
met_table <- met_table[colnames(otu_table), met_order]
met_table %>% head
met_table %>% dim # 30 * 10

## 检查样本顺序
if (length(unique(colnames(otu_table) == rownames(gro_table))) > 1 || 
    !unique(colnames(otu_table) == rownames(gro_table))) stop("sample order error!") 
if (length(unique(rownames(met_table) == rownames(gro_table))) > 1 || 
    !unique(rownames(met_table) == rownames(gro_table))) stop("sample order error!") 
if (length(unique(rownames(tax_table) == rownames(otu_table))) > 1 || 
    !unique(rownames(tax_table) == rownames(otu_table))) stop("asv order error!") 
if (length(unique(rownames(tax_table) == rownames(kegg_copies))) > 1 || 
    !unique(rownames(tax_table) == rownames(kegg_copies))) stop("asv order error!") 

## 设置可视化参数
args <- list(
    taxa_type_full_order = c('all', 'abundant','rare', 'hyper_rare'),
    taxa_type_part_order = c('abundant','rare', 'hyper_rare'),
    elev_order = c('2900','3102', '3194'),
    slop_order = c('NS','CS'),
    met_order = c("pH", "MC", "ST", "CD", "TOC", "TN", "NH", "NO", "TP", "AP"),
    taxa_type_full_color = RColorBrewer::brewer.pal(9, "Set1")[c(4,3,2,1)],
    taxa_type_part_color = RColorBrewer::brewer.pal(9, "Set1")[c(3,2,1)],
    elev_color = RColorBrewer::brewer.pal(8, "Set1")[1:3],
    slop_color = RColorBrewer::brewer.pal(8, "Set2")[1:2]
)

## 保存数据集
raw_data = list(otu = otu_table, tax = tax_table, gro = gro_table, met = met_table, cop = kegg_copies, args = args) 
saveRDS(object = raw_data, file = file.path(outdir_01, "16SrRNA_dataset.Rds"))

# 2. 拆分ASV数据集
## 创建目录
outdir_02 <- "results/02.split-data"
if (!outdir_02 %>% dir.exists) dir.create(outdir_02, recursive = T)

## 将ASV表拆分为3类: `hyper_rare`, `rare` 和 `abundant`
otu_table_trans <- (otu_table/sum(otu_table)) * 100
mean_r_abund <- lapply(X = otu_table_trans %>% nrow %>% seq, function(x){
    otu_table_trans[x,] %>% as.numeric %>% mean
}) %>% unlist()
mean_r_abund_df <- data.frame(OTUID = rownames(otu_table_trans), 
                              abund_perct = mean_r_abund)
mean_r_abund_df$type <- cut(mean_r_abund_df$abund_perct, 
                            breaks = c(-Inf, 0.0001, 0.001, Inf), label = c("hyper_rare", "rare", "abundant")) 
type_count_table <- table(mean_r_abund_df$type) %>% as.data.frame
otu_tables <- list(hyper_rare = otu_table[mean_r_abund_df[which(mean_r_abund_df$type == 'hyper_rare'),]$OTUID,] %>% t, 
                   rare = otu_table[mean_r_abund_df[which(mean_r_abund_df$type == 'rare'),]$OTUID,] %>% t, 
                   abundant = otu_table[mean_r_abund_df[which(mean_r_abund_df$type == 'abundant'),]$OTUID,] %>% t, 
                   all = otu_table %>% t)
saveRDS(object = otu_tables, file = file.path(outdir_02, "ranked-otu-tables.Rds"))

## 计算每个样本中每一类ASV的数量和相对丰度
perct_rst <- lapply(X = otu_tables$all %>% nrow %>% seq, function(x){
    sample <- rownames(otu_tables$all)[x]
    elev <- gro_table[which(gro_table$ID == sample),]$Elevation
    slop <- gro_table[which(gro_table$ID == sample),]$SlopeType   
    dat.abundant <- otu_tables$all[x, mean_r_abund_df[which(mean_r_abund_df$type == 'abundant'),]$OTUID]
    dat.rare <- otu_tables$all[x, mean_r_abund_df[which(mean_r_abund_df$type == 'rare'),]$OTUID]
    dat.hyper_rare <- otu_tables$all[x, mean_r_abund_df[which(mean_r_abund_df$type == 'hyper_rare'),]$OTUID]
    dat.abundant.1 <- dat.abundant; dat.rare.1 <- dat.rare; dat.hyper_rare.1 <- dat.hyper_rare
    dat.abundant[dat.abundant > 0] <- 1
    dat.rare[dat.rare > 0] <- 1
    dat.hyper_rare[dat.hyper_rare > 0] <- 1
    data.frame(sample = rep(sample, 3), elev = rep(elev, 3), slop = rep(slop, 3), 
               freq = c(dat.abundant %>% sum, dat.rare %>% sum, dat.hyper_rare %>% sum),
               abund = c(dat.abundant.1 %>% sum, dat.rare.1 %>% sum, dat.hyper_rare.1 %>% sum),
               type = c('abundant', 'rare', 'hyper_rare'))
}) %>% do.call('rbind', .) 
perct_rst$slop <- factor(perct_rst$slop, levels = args$slop_order)
p.perc.2900.freq <- ggplot(perct_rst[which(perct_rst$elev == '2900'),], aes(x = sample, y = freq, fill = type)) + 
        geom_bar(stat = 'identity', position = 'stack') + 
        ylab("ASV number") + facet_grid(. ~ slop, scales = 'free', space = 'free') + 
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
              axis.title.x = element_blank(), legend.position = 'right', legend.title = element_blank())
p.perc.2900.abud <- ggplot(perct_rst[which(perct_rst$elev == '2900'),], aes(x = sample, y = abund, fill = type)) + 
    geom_bar(stat = 'identity', position = 'fill') + 
    ylab("Relative abundance (%)") + scale_y_continuous(labels = function(x) x * 100) + 
    facet_grid(. ~ slop, scales = 'free', space = 'free') + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.title.x = element_blank(), legend.position = 'right', legend.title = element_blank())
p.perc.3102.freq <- ggplot(perct_rst[which(perct_rst$elev == '3102'),], aes(x = sample, y = freq, fill = type)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    ylab("ASV number") + facet_grid(. ~ slop, scales = 'free', space = 'free') + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.title.x = element_blank(), legend.position = 'right', legend.title = element_blank())
p.perc.3102.abud <- ggplot(perct_rst[which(perct_rst$elev == '3102'),], aes(x = sample, y = abund, fill = type)) + 
    geom_bar(stat = 'identity', position = 'fill') + 
    ylab("Relative abundance (%)") + scale_y_continuous(labels = function(x) x * 100) + 
    facet_grid(. ~ slop, scales = 'free', space = 'free') + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.title.x = element_blank(), legend.position = 'right', legend.title = element_blank())
p.perc.3194.freq <- ggplot(perct_rst[which(perct_rst$elev == '3194'),], aes(x = sample, y = freq, fill = type)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    ylab("ASV number") + facet_grid(. ~ slop, scales = 'free', space = 'free') + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.title.x = element_blank(), legend.position = 'right', legend.title = element_blank())
p.perc.3194.abud <- ggplot(perct_rst[which(perct_rst$elev == '3194'),], aes(x = sample, y = abund, fill = type)) + 
    geom_bar(stat = 'identity', position = 'fill') + 
    ylab("Relative abundance (%)") + scale_y_continuous(labels = function(x) x * 100) + 
    facet_grid(. ~ slop, scales = 'free', space = 'free') + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.title.x = element_blank(), legend.position = 'right', legend.title = element_blank())
(p.perc.freq <- ggpubr::ggarrange(p.perc.2900.freq, p.perc.3102.freq, p.perc.3194.freq, 
                                  ncol = 1, labels = c("A", "B", "C"), align = 'hv', common.legend = T, legend = 'right'))
(p.perc.abud <- ggpubr::ggarrange(p.perc.2900.abud, p.perc.3102.abud, p.perc.3194.abud, 
                                  ncol = 1, labels = c("A", "B", "C"), align = 'hv', common.legend = T, legend = 'right'))
saved <- savePDF(object = p.perc.freq, path = outdir_02, filename = "每一个样本中三类ASV的种类数量.pdf", width = 8.52, height = 12.89)
saved <- savePDF(object = p.perc.abud, path = outdir_02, filename = "每一个样本中三类ASV的相对丰度.pdf", width = 8.52, height = 12.89)

## 计算每一类ASV的分布情况
otu_occupies <- lapply(X = otu_tables %>% length %>% seq, function(x){
    otu_tab_abu <- otu_tables[[x]] %>% t
    otu_tab_bin <- otu_tables[[x]] %>% t
    otu_tab_bin[otu_tab_bin > 0] <- 1
    otu_nam <- names(otu_tables)[x]
    occupies <- rowSums(otu_tab_bin) %>% as.data.frame
    reads_numbers <- rowSums(otu_tab_abu) %>% as.data.frame
    dat <- data.frame(otu = rownames(occupies), occupies = occupies, reads_numbers = reads_numbers, 
                      type = rep(otu_nam, nrow(occupies)))
    colnames(dat) <- c('otu', 'occupies', 'reads_numbers', 'type')
    dat
}) %>% do.call('rbind', .)
otu_occupies$type <- factor(otu_occupies$type, levels = args$taxa_type_part_order)
type_count_table$Var1 <- factor(type_count_table$Var1, levels = args$taxa_type_part_order)
mean_r_abund_df$type <- factor(mean_r_abund_df$type, levels = args$taxa_type_part_order)
p1.1 <- ggplot(otu_occupies[which(otu_occupies$type != 'all'),], aes(x = occupies, y = log2(reads_numbers))) + 
    geom_point(aes(fill = type, color = type, size = reads_numbers), shape = 21,  alpha = 0.7) + 
    theme_gray() + xlab("Occupied sites") + ylab("log2(read number)") + 
    scale_fill_manual(values = args$taxa_type_part_color) + 
    scale_color_manual(values = args$taxa_type_part_color) + 
    theme(legend.position = 'none', axis.title = element_text(size = 12), axis.text = element_text(size = 10)) 
min(otu_occupies[which(otu_occupies$type == 'hyper_rare'),]$occupies)
max(otu_occupies[which(otu_occupies$type == 'hyper_rare'),]$occupies)
min(otu_occupies[which(otu_occupies$type == 'rare'),]$occupies)
max(otu_occupies[which(otu_occupies$type == 'rare'),]$occupies)
min(otu_occupies[which(otu_occupies$type == 'abundant'),]$occupies)
max(otu_occupies[which(otu_occupies$type == 'abundant'),]$occupies)
min(mean_r_abund_df[which(mean_r_abund_df$type == 'hyper_rare'),]$abund_perct)
max(mean_r_abund_df[which(mean_r_abund_df$type == 'hyper_rare'),]$abund_perct)
min(mean_r_abund_df[which(mean_r_abund_df$type == 'rare'),]$abund_perct)
max(mean_r_abund_df[which(mean_r_abund_df$type == 'rare'),]$abund_perct)
min(mean_r_abund_df[which(mean_r_abund_df$type == 'abundant'),]$abund_perct)
max(mean_r_abund_df[which(mean_r_abund_df$type == 'abundant'),]$abund_perct)
p1.2 <- ggplot(type_count_table, aes(x = Freq, y = Var1)) + 
    geom_bar(stat = 'identity', position = 'dodge', aes(fill = Var1)) +
    geom_text(aes(label = Freq), hjust = 1, size = 3) + 
    scale_fill_manual(values = args$taxa_type_part_color) + 
    theme(legend.position = 'none', axis.line = element_blank(), axis.text.x = element_blank(), 
          axis.title = element_blank(), axis.text.y = element_text(size = 10), axis.ticks.x = element_blank())
p1 <- p1.1 + annotation_custom(ggplotGrob(p1.2), xmin = 6, xmax = 31.86, ymin = -0.85, ymax = 2.4)   
p2.1 <- ggplot(mean_r_abund_df, aes(x = type, y = abund_perct)) + 
    geom_boxplot(aes(fill = type), outlier.shape = NA)  + 
    geom_hline(yintercept = 0.0001, linetype = 2, color = "gray50") + 
    geom_hline(yintercept = 0.001, linetype = 2, color = "gray50") + 
    scale_fill_manual(values = args$taxa_type_part_color) + 
    ylim(0, 0.004) + theme_gray() + xlab("ASV type") + ylab("Average relative abundance (%)") + 
    theme(legend.position = 'none', axis.title = element_text(size = 12), axis.text = element_text(size = 10))
p2.2 <- ggplot(mean_r_abund_df, aes(x = type, y = abund_perct)) + 
    geom_boxplot(aes(fill = type), outlier.size = 0.5) +
    xlab("ASV type") + ylab("Average relative abundance (%)") + 
    scale_fill_manual(values = args$taxa_type_part_color) + 
    theme(legend.position = 'none', axis.title = element_text(size = 10), axis.text = element_text(size = 8), axis.title.x = element_blank())
p2 <- p2.1 + annotation_custom(ggplotGrob(p2.2), xmin = 1.5, xmax = 3.635, ymin = 0.0015, ymax = 0.00425)
(p_asv_distri1 <- cowplot::plot_grid(p1, p2, nrow = 2, labels = c("A", "B"), label_size = 18, rel_heights = c(4, 3), align = 'hv'))
(p_asv_distri2 <- cowplot::plot_grid(p1, p2.2, nrow = 2, labels = c("A", "B"), label_size = 18, rel_heights = c(4, 3), align = 'hv'))
saved <- savePDF(object = p_asv_distri1, path = outdir_02, filename = "不同优势类群在样本中的分布情况-1.pdf", width = 5.5, height = 9.64)
saved <- savePDF(object = p_asv_distri2, path = outdir_02, filename = "不同优势类群在样本中的分布情况-2.pdf", width = 5.5, height = 9.64)

## 创建`microgeo` 数据集
datasets <- lapply(otu_tables, function(tab){
    dataset <- microtable$new(otu_table = t(tab) %>% as.data.frame, sample_table = gro_table, tax_table = tax_table)
    dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
    dataset$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]
    dataset$tax_table %<>% tidy_taxonomy()
    dataset$tidy_dataset()
    dataset$cal_abund()
    dataset$cal_alphadiv(PD = FALSE)
    dataset$cal_betadiv(unifrac = FALSE)
    dataset
})
names(datasets) <- names(otu_tables)
saveRDS(object = datasets, file = file.path(outdir_02, "microgeo-datasets.Rds"))

# 3. 环境参数分析
## 创建结果目录
outdir_03 <- "results/03.env-summary"
if (!outdir_03 %>% dir.exists) dir.create(outdir_03, recursive = T)

## 环境参数PCA分析、双因素PERMANOVA，两两PERMANOVA
sol_eucl_dist <- met_table %>% vegan::decostand(method = "standardize") %>% dist %>% as.matrix
if (length(unique(rownames(datasets$all$sample_table) == rownames(sol_eucl_dist))) > 1 ||
    !unique(rownames(datasets$all$sample_table) == rownames(sol_eucl_dist))) stop("sample id order error!")
if (length(unique(rownames(datasets$all$sample_table) == colnames(sol_eucl_dist))) > 1 ||
    !unique(rownames(datasets$all$sample_table) == colnames(sol_eucl_dist))) stop("sample id order error!")
datasets$all$beta_diversity$sol_eucl <- sol_eucl_dist
pca_sol <- trans_beta$new(dataset = datasets$all, measure = "sol_eucl", group = "Group")
pca_sol$cal_ordination(ordination = "PCA", trans_otu = T, scale_species = T)
(p_pca_sol <- pca_sol$plot_ordination(plot_color = "Elevation", plot_shape = "SlopeType", 
                                        plot_type = c("point"), point_size = 5, color_values = args$elev_color) + 
        theme(legend.position = c(0.5, 0.5)))
saved <- savePDF(object = p_pca_sol, path = outdir_03, filename = '土壤参数欧式距离PCA分析.pdf', 
                 width = 6.29, height = 5.90)
sol_permanova <- vegan::adonis(datasets$all$beta_diversity$sol_eucl %>% as.dist ~ Elevation + SlopeType, 
                               data = gro_table, permutations = 999) # Two-way PERMANOVA
sol_permanova_rst <- data.frame(Factor = c('Elevation', 'SlopeType'),
                                F.Model = c(sol_permanova$aov.tab$F.Model[1], sol_permanova$aov.tab$F.Model[2]),
                                R2 = c(sol_permanova$aov.tab$R2[1], sol_permanova$aov.tab$R2[2]),
                                P = c(sol_permanova$aov.tab$`Pr(>F)`[1], sol_permanova$aov.tab$`Pr(>F)`[2]))
pca_sol$cal_manova(group = "Group", manova_all = FALSE, p_adjust_method = 'fdr') # Paired PERMANOVA
saveMutipleXlsx(objectList = list(TwoWayPERMANOVA = sol_permanova_rst, PairedPERMANOVA = pca_sol$res_manova),
                path = outdir_03, filename = '土壤参数欧式距离PERMANOVA检验.xlsx', row.names = F, col.names = T)

## 土壤参数差异显著性检验之前的假设检验
### 准备数据：NS和CS
grop_ns <- gro_table[which(gro_table$SlopeType == 'NS'),]
grop_cs <- gro_table[which(gro_table$SlopeType == 'CS'),]
meta_ns_raw <- met_table[rownames(grop_ns),] 
meta_cs_raw <- met_table[rownames(grop_cs),]
meta_ns <- data.frame(elev = grop_ns$Elevation, meta_ns_raw) %>% reshape2::melt()
meta_cs <- data.frame(elev = grop_cs$Elevation, meta_cs_raw) %>% reshape2::melt() 
meta_ns$elev <- factor(meta_ns$elev, levels = args$elev_order)
meta_cs$elev <- factor(meta_cs$elev, levels = args$elev_order)
meta_ns$variable <- factor(meta_ns$variable, levels = args$met_order)
meta_cs$variable <- factor(meta_cs$variable, levels = args$met_order)

### 准备数据：2900, 3102和3194
grop_2900 <- gro_table[which(gro_table$Elevation == '2900'),]
grop_3102 <- gro_table[which(gro_table$Elevation == '3102'),]
grop_3194 <- gro_table[which(gro_table$Elevation == '3194'),]
meta_2900_raw <- met_table[rownames(grop_2900),] 
meta_3102_raw <- met_table[rownames(grop_3102),]
meta_3194_raw <- met_table[rownames(grop_3194),]
meta_2900 <- data.frame(slop = grop_2900$SlopeType, meta_2900_raw) %>% reshape2::melt()
meta_3102 <- data.frame(slop = grop_3102$SlopeType, meta_3102_raw) %>% reshape2::melt() 
meta_3194 <- data.frame(slop = grop_3194$SlopeType, meta_3194_raw) %>% reshape2::melt() 
meta_2900$slop <- factor(meta_2900$slop, levels = args$slop_order)
meta_3102$slop <- factor(meta_3102$slop, levels = args$slop_order)
meta_3194$slop <- factor(meta_3194$slop, levels = args$slop_order)
meta_2900$variable <- factor(meta_2900$variable, levels = args$met_order)
meta_3102$variable <- factor(meta_3102$variable, levels = args$met_order)
meta_3194$variable <- factor(meta_3194$variable, levels = args$met_order)

### 正态性检验
ns_shapiro.test <- lapply(X = meta_ns_raw %>% ncol %>% seq, function(x){
    test_rst <- shapiro.test(meta_ns_raw[,x])
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
cs_shapiro.test <- lapply(X = meta_cs_raw %>% ncol %>% seq, function(x){
    test_rst <- shapiro.test(meta_cs_raw[,x])
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e2900_shapiro.test <- lapply(X = meta_2900_raw %>% ncol %>% seq, function(x){
    test_rst <- shapiro.test(meta_2900_raw[,x])
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e3102_shapiro.test <- lapply(X = meta_3102_raw %>% ncol %>% seq, function(x){
    test_rst <- shapiro.test(meta_3102_raw[,x])
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e3194_shapiro.test <- lapply(X = meta_3194_raw %>% ncol %>% seq, function(x){
    test_rst <- shapiro.test(meta_3194_raw[,x])
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
rownames(ns_shapiro.test) <- colnames(meta_ns_raw)
rownames(cs_shapiro.test) <- colnames(meta_cs_raw)
rownames(e2900_shapiro.test) <- colnames(meta_2900_raw)
rownames(e3102_shapiro.test) <- colnames(meta_3102_raw)
rownames(e3194_shapiro.test) <- colnames(meta_3194_raw)
ns_shapiro.test <- data.frame(variable = rownames(ns_shapiro.test), ns_shapiro.test, slope = rep('NS', nrow(ns_shapiro.test)))
cs_shapiro.test <- data.frame(variable = rownames(cs_shapiro.test), cs_shapiro.test, slope = rep('CS', nrow(cs_shapiro.test)))
e2900_shapiro.test <- data.frame(variable = rownames(e2900_shapiro.test), e2900_shapiro.test, elev = rep('2900', nrow(e2900_shapiro.test)))
e3102_shapiro.test <- data.frame(variable = rownames(e3102_shapiro.test), e3102_shapiro.test, elev = rep('3102', nrow(e3102_shapiro.test)))
e3194_shapiro.test <- data.frame(variable = rownames(e3194_shapiro.test), e3194_shapiro.test, elev = rep('3194', nrow(e3194_shapiro.test)))
shapiro.test_list <- list(NS = ns_shapiro.test, CS = cs_shapiro.test,
                          E2900 = e2900_shapiro.test, E3102 = e3102_shapiro.test, E3194 = e3194_shapiro.test)
saveMutipleXlsx(objectList = shapiro.test_list, path = outdir_03, 
                filename = '土壤参数shapiro.test正态性检验.xlsx', row.names = F, col.names = T)

### 方差齐性检验
ns_bartlett.test <- lapply(X = meta_ns_raw %>% ncol %>% seq, function(x){
    test_rst <- bartlett.test(list(meta_ns_raw[grop_ns[which(grop_ns$Elevation == 2900),]$ID,x],
                                   meta_ns_raw[grop_ns[which(grop_ns$Elevation == 3102),]$ID,x],
                                   meta_ns_raw[grop_ns[which(grop_ns$Elevation == 3194),]$ID,x]))
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
cs_bartlett.test <- lapply(X = meta_cs_raw %>% ncol %>% seq, function(x){
    test_rst <- bartlett.test(list(meta_cs_raw[grop_cs[which(grop_cs$Elevation == 2900),]$ID,x],
                                   meta_cs_raw[grop_cs[which(grop_cs$Elevation == 3102),]$ID,x],
                                   meta_cs_raw[grop_cs[which(grop_cs$Elevation == 3194),]$ID,x]))
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e2900_bartlett.test <- lapply(X = meta_2900_raw %>% ncol %>% seq, function(x){
    test_rst <- bartlett.test(list(meta_2900_raw[grop_2900[which(grop_2900$SlopeType == 'NS'),]$ID,x],
                                   meta_2900_raw[grop_2900[which(grop_2900$SlopeType == 'CS'),]$ID,x]))
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e3102_bartlett.test <- lapply(X = meta_3102_raw %>% ncol %>% seq, function(x){
    test_rst <- bartlett.test(list(meta_3102_raw[grop_3102[which(grop_3102$SlopeType == 'NS'),]$ID,x],
                                   meta_3102_raw[grop_3102[which(grop_3102$SlopeType == 'CS'),]$ID,x]))
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e3194_bartlett.test <- lapply(X = meta_3194_raw %>% ncol %>% seq, function(x){
    test_rst <- bartlett.test(list(meta_3194_raw[grop_3194[which(grop_3194$SlopeType == 'NS'),]$ID,x],
                                   meta_3194_raw[grop_3194[which(grop_3194$SlopeType == 'CS'),]$ID,x]))
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
rownames(ns_bartlett.test) <- colnames(meta_ns_raw)
rownames(cs_bartlett.test) <- colnames(meta_cs_raw)
rownames(e2900_bartlett.test) <- colnames(meta_2900_raw)
rownames(e3102_bartlett.test) <- colnames(meta_3102_raw)
rownames(e3194_bartlett.test) <- colnames(meta_3194_raw)
ns_bartlett.test <- data.frame(variable = rownames(ns_bartlett.test), ns_bartlett.test, slope = rep('NS', nrow(ns_bartlett.test)))
cs_bartlett.test <- data.frame(variable = rownames(cs_bartlett.test), cs_bartlett.test, slope = rep('CS', nrow(cs_bartlett.test)))
e2900_bartlett.test <- data.frame(variable = rownames(e2900_bartlett.test), e2900_bartlett.test, elev = rep('2900', nrow(e2900_bartlett.test)))
e3102_bartlett.test <- data.frame(variable = rownames(e3102_bartlett.test), e3102_bartlett.test, elev = rep('3102', nrow(e3102_bartlett.test)))
e3194_bartlett.test <- data.frame(variable = rownames(e3194_bartlett.test), e3194_bartlett.test, elev = rep('3194', nrow(e3194_bartlett.test)))
bartlett.test_list <- list(NS = ns_bartlett.test, CS = cs_bartlett.test,
                          E2900 = e2900_bartlett.test, E3102 = e3102_bartlett.test, E3194 = e3194_bartlett.test)
saveMutipleXlsx(objectList = bartlett.test_list, path = outdir_03, 
                filename = '土壤参数bartlett.test方差齐性检验.xlsx', row.names = F, col.names = T)

### Kruskal-Wallis检验
ns_kruskal.test <- lapply(X = meta_ns_raw %>% ncol %>% seq, function(x){
    test_rst <- kruskal.test(meta_ns_raw[,x], grop_ns$Elevation)
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
cs_kruskal.test <- lapply(X = meta_cs_raw %>% ncol %>% seq, function(x){
    test_rst <- kruskal.test(meta_cs_raw[,x], grop_cs$Elevation)
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e2900_kruskal.test <- lapply(X = meta_2900_raw %>% ncol %>% seq, function(x){
    test_rst <- kruskal.test(meta_2900_raw[,x], grop_2900$SlopeType)
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e3102_kruskal.test <- lapply(X = meta_3102_raw %>% ncol %>% seq, function(x){
    test_rst <- kruskal.test(meta_3102_raw[,x], grop_3102$SlopeType)
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
e3194_kruskal.test <- lapply(X = meta_3194_raw %>% ncol %>% seq, function(x){
    test_rst <- kruskal.test(meta_3194_raw[,x], grop_3194$SlopeType)
    data.frame(statistic = test_rst$statistic, p = test_rst$p.value)
}) %>% do.call('rbind', .)
rownames(ns_kruskal.test) <- colnames(meta_ns_raw)
rownames(cs_kruskal.test) <- colnames(meta_cs_raw)
rownames(e2900_kruskal.test) <- colnames(meta_2900_raw)
rownames(e3102_kruskal.test) <- colnames(meta_3102_raw)
rownames(e3194_kruskal.test) <- colnames(meta_3194_raw)
ns_kruskal.test <- data.frame(variable = rownames(ns_kruskal.test), ns_kruskal.test, slope = rep('NS', nrow(ns_kruskal.test)))
cs_kruskal.test <- data.frame(variable = rownames(cs_kruskal.test), cs_kruskal.test, slope = rep('CS', nrow(cs_kruskal.test)))
e2900_kruskal.test <- data.frame(variable = rownames(e2900_kruskal.test), e2900_kruskal.test, elev = rep('2900', nrow(e2900_kruskal.test)))
e3102_kruskal.test <- data.frame(variable = rownames(e3102_kruskal.test), e3102_kruskal.test, elev = rep('3102', nrow(e3102_kruskal.test)))
e3194_kruskal.test <- data.frame(variable = rownames(e3194_kruskal.test), e3194_kruskal.test, elev = rep('3194', nrow(e3194_kruskal.test)))
kruskal.test_list <- list(NS = ns_kruskal.test, CS = cs_kruskal.test,
                           E2900 = e2900_kruskal.test, E3102 = e3102_kruskal.test, E3194 = e3194_kruskal.test)
saveMutipleXlsx(objectList = kruskal.test_list, path = outdir_03, 
                filename = '土壤参数kruskal.test组间差异检验.xlsx', row.names = F, col.names = T)

## 绘制土壤参数柱状图并进行两两比较
compare_list_elev <- list(c('2900', '3102'), c('3102', '3194'), c('2900', '3194'))
meta_2900$elev <- rep(2900, nrow(meta_2900))
meta_3102$elev <- rep(3102, nrow(meta_3102))
meta_3194$elev <- rep(3194, nrow(meta_3194))
meta_elev <- rbind(meta_2900, meta_3102, meta_3194)
(p.env.ns <- ggpubr::ggbarplot(meta_ns, x = "elev", y = "value", add = c("mean_se"), width = 0.7, 
                               ylab = "Observed values", fill = 'elev', position = position_dodge(0.7), 
                               facet.by = 'variable', scales = 'free_y', nrow = 2) +
        theme_gray() + 
        scale_fill_manual(values = args$elev_color) + 
        theme(legend.position = 'none', axis.title = element_text(size = 12), axis.text = element_text(size = 10),
              strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
        ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = compare_list_elev, 
                                   label = "p.signif", hide.ns = TRUE, size = 6, color = 'black', vjust = 0.5))
(p.env.cs <- ggpubr::ggbarplot(meta_cs, x = "elev", y = "value", add = c("mean_se"), width = 0.7, 
                               ylab = "Observed values", fill = 'elev', position = position_dodge(0.7), 
                               facet.by = 'variable', scales = 'free_y', nrow = 2) +
        theme_gray() + 
        scale_fill_manual(values = args$elev_color) + 
        theme(legend.position = 'none', axis.title = element_text(size = 12), axis.text = element_text(size = 10),
              strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
        ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = compare_list_elev, 
                                   label = "p.signif", hide.ns = TRUE, size = 6, color = 'black', vjust = 0.5))
(p.env.elev <- ggpubr::ggbarplot(meta_elev, x = "elev", y = "value", add = c("mean_se"), width = 0.7, xlab = "Elevation (m)", 
                                 ylab = "Observed values", fill = 'slop', position = position_dodge(0.7), 
                                 facet.by = 'variable', scales = 'free_y', nrow = 2) +
        theme_gray() + 
        scale_fill_manual(values = args$slop_color) + 
        theme(legend.position = c(0.9, 0.3), axis.title = element_text(size = 12), axis.text = element_text(size = 10),
              strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
        ggpubr::stat_compare_means(method = 'wilcox.test', aes(group = slop), 
                                   label = "p.signif", hide.ns = TRUE, size = 6, color = 'black', vjust = 0.5))
saved <- savePDF(object = p.env.ns, path = outdir_03, filename = "NS土壤参数在海拔上的变化.pdf", width = 11.04, height = 6.91)
saved <- savePDF(object = p.env.cs, path = outdir_03, filename = "CS土壤参数在海拔上的变化.pdf", width = 11.04, height = 6.91)
saved <- savePDF(object = p.env.elev, path = outdir_03, filename = "海拔土壤境参数在坡型间的变化.pdf", width = 11.04, height = 6.91)

# 4. 群物种组成
## 创建结果目录
outdir_04 <- "results/04.taxa_compos"
if (!outdir_04 %>% dir.exists) dir.create(outdir_04, recursive = T)

## 物种组成差异分析
abund_rst <- lapply(seq(from = 2, to = datasets %>% length), function(x){
    
    ## 提取要用的数据集
    dataset <- datasets[[x]]
    dataset_name <- names(datasets)[x]
    
    ## 差异分析: Lefse分析
    lefse_rst <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", taxa_level = "Genus", alpha = 0.01)
    p_lefse <- lefse_rst$plot_diff_bar(threshold = 4) + theme_gray() + 
        theme(legend.position = 'right', axis.title = element_text(size= 12), axis.text = element_text(size = 10), 
              legend.title = element_blank(), legend.text = element_text(size = 12))
    
    ## 提取差异的属，并将差异属下面的ASV编号找出来
    sig_genus <- lefse_rst$res_diff[which(lefse_rst$res_diff$P.adj <= 0.01 & lefse_rst$res_diff$LDA >= 4),]$Taxa
    genus2asv <- lapply(sig_genus, function(g){
        gg <- strsplit(g, split = '|', fixed = TRUE) %>% unlist; gg <- gg[gg %>% length]
        asv <- dataset$tax_table[which(dataset$tax_table$Genus == gg),] %>% rownames
        data.frame(genus = rep(gg, asv %>% length), asv = asv)
    }) %>% do.call('rbind', .)
    
    ## 相关分析
    sig_abund_dat <- dataset$taxa_abund$Genus[sig_genus,]
    genus_names <- lapply(rownames(sig_abund_dat), function(nam){
        part_nam_split <- strsplit(nam, 'g__', fixed = T) %>% unlist
        part_name <- part_nam_split[2]
    }) %>% unlist()
    rownames(sig_abund_dat) <- genus_names
    sig_abund_dat <- t(sig_abund_dat)
    sig_abund_dat <- sig_abund_dat[met_table %>% rownames, ]
    corr_rst <- psych::corr.test(sig_abund_dat, met_table, method = 'spearman', adjust = 'fdr')
    r_df <- reshape2::melt(corr_rst$r)
    p_df <- reshape2::melt(corr_rst$p.adj)
    colnames(r_df) <- c('genus', 'env', 'r')
    r_df$p <- p_df$value
    r_df$type <- rep(dataset_name, nrow(r_df))

    ## 准备差异属热图绘图数据
    sig_abund <- dataset$taxa_abund$Genus[sig_genus,]
    genus_names <- lapply(rownames(sig_abund), function(nam){
        full_nam_split <- strsplit(nam, '|', fixed = T) %>% unlist
        part_nam_split <- strsplit(nam, 'g__', fixed = T) %>% unlist
        full_name <- full_nam_split[length(full_nam_split)]
        part_name <- part_nam_split[2]
        asv_number <- genus2asv[which(genus2asv$genus == full_name),]$asv %>% length
        paste0(part_name, ' [', asv_number, ']')
    }) %>% unlist()
    sig_abund <- data.frame(genus = genus_names, sig_abund) %>% reshape2::melt()
    sig_abund$elev <- lapply(sig_abund$variable, function(v) {
        dataset$sample_table[which(dataset$sample_table$ID == v),]$Elevation
    }) %>% unlist
    sig_abund$slop <- lapply(sig_abund$variable, function(v) {
        dataset$sample_table[which(dataset$sample_table$ID == v),]$SlopeType
    }) %>% unlist
    sig_abund$group <- paste0(sig_abund$elev, '_', sig_abund$slop)
    sig_abund$type <- rep(dataset_name, nrow(sig_abund))
    
    ## 将差异的属提取出来，做功能预测，看看这些属主要负责哪些功能
    dataset_new <- microtable$new(otu_table = dataset$otu_table[genus2asv$asv,], 
                                  sample_table = dataset$sample_table, tax_table = dataset$tax_table) ## 重构数据集
    dataset_new$tax_table %<>% tidy_taxonomy()
    dataset_new$tidy_dataset()
    func_pred <- trans_func$new(dataset_new)
    func_pred$cal_spe_func(prok_database = "FAPROTAX") ## 功能预测
    func_pred_tab <- func_pred$res_spe_func
    func_pred_tab <- func_pred_tab[genus2asv$asv,]
    if (length(unique(rownames(func_pred_tab) == genus2asv$asv)) > 1 ||
        !unique(rownames(func_pred_tab) == genus2asv$asv)) stop('ASV order error!')
    func_pred_tab <- data.frame(taxa = genus2asv$genus, asv = genus2asv$asv, func_pred_tab)
    func_pred_tab_df <- func_pred_tab %>% reshape2::melt()
    func_pred_tab_df.agg <- aggregate(func_pred_tab_df$value, 
                                      by = list(func_pred_tab_df$taxa, func_pred_tab_df$variable), FUN = sum)
    colnames(func_pred_tab_df.agg) <- c('genus', 'func', 'freq')
    func_pred_tab_df.agg$genus <- lapply(func_pred_tab_df.agg$genus, function(genus){
        genus_split <- strsplit(genus, 'g__', fixed = T) %>% unlist(); genus_split[2]
    }) %>% unlist
    func_pred_tab_df.agg$type <- rep(dataset_name, nrow(func_pred_tab_df.agg))
    
    ## 返回数据
    list(lefse = p_lefse, cor = r_df, abund = sig_abund, func = func_pred_tab_df.agg)
})
names(abund_rst) <- names(datasets)[2:4]
(p_lefse <- cowplot::plot_grid(abund_rst$all$lefse, abund_rst$abundant$lefse, 
                               abund_rst$rare$lefse, ncol = 1, 
                               align = 'hv', rel_heights = c(8, 15, 5), labels = c('A', 'B', 'C')))
sig_taxa_cor <- lapply(abund_rst, function(obj) obj$cor) %>% do.call('rbind', .)
sig_taxa_cor$signif <- cut(sig_taxa_cor$p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", '')) 
sig_taxa_cor$coef <- ifelse(sig_taxa_cor$p <= 0.05, round(sig_taxa_cor$r, 2), '')
sig_taxa_cor$type <- factor(sig_taxa_cor$type, levels = args$taxa_type_full_order)
sig_taxa_cor$env <- factor(sig_taxa_cor$env, levels = args$met_order)
sig_abund_all <- lapply(abund_rst, function(obj) obj$abund) %>% do.call('rbind', .)
sig_abund_all$group <- factor(sig_abund_all$group, levels = c('2900_NS', '2900_CS', '3102_NS', '3102_CS', '3194_NS', '3194_CS'))
sig_abund_all$type <- factor(sig_abund_all$type, levels = args$taxa_type_full_order)
func_pred_tab_df.agg_all <- lapply(abund_rst, function(obj) obj$func) %>% do.call('rbind', .)
func_pred_tab_df.agg_all <- func_pred_tab_df.agg_all[which(func_pred_tab_df.agg_all$freq > 0),]
carbon_cycle_terms <- c("chitinolysis", "aerobic_chemoheterotrophy", "chemoheterotrophy",
                        "dark_hydrogen_oxidation", "photosynthetic_cyanobacteria",
                        "oxygenic_photoautotrophy", "photoautotrophy", "phototrophy")
nitrogen_cycle_terms <- c("aerobic_ammonia_oxidation", "nitrification", "nitrogen_fixation")
func_pred_tab_df.agg_all$ele_cycle <- ifelse(func_pred_tab_df.agg_all$func %in% carbon_cycle_terms, 'Carbon cycling', 'Nitrogen cycling')
func_pred_tab_df.agg_all$ele_cycle <- factor(func_pred_tab_df.agg_all$ele_cycle , levels = c('Carbon cycling', 'Nitrogen cycling'))
func_pred_tab_df.agg_all$type <- factor(func_pred_tab_df.agg_all$type, levels = args$taxa_type_full_order)
p_diff_genus_abund <- ggplot(sig_abund_all, aes(x = variable, y = genus, fill = log(value))) + geom_tile(color = 'black') + 
    scale_fill_gradientn(name = 'log(RB)', colors = RColorBrewer::brewer.pal(9, "Spectral") %>% rev, na.value = "white") + 
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_text(size = 10), 
          strip.text = element_text(size = 12), legend.position = 'right', 
          legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.ticks = element_blank()) + 
    facet_grid(type~ group, scales = 'free', space = 'free')
p_diff_genus_func <- ggplot(func_pred_tab_df.agg_all, aes(x = genus, y = func)) + 
    geom_tile(aes(fill = freq), color = 'black') + geom_text(aes(label = freq), size = 3) + ylab('Genus') + 
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Spectral") %>% rev, na.value = "white") + theme_gray() + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text = element_text(size = 10), axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          strip.text = element_text(size = 12), legend.title = element_blank(), legend.text = element_text(size = 12)) + 
    facet_grid(ele_cycle ~ type, scales = 'free', space = 'free')
p_diff_tax_cor_non_func <- ggplot(sig_taxa_cor[which(!sig_taxa_cor$genus %in% func_pred_tab_df.agg_all$genus),], 
                              aes(x = env, y = genus, fill = r)) + geom_tile(color = 'black') + 
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Spectral") %>% rev, na.value = "white") + 
    geom_text(aes(label = coef), vjust = 0.2) + geom_text(aes(label = signif), vjust = 1.5) + 
    theme(legend.position = 'none', axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12), axis.text = element_text(size = 10),
          strip.text = element_text(size = 12)) + facet_grid(type~., scales = 'free_y', space = 'free_y') 
p_diff_tax_cor_func <- ggplot(sig_taxa_cor[which(sig_taxa_cor$genus %in% func_pred_tab_df.agg_all$genus),], 
                              aes(x = env, y = genus, fill = r)) + geom_tile(color = 'black') + 
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Spectral") %>% rev, na.value = "white") + 
    geom_text(aes(label = coef), vjust = 0.2) + geom_text(aes(label = signif), vjust = 1.5) + 
    theme(legend.position = 'none', axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12), axis.text = element_text(size = 10),
          strip.text = element_text(size = 12)) + facet_grid(type~., scales = 'free_y', space = 'free_y') 
(p_comms <- cowplot::plot_grid(p_diff_genus_abund, p_diff_genus_func, nrow = 2, align = 'hv', labels = c('A', 'B'), rel_heights = c(5,3)))
saved <- savePDF(object = p_comms, path = outdir_04, filename = '差异属及其功能.pdf', width = 10.20, height = 10.20)
saved <- savePDF(object = p_lefse, path = outdir_04, filename = '属水平lefse分析结果.pdf', width = 10.20, height = 9.89)
saved <- savePDF(object = p_diff_tax_cor_non_func, path = outdir_04, filename = '未知功能差异属spearman相关.pdf', width = 6.18, height = 5.96)
saved <- savePDF(object = p_diff_tax_cor_func, path = outdir_04, filename = '已知功能差异属spearman相关.pdf', width = 8.35, height = 9.89)

# 5. alpha多样性分析
## 创建结果目录
outdir_05 <- "results/05.alpha-div"
if (!outdir_05 %>% dir.exists) dir.create(outdir_05, recursive = T)

## 整理alpha多样性绘图数据
alpha_rst_raw <- lapply(X = datasets %>% length %>% seq, function(x){
    
    ## 提取每一个数据集的Shannon指数
    dataset <- datasets[[x]]
    dataset_name <- names(datasets)[x]
    dat <- data.frame(row.names = rownames(dataset$alpha_diversity), 
                      sample_ids = rownames(dataset$alpha_diversity), 
                      alpha = dataset$alpha_diversity$Shannon)
    dat <- base::merge(x = dat, y = gro_table, by.x = 'sample_ids', by.y = 'ID', all.x = T, all.y = T)
    rownames(dat) <- dat$sample_ids
    dat$type <- rep(dataset_name, nrow(dat))
    
    ## 正态性检验
    ns_shapiro.test <- shapiro.test(dat[which(dat$SlopeType == 'NS'),]$alpha)
    cs_shapiro.test <- shapiro.test(dat[which(dat$SlopeType == 'CS'),]$alpha)
    e2900_shapiro.test <- shapiro.test(dat[which(dat$Elevation == '2900'),]$alpha)
    e3102_shapiro.test <- shapiro.test(dat[which(dat$Elevation == '3102'),]$alpha)
    e3194_shapiro.test <- shapiro.test(dat[which(dat$Elevation == '3194'),]$alpha)
    
    ## 方差齐性检验
    ns_bartlett.test <- bartlett.test(list(
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$alpha))
    cs_bartlett.test <- bartlett.test(list(
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$alpha))
    e2900_bartlett.test <- bartlett.test(list(
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$alpha))
    e3102_bartlett.test <- bartlett.test(list(
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$alpha))
    e3194_bartlett.test <- bartlett.test(list(
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$alpha,
        dat[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$alpha))
    
    ## Kruskal-Wallis检验  
    ns_kruskal.test <- kruskal.test(dat[which(dat$SlopeType == 'NS'),]$alpha, dat[which(dat$SlopeType == 'NS'),]$Elevation)
    cs_kruskal.test <- kruskal.test(dat[which(dat$SlopeType == 'CS'),]$alpha, dat[which(dat$SlopeType == 'CS'),]$Elevation)
    e2900_kruskal.test <- kruskal.test(dat[which(dat$Elevation == '2900'),]$alpha, dat[which(dat$Elevation == '2900'),]$SlopeType)
    e3102_kruskal.test <- kruskal.test(dat[which(dat$Elevation == '3102'),]$alpha, dat[which(dat$Elevation == '3102'),]$SlopeType)
    e3194_kruskal.test <- kruskal.test(dat[which(dat$Elevation == '3194'),]$alpha, dat[which(dat$Elevation == '3194'),]$SlopeType)
    
    ## 整合检验结果
    test_rst <- data.frame(type = rep(dataset_name, 5),
                           group = c('NS', 'CS', '2900', '3102', '3194'),
                           shapiro.test.statistic = c(ns_shapiro.test$statistic %>% as.numeric, cs_shapiro.test$statistic %>% as.numeric, 
                                                      e2900_shapiro.test$statistic %>% as.numeric, e3102_shapiro.test$statistic %>% as.numeric, e3194_shapiro.test$statistic %>% as.numeric),
                           shapiro.test.p.value = c(ns_shapiro.test$p.value %>% as.numeric, cs_shapiro.test$p.value %>% as.numeric, 
                                                    e2900_shapiro.test$p.value %>% as.numeric, e3102_shapiro.test$p.value %>% as.numeric, e3194_shapiro.test$p.value %>% as.numeric),
                           bartlett.test.statistic = c(ns_bartlett.test$statistic %>% as.numeric, cs_bartlett.test$statistic %>% as.numeric, 
                                                       e2900_bartlett.test$statistic %>% as.numeric, e3102_bartlett.test$statistic %>% as.numeric, e3194_bartlett.test$statistic %>% as.numeric),
                           bartlett.test.p.value = c(ns_bartlett.test$p.value %>% as.numeric, cs_bartlett.test$p.value %>% as.numeric, 
                                                     e2900_bartlett.test$p.value %>% as.numeric, e3102_bartlett.test$p.value %>% as.numeric, e3194_bartlett.test$p.value %>% as.numeric),
                           kruskal.test.statistic = c(ns_kruskal.test$statistic %>% as.numeric, cs_kruskal.test$statistic %>% as.numeric, 
                                                      e2900_kruskal.test$statistic %>% as.numeric, e3102_kruskal.test$statistic %>% as.numeric, e3194_kruskal.test$statistic %>% as.numeric),
                           kruskal.test.p.value = c(ns_kruskal.test$p.value %>% as.numeric, cs_kruskal.test$p.value %>% as.numeric, 
                                                    e2900_kruskal.test$p.value %>% as.numeric, e3102_kruskal.test$p.value %>% as.numeric, e3194_kruskal.test$p.value %>% as.numeric))
    list(dat = dat, tes = test_rst)
}) 
alpha_rst <- lapply(X = alpha_rst_raw %>% length %>% seq, function(x) alpha_rst_raw[[x]]$dat) %>% do.call('rbind', .)
test_alpha <- lapply(X = alpha_rst_raw %>% length %>% seq, function(x) alpha_rst_raw[[x]]$tes); names(test_alpha) <- names(datasets)
saveMutipleXlsx(objectList = test_alpha, path = outdir_05, filename = 'alpha多样性统计检验.xlsx', row.names = F, col.names = T)
alpha_rst$type <- factor(alpha_rst$type, levels = args$taxa_type_full_order)
alpha_rst$Elevation <- factor(alpha_rst$Elevation, levels = args$elev_order)
alpha_rst$SlopeType <- factor(alpha_rst$SlopeType, levels = args$slop_order)

## 绘制alpha多样性图
compare_list_elev <- list(c('2900', '3102'), c('3102', '3194'), c('2900', '3194'))
p.alpha.all <- ggplot(alpha_rst, aes(x = Elevation, y = alpha)) + geom_boxplot(aes(fill = SlopeType), outlier.size = 0.5) + 
    xlab('Elevation (m)') + ylab("Shannon-Weiner index") +  scale_fill_manual(values = args$slop_color) + 
    ggpubr::stat_compare_means(aes(group = SlopeType), label = "p.signif", method = "wilcox.test",
                               hide.ns = FALSE, size = 6, color = 'black', vjust = 1) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.2, 0.62), 
          legend.title = element_blank(), strip.text = element_text(size = 12)) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
p.alpha.ns <- ggplot(alpha_rst[which(alpha_rst$SlopeType == 'NS'),], aes(x = Elevation, y = alpha)) + 
    geom_boxplot(aes(fill = Elevation), outlier.size = 0.5) + 
    theme_bw() + ylab("Shannon-Weiner index") +  scale_fill_manual(values = args$elev_color) + 
    ggpubr::stat_compare_means(comparisons = compare_list_elev, label = "p.signif", method = "wilcox.test",
                               hide.ns = FALSE, size = 6, color = 'black', vjust = 0.7) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.2, 0.62), 
          legend.title = element_blank(), strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
p.alpha.cs <- ggplot(alpha_rst[which(alpha_rst$SlopeType == 'CS'),], aes(x = Elevation, y = alpha)) + 
    geom_boxplot(aes(fill = Elevation), outlier.size = 0.5) + 
    theme_bw() + ylab("Shannon-Weiner index") +  scale_fill_manual(values = args$elev_color) + 
    ggpubr::stat_compare_means(comparisons = compare_list_elev, label = "p.signif", method = "wilcox.test",
                               hide.ns = FALSE, size = 6, color = 'black', vjust = 0.7) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.2, 0.62), 
          legend.title = element_blank(), strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
(p.alpha.slope <- cowplot::plot_grid(p.alpha.ns, p.alpha.cs, align = 'hv', labels = c('A', 'B')))
saved <- savePDF(object = p.alpha.slope, path = outdir_05, filename = 'NS和CS中alpha多样性随海拔的变化趋势.pdf', width = 8.24, height = 9.89)

## 计算所有ASV多样性和不同ASV多样性的相关性
all_alpha_data <- alpha_rst[which(alpha_rst$type == 'all'),]
colnames(all_alpha_data)[2] <- "alpha_all"
all_alpha_data <- all_alpha_data[,-which(colnames(all_alpha_data) == 'type')]
split_dat_names <- args$taxa_type_part_order
corr_df_data <- lapply(split_dat_names, function(name){
    tmp_alp_dat <- alpha_rst[which(alpha_rst$type == name),]
    if (length(unique(all_alpha_data$sample_ids == tmp_alp_dat$sample_ids)) > 1 || 
        !unique(all_alpha_data$sample_ids == tmp_alp_dat$sample_ids)){
        stop('error!')
    }
    res <- data.frame(all_alpha_data, alpha_split = tmp_alp_dat$alpha, type = rep(name, nrow(tmp_alp_dat)))
}) %>% do.call('rbind', .)
corr_df_data$type <- factor(corr_df_data$type, levels = args$taxa_type_part_order)
corr_df_data$SlopeType <- factor(corr_df_data$SlopeType, levels = args$slop_order)

## 绘制相关性散点图并和alpha多样性变化趋势图合并
p.alpha.2 <- ggplot(corr_df_data, aes(x = alpha_all, y = alpha_split, group = SlopeType)) + 
    geom_point(aes(fill = SlopeType), color = 'gray70', shape = 21) + 
    geom_smooth(method = 'lm', aes(color = SlopeType)) + 
    scale_fill_manual(values = args$slop_color) + 
    scale_color_manual(values = args$slop_color) + 
    ggpubr::stat_cor(aes(color = SlopeType), method = 'pearson') + 
    xlab("Shannon-Weiner index of all taxa") + ylab("Shannon-Weiner index of ranked taxa") + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), 
          strip.text = element_text(size = 12), legend.position = c(0.8, 0.82)) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
(p.alpha <- cowplot::plot_grid(p.alpha.all, p.alpha.2, nrow = 1, labels = c("A", "B"), label_size = 18, align ='hv'))
saved <- savePDF(object = p.alpha, path = outdir_05, filename = "alpha多样性变化趋势图.pdf", width = 8.28, height = 11.27)

## 计算alpha多样性驱动因子: Shannon指数和ASV numbers
alpha_cor_rst <- lapply(X = datasets %>% length %>% seq, function(x){
    dataset <- datasets[[x]]
    dataset_name <- names(datasets)[x]
    grop_ns <- gro_table[which(gro_table$SlopeType == 'NS'),]
    grop_cs <- gro_table[which(gro_table$SlopeType == 'CS'),]
    meta_ns <- met_table[rownames(grop_ns),] 
    meta_cs <- met_table[rownames(grop_cs),]
    alpha_div_ns <- dataset$alpha_diversity[rownames(grop_ns),]
    alpha_div_cs <- dataset$alpha_diversity[rownames(grop_cs),]
    alpha_ns <- data.frame(row.names = rownames(alpha_div_ns), shannon = alpha_div_ns$Shannon, observed = alpha_div_ns$Observed)
    alpha_cs <- data.frame(row.names = rownames(alpha_div_cs), shannon = alpha_div_cs$Shannon, observed = alpha_div_cs$Observed)
    if (length(unique(rownames(meta_ns) == rownames(alpha_ns))) > 1 || !unique(rownames(meta_ns) == rownames(alpha_ns))){
        stop('error-01')
    }
    if (length(unique(rownames(meta_cs) == rownames(alpha_cs))) > 1 || !unique(rownames(meta_cs) == rownames(alpha_cs))){
        stop('error-02')
    }
    corr_rst_ns <- psych::corr.test(x = alpha_ns, y = meta_ns, method = 'spearman', adjust = 'fdr')
    corr_rst_cs <- psych::corr.test(x = alpha_cs, y = meta_cs, method = 'spearman', adjust = 'fdr')
    res <- data.frame(r = c(corr_rst_ns$r[1,], corr_rst_ns$r[2,], corr_rst_cs$r[1,], corr_rst_cs$r[2,]), 
                      p = c(corr_rst_ns$p[1,], corr_rst_ns$p[2,], corr_rst_cs$p[1,], corr_rst_cs$p[2,]),
                      var = c(corr_rst_ns$r %>% colnames, corr_rst_ns$r %>% colnames, corr_rst_cs$r %>% colnames, corr_rst_cs$r %>% colnames),
                      slope = c(rep("NS", (corr_rst_ns$r %>% ncol) * 2), rep("CS", (corr_rst_cs$r %>% ncol) * 2)),
                      index = c(rep('shannon', corr_rst_ns$r %>% ncol), rep('observed', corr_rst_ns$r %>% ncol),
                                rep('shannon', corr_rst_cs$r %>% ncol), rep('observed', corr_rst_cs$r %>% ncol)),
                      type = rep(dataset_name, (corr_rst_ns$r %>% ncol) * 4))  
}) %>% do.call('rbind', .)
alpha_cor_rst$signif <- cut(alpha_cor_rst$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
alpha_cor_rst$var <- factor(alpha_cor_rst$var, levels = args$met_order)
alpha_cor_rst$type <- factor(alpha_cor_rst$type, levels = args$taxa_type_full_order)
alpha_cor_rst$slope <- factor(alpha_cor_rst$slope, levels = args$slop_order)
(p_alpha_drivers_shannon <- ggplot(alpha_cor_rst[which(alpha_cor_rst$index == 'shannon'),], aes(x = var, y = type, fill = r)) + 
        geom_tile() + geom_text(aes(label = paste0(round(r, 2), "\n", signif)), vjust = 0.6, size = 4) + 
        scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Spectral"), na.value = "white") + theme_gray() + 
        theme(legend.position = 'right', axis.title = element_blank(), axis.text = element_text(size = 10),
              strip.text = element_text(size = 12), legend.title = element_text(size = 12), 
              legend.text = element_text(size = 12)) + 
        facet_wrap(~slope, nrow = 2))
(p_alpha_drivers_observed <- ggplot(alpha_cor_rst[which(alpha_cor_rst$index == 'observed'),], aes(x = var, y = type, fill = r)) + 
        geom_tile() + geom_text(aes(label = paste0(round(r, 2), "\n", signif)), vjust = 0.6, size = 4) + 
        scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Spectral"), na.value = "white") + theme_gray() + 
        theme(legend.position = 'right', axis.title = element_blank(), axis.text = element_text(size = 10),
              strip.text = element_text(size = 12), legend.title = element_text(size = 12), 
              legend.text = element_text(size = 12)) + 
        facet_wrap(~slope, nrow = 2))
saved <- savePDF(object = p_alpha_drivers_shannon, path = outdir_05, filename = "alpha多样性变化驱动因子-shannon.pdf", width = 11.45, height = 8.64)
saved <- savePDF(object = p_alpha_drivers_observed, path = outdir_05, filename = "alpha多样性变化驱动因子-observed.pdf", width = 11.45, height = 8.64)

# 6. beta多样性分析
## 创建结果目录
outdir_06 <- "results/06.beta-div"
if (!outdir_06 %>% dir.exists) dir.create(outdir_06, recursive = T)

## 加一个聚类分析
cluster_rst <- lapply(X = datasets %>% length %>% seq, function(x){
    dataset <- clone(datasets[[x]])
    t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = 'bray')
    t1$plot_clustering(group = "Group", measure = 'bray') + theme_gray() + 
        theme(legend.position = 'right', axis.text.y = element_blank(), 
              axis.title.y = element_blank(), axis.ticks.y = element_blank())
})
names(cluster_rst) <- names(datasets)
(p.beta.clustering <- cowplot::plot_grid(cluster_rst$all, cluster_rst$abundant, 
                                         cluster_rst$rare, cluster_rst$hyper_rare, 
                                         nrow = 2, labels = c("A", "B", "C", "D"), label_size = 18))
saved <- savePDF(object = p.beta.clustering, path = outdir_06, filename = '聚类分析.pdf', width = 8.35, height = 9.89)

## NMDS等分析
beta_data_rst <- lapply(X = datasets %>% length %>% seq, function(x){
    
    ## 提取数据集构建NMDS模型
    dataset <- datasets[[x]]
    t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
    t1$cal_ordination(ordination = "NMDS")
    p.nmds <- t1$plot_ordination(plot_color = "Elevation", plot_shape = "SlopeType", 
                       color_values = args$elev_color, plot_type = c("point", "ellipse"), 
                       point_size = 3) + theme(legend.position = c(0.2, 0.4))
    ## 绘制sheppard图
    nmds_coords <- t1$res_ordination$model$points # 提取 NMDS 的坐标数据
    nmds_coords_1d <- nmds_coords[, 1] # 计算每个样本的 NMDS 坐标在第一维的位置
    rank_percent <- ((1:length(nmds_coords_1d)) - 0.5) / length(nmds_coords_1d) # 计算标准化秩百分位数
    norm_percent <- pnorm(nmds_coords_1d) # 计算对应的正态分布累积概率
    df <- data.frame(rank_percent, norm_percent) # 创建数据框
    p.sp <- ggplot(df, aes(x = rank_percent, y = norm_percent)) + 
        geom_point(fill = 'white', color = 'black', shape = 21, size = 4) + 
        xlim(0, 1) + ylim(0, 1) + 
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + 
        theme_bw() + xlab("Standardized rank percentile") + ylab("Normal distribution cumulative probability") + 
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), 
              strip.text = element_text(size = 12), legend.position = c(0.8, 0.82)) + 
        coord_fixed()
    
    ## Two-way PERMANOVA
    result  <- vegan::adonis(dataset$beta_diversity$bray %>% as.dist ~ Elevation + SlopeType, data = dataset$sample_table)
    fv.elev <- result$aov.tab$F.Model[1]
    fv.slop <- result$aov.tab$F.Model[2]
    r2.elev <- result$aov.tab$R2[1]
    r2.slop <- result$aov.tab$R2[2]
    p.elev <- result$aov.tab$`Pr(>F)`[1]
    p.slop <- result$aov.tab$`Pr(>F)`[2]
    tw_adonis <- data.frame(type = rep(names(datasets)[x], 2), 
                            f = c(fv.elev, fv.slop), r2 = c(r2.elev, r2.slop), 
                            p = c(p.elev , p.slop), group = c("Elevation", "SlopeType"))
    
    ## Paired PERMANOVA
    t1$cal_manova(group = "Group", manova_all = FALSE, p_adjust_method = 'fdr') # Paired PERMANOVA
    res_manova <- t1$res_manova
    res_manova$type <- rep(names(datasets)[x], nrow(res_manova))
    
    ## RDA分析
    rda_res <- trans_env$new(dataset = dataset, add_data = met_table %>% vegan::decostand(method = "standardize"))
    rda_res$cal_ordination(method = "dbRDA", use_measure = "bray")
    rda_res$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
    rda_res$res_ordination_trans$df_sites$Elevation <- factor(rda_res$res_ordination_trans$df_sites$Elevation, levels = args$elev_order)
    rda_res$res_ordination_trans$df_sites$SlopeType <- factor(rda_res$res_ordination_trans$df_sites$SlopeType, levels = args$slop_order)
    p_rda <- rda_res$plot_ordination(plot_color = "Elevation", plot_shape = "SlopeType", plot_type = c( "point"),
                                     point_size = 5, color_values = args$elev_color) + theme_gray() +
            theme(legend.position = c(0.9, 0.23))
    rda_res.envfit <- vegan::envfit(rda_res$res_ordination, rda_res$data_env)
    rda_res.envfit_df <- data.frame(Var = names(rda_res.envfit$vectors$r), R2 = rda_res.envfit$vectors$r, P = rda_res.envfit$vectors$pvals)
    rda_res.envfit_df <- rda_res.envfit_df[order(rda_res.envfit_df$R2, decreasing = T),]
    rda_res.envfit_df$Var <- factor(rda_res.envfit_df$Var, levels = as.character(rda_res.envfit_df$Var))
    rda_res.envfit_df$Signif <- cut(rda_res.envfit_df$P, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", '')) 
    p_envfit <- ggplot(rda_res.envfit_df, aes(x = Var, y = R2)) +
            geom_bar(stat = 'identity', position = 'dodge', fill = 'white', color = 'black') +
            geom_text(aes(label = paste0(round(R2, 2), " ", Signif)), vjust = -0.5) +
            ylab("R2 values") +
            theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = 'none',
                  axis.title.x = element_blank())
    var_sel <- rda_res.envfit_df[which(rda_res.envfit_df$P <= 0.05),] %>% rownames
    rda_res_new <-  trans_env$new(dataset = dataset, add_data = met_table[,var_sel] %>% vegan::decostand(method = "standardize"))
    rda_res_new$cal_ordination(method = "dbRDA", use_measure = "bray")
    rda_res_new$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
    rda_res_new$res_ordination_trans$df_sites$Elevation <- factor(rda_res_new$res_ordination_trans$df_sites$Elevation, levels = args$elev_order)
    rda_res_new$res_ordination_trans$df_sites$SlopeType <- factor(rda_res_new$res_ordination_trans$df_sites$SlopeType, levels = args$slop_order)
    p_rda_new <- rda_res_new$plot_ordination(plot_color = "Elevation", plot_shape = "SlopeType", plot_type = c( "point"),
                                     point_size = 5, color_values = args$elev_color) + theme_gray() +
        theme(legend.position = c(0.9, 0.23))
    
    ## 返回数据
    res <- list(p.nmds = p.nmds, p.sp = p.sp, tw_adonis = tw_adonis, 
                res_manova = res_manova, p_rda = p_rda, p_envfit = p_envfit, p_rda_new = p_rda_new)
})
names(beta_data_rst) <- names(datasets)
p.beta.nmds <- cowplot::plot_grid(beta_data_rst$all$p.nmds, beta_data_rst$abundant$p.nmds, 
                                  beta_data_rst$rare$p.nmds, beta_data_rst$hyper_rare$p.nmds, 
                                  nrow = 2, labels = c("A", "B", "C", "D"), label_size = 18)
p.beta.sp <- cowplot::plot_grid(beta_data_rst$all$p.sp, beta_data_rst$abundant$p.sp, 
                                beta_data_rst$rare$p.sp, beta_data_rst$hyper_rare$p.sp, 
                                nrow = 2, labels = c("A", "B", "C", "D"), label_size = 18)
p.rda <- cowplot::plot_grid(beta_data_rst$all$p_rda, beta_data_rst$abundant$p_rda, 
                            beta_data_rst$rare$p_rda, beta_data_rst$hyper_rare$p_rda, 
                            nrow = 2, labels = c("A", "B", "C", "D"), label_size = 18)
p.envfit <- cowplot::plot_grid(beta_data_rst$all$p_envfit, beta_data_rst$abundant$p_envfit, 
                               beta_data_rst$rare$p_envfit, beta_data_rst$hyper_rare$p_envfit, 
                               nrow = 2, labels = c("A", "B", "C", "D"), label_size = 18)
p.rda.new <- cowplot::plot_grid(beta_data_rst$all$p_rda_new, beta_data_rst$abundant$p_rda_new, 
                                beta_data_rst$rare$p_rda_new, beta_data_rst$hyper_rare$p_rda_new, 
                                nrow = 2, labels = c("A", "B", "C", "D"), label_size = 18)
saved <- savePDF(object = p.beta.sp, path = outdir_06, filename = 'sheppard-plot.pdf', width = 10.82, height = 9.89)
saved <- savePDF(object = p.rda, path = outdir_06, filename = 'rda-raw.pdf', width = 10.75, height = 9.89)
saved <- savePDF(object = p.envfit, path = outdir_06, filename = 'rda-envfit.pdf', width = 10.75, height = 9.89)
saved <- savePDF(object = p.rda.new, path = outdir_06, filename = 'rda-sels.pdf', width = 10.75, height = 9.89)
tw_adonis_rst <- lapply(beta_data_rst, function(obj) obj$tw_adonis) %>% do.call('rbind', .)
res_manova_rst<- lapply(beta_data_rst, function(obj) obj$res_manova) %>% do.call('rbind', .)
saveMutipleXlsx(objectList = list(tw_adonis_rst = tw_adonis_rst, res_manova_rst = res_manova_rst), path = outdir_06, 
                filename = 'adonis.xlsx', row.names = F, col.names = T)

## 在NS和CS里面分别进行RDA分析
rda_ns_cs_rst <- lapply(X = datasets %>% length %>% seq, function(x){
    dataset <- datasets[[x]]
    sample_ns <- dataset$sample_table[which(dataset$sample_table$SlopeType == 'NS'),]$ID
    sample_cs <- dataset$sample_table[which(dataset$sample_table$SlopeType == 'CS'),]$ID
    dataset_ns <- microtable$new(otu_table = dataset$otu_table[,sample_ns], 
                                 sample_table = dataset$sample_table[sample_ns,], 
                                 tax_table = dataset$tax_table)
    dataset_ns$tax_table %<>% tidy_taxonomy()
    dataset_ns$tidy_dataset()
    dataset_ns$cal_betadiv(unifrac = FALSE)
    dataset_cs <- microtable$new(otu_table = dataset$otu_table[,sample_cs], 
                                 sample_table = dataset$sample_table[sample_cs,], 
                                 tax_table = dataset$tax_table)
    dataset_cs$tax_table %<>% tidy_taxonomy()
    dataset_cs$tidy_dataset()
    dataset_cs$cal_betadiv(unifrac = FALSE)
    rda_res_ns <- trans_env$new(dataset = dataset_ns, add_data = met_table[sample_ns,] %>% vegan::decostand(method = "standardize"))  
    rda_res_cs <- trans_env$new(dataset = dataset_cs, add_data = met_table[sample_cs,] %>% vegan::decostand(method = "standardize")) 
    rda_res_ns$cal_ordination(method = "dbRDA", use_measure = "bray")
    rda_res_cs$cal_ordination(method = "dbRDA", use_measure = "bray")
    rda_res_ns$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
    rda_res_cs$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
    rda_res_ns$res_ordination_trans$df_sites$Elevation <- factor(rda_res_ns$res_ordination_trans$df_sites$Elevation, levels = args$elev_order)
    rda_res_cs$res_ordination_trans$df_sites$Elevation <- factor(rda_res_cs$res_ordination_trans$df_sites$Elevation, levels = args$elev_order)
    p_rda_ns <- rda_res_ns$plot_ordination(plot_color = "Elevation", plot_type = c( "point", "chull"),
                                     point_size = 2, color_values = args$elev_color) + theme_gray() +
        theme(legend.position = c(0.9, 0.23))
    p_rda_cs <- rda_res_cs$plot_ordination(plot_color = "Elevation", plot_type = c( "point", "chull"),
                                           point_size = 2, color_values = args$elev_color) + theme_gray() +
        theme(legend.position = c(0.9, 0.23))
    rda_res.envfit_ns <- vegan::envfit(rda_res_ns$res_ordination, rda_res_ns$data_env)
    rda_res.envfit_cs <- vegan::envfit(rda_res_cs$res_ordination, rda_res_cs$data_env)
    rda_res.envfit_df_ns <- data.frame(Var = names(rda_res.envfit_ns$vectors$r), R2 = rda_res.envfit_ns$vectors$r, P = rda_res.envfit_ns$vectors$pvals)
    rda_res.envfit_df_cs <- data.frame(Var = names(rda_res.envfit_cs$vectors$r), R2 = rda_res.envfit_cs$vectors$r, P = rda_res.envfit_cs$vectors$pvals)
    rda_res.envfit_df_ns <- rda_res.envfit_df_ns[order(rda_res.envfit_df_ns$R2, decreasing = T),]
    rda_res.envfit_df_cs <- rda_res.envfit_df_cs[order(rda_res.envfit_df_cs$R2, decreasing = T),]
    rda_res.envfit_df_ns$Var <- factor(rda_res.envfit_df_ns$Var, levels = as.character(rda_res.envfit_df_ns$Var))
    rda_res.envfit_df_cs$Var <- factor(rda_res.envfit_df_cs$Var, levels = as.character(rda_res.envfit_df_cs$Var))
    rda_res.envfit_df_ns$Signif <- cut(rda_res.envfit_df_ns$P, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", '')) 
    rda_res.envfit_df_cs$Signif <- cut(rda_res.envfit_df_cs$P, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", '')) 
    p_envfit_ns <- ggplot(rda_res.envfit_df_ns, aes(x = Var, y = R2)) +
        geom_bar(stat = 'identity', position = 'dodge', fill = 'white', color = 'black') +
        geom_text(aes(label = paste0(round(R2, 2), " ", Signif)), vjust = -0.5) +
        ylab("R2 values") +
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = 'none',
              axis.title.x = element_blank())
    p_envfit_cs <- ggplot(rda_res.envfit_df_cs, aes(x = Var, y = R2)) +
        geom_bar(stat = 'identity', position = 'dodge', fill = 'white', color = 'black') +
        geom_text(aes(label = paste0(round(R2, 2), " ", Signif)), vjust = -0.5) +
        ylab("R2 values") +
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = 'none',
              axis.title.x = element_blank())
    var_sel_ns <- rda_res.envfit_df_ns[which(rda_res.envfit_df_ns$P <= 0.05),] %>% rownames
    var_sel_cs <- rda_res.envfit_df_cs[which(rda_res.envfit_df_cs$P <= 0.05),] %>% rownames
    rda_res_new_ns <-  trans_env$new(dataset = dataset_ns, add_data = met_table[sample_ns,var_sel_ns] %>% vegan::decostand(method = "standardize"))
    rda_res_new_ns$cal_ordination(method = "dbRDA", use_measure = "bray")
    rda_res_new_ns$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
    rda_res_new_ns$res_ordination_trans$df_sites$Elevation <- factor(rda_res_new_ns$res_ordination_trans$df_sites$Elevation, levels = args$elev_order)
    rda_res_new_cs <-  trans_env$new(dataset = dataset_cs, add_data = met_table[sample_cs,var_sel_cs] %>% vegan::decostand(method = "standardize"))
    rda_res_new_cs$cal_ordination(method = "dbRDA", use_measure = "bray")
    rda_res_new_cs$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
    rda_res_new_cs$res_ordination_trans$df_sites$Elevation <- factor(rda_res_new_cs$res_ordination_trans$df_sites$Elevation, levels = args$elev_order)
    p_rda_new_ns <- rda_res_new_ns$plot_ordination(plot_color = "Elevation", plot_type = c( "point", "chull"),
                                             point_size = 2, color_values = args$elev_color) + theme_gray() +
        theme(legend.position = c(0.9, 0.23))
    p_rda_new_cs <- rda_res_new_cs$plot_ordination(plot_color = "Elevation", plot_type = c( "point", "chull"),
                                                   point_size = 2, color_values = args$elev_color) + theme_gray() +
        theme(legend.position = c(0.9, 0.23))
    res <- list(raw_rda = list(ns = p_rda_ns, cs = p_rda_cs), envfit = list(ns = p_envfit_ns, cs = p_envfit_cs),
                sel_rda = list(ns = p_rda_new_ns, cs = p_rda_new_cs))
})
names(rda_ns_cs_rst) <- names(datasets)
p.rda_ns_cs_raw <- cowplot::plot_grid(
    rda_ns_cs_rst$all$raw_rda$ns, rda_ns_cs_rst$abundant$raw_rda$ns, 
    rda_ns_cs_rst$rare$raw_rda$ns, rda_ns_cs_rst$hyper_rare$raw_rda$ns,
    rda_ns_cs_rst$all$raw_rda$cs, rda_ns_cs_rst$abundant$raw_rda$cs, 
    rda_ns_cs_rst$rare$raw_rda$cs, rda_ns_cs_rst$hyper_rare$raw_rda$cs, ncol = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"), label_size = 18
)
p.envfit_ns_cs_raw <- cowplot::plot_grid(
    rda_ns_cs_rst$all$envfit$ns, rda_ns_cs_rst$abundant$envfit$ns, 
    rda_ns_cs_rst$rare$envfit$ns, rda_ns_cs_rst$hyper_rare$envfit$ns,
    rda_ns_cs_rst$all$envfit$cs, rda_ns_cs_rst$abundant$envfit$cs, 
    rda_ns_cs_rst$rare$envfit$cs, rda_ns_cs_rst$hyper_rare$envfit$cs, ncol = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H"), label_size = 18
)
p.rda_ns_cs_new <- cowplot::plot_grid(
    rda_ns_cs_rst$all$sel_rda$ns, rda_ns_cs_rst$abundant$sel_rda$ns, 
    rda_ns_cs_rst$rare$sel_rda$ns, rda_ns_cs_rst$hyper_rare$sel_rda$ns,
    rda_ns_cs_rst$all$sel_rda$cs, rda_ns_cs_rst$abundant$sel_rda$cs, 
    rda_ns_cs_rst$rare$sel_rda$cs, rda_ns_cs_rst$hyper_rare$sel_rda$cs, ncol = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"), label_size = 18
)
saved <- savePDF(object = p.rda_ns_cs_raw, path = outdir_06, 
                 filename = 'NS和CS的RDA-原始模型.pdf', width = 13.14, height = 9.64)
saved <- savePDF(object = p.envfit_ns_cs_raw, path = outdir_06, 
                 filename = 'NS和CS的envfit.pdf', width = 13.14, height = 9.64)
saved <- savePDF(object = p.rda_ns_cs_new, path = outdir_06, 
                 filename = 'NS和CS的RDA-显著变量.pdf', width = 13.14, height = 9.64)

## 计算距离矩阵之间的关系
dist_rst <- lapply(X = datasets %>% length %>% seq, function(x){
    dataset <- datasets[[x]]
    dataset_name <- names(datasets)[x]
    bray <- dataset$beta_diversity$bray 
    ns_samples <- gro_table[which(gro_table$SlopeType == 'NS'),]$ID
    cs_samples <- gro_table[which(gro_table$SlopeType == 'CS'),]$ID
    ns_bray <- bray[ns_samples,ns_samples] %>% as.dist %>% as.vector
    cs_bray <- bray[cs_samples,cs_samples] %>% as.dist %>% as.vector   
    res <- data.frame(bray = c(ns_bray, cs_bray), 
                      slope = c(rep("NS", length(ns_bray)), rep("CS", length(cs_bray))),
                      type = rep(dataset_name, length(ns_bray) + length(cs_bray)))
}) %>%do.call('rbind', .)
dist_all_ns <- dist_rst[which(dist_rst$type == 'all' & dist_rst$slope == 'NS'),]
dist_ran_ns <- dist_rst[which(dist_rst$type != 'all' & dist_rst$slope == 'NS'),]
dist_all_cs <- dist_rst[which(dist_rst$type == 'all' & dist_rst$slope == 'CS'),]
dist_ran_cs <- dist_rst[which(dist_rst$type != 'all' & dist_rst$slope == 'CS'),]
dist_plotdat_ns <- data.frame(all_bray = rep(dist_all_ns$bray, 3), 
                              ran_bray = dist_ran_ns$bray, 
                              ran_type = dist_ran_ns$type,
                              slope = rep("NS", length(dist_ran_ns$bray)))
dist_plotdat_cs <- data.frame(all_bray = rep(dist_all_cs$bray, 3), 
                              ran_bray = dist_ran_cs$bray, 
                              ran_type = dist_ran_cs$type,
                              slope = rep("CS", length(dist_ran_cs$bray)))
dist_plotdat <- rbind(dist_plotdat_ns, dist_plotdat_cs)
dist_plotdat$ran_type <- factor(dist_plotdat$ran_type, levels = args$taxa_type_part_order)
dist_plotdat$slope <- factor(dist_plotdat$slope, levels = args$slop_order)
p.bray.cor <- ggplot(dist_plotdat, aes(x = all_bray, y = ran_bray, group = slope)) + 
    geom_point(aes(fill = slope), color = 'gray70', shape = 21) + 
    geom_smooth(method = 'lm', aes(color = slope)) + 
    scale_fill_manual(values = args$slop_color) + 
    scale_color_manual(values = args$slop_color) + 
    ggpubr::stat_cor(aes(color = slope), method = 'pearson') + 
    theme_bw() + xlab("Bray-Curtis distance of all taxa") + ylab("Bray-Curtis distance of ranked taxa") + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.8, 0.82)) + 
    facet_wrap(~ran_type, ncol = 1, scales = 'free_y')
saved <- savePDF(object = p.bray.cor, path = outdir_06, filename = "不同类群和整体beta多样性的相关分析.pdf", width = 5.84, height = 9.89)
(p.beta <- ggpubr::ggarrange(p.beta.nmds, p.bray.cor, ncol = 2, labels = c("", "E"), widths = c(5, 2),
                             font.label = list(size = 18, color = "black", face = "bold", family = NULL)))
saved <- savePDF(object = p.beta, path = outdir_06, filename = "beta多样性变化趋势图.pdf", width = 11.98, height = 8.17)

# 7. 功能冗余分析
## 创建结果目录
outdir_07 <- "results/07.func_red"
if (!outdir_07 %>% dir.exists) dir.create(outdir_07, recursive = T)

## 计算功能冗余
fun_red_res_file <- file.path(outdir_07, 'fun-red-results.Rds')
if (!file.exists(fun_red_res_file)){
    fr_rst <- lapply(X = datasets %>% length %>% seq, function(x){
        dataset <- datasets[[x]]
        otu_table_use <- dataset$otu_table
        kegg_copies_use <- kegg_copies[rownames(otu_table_use),]
        kegg_copies_use <- kegg_copies_use[,apply(kegg_copies_use, 2, sum) > 0]
        func_dist <- parallelDistCal(kegg_copies_use, threads = 90)
        fr <- FunctionalRedParallelDist(comm = otu_table_use, dis = func_dist, threads = 90)
    })
    names(fr_rst) <- names(datasets)
    saveRDS(object = fr_rst, file = fun_red_res_file)
}else{
    fr_rst <- readRDS(fun_red_res_file)
}

## 提取功能冗余数据用于绘图
fred_raw <- lapply(fr_rst %>% length %>% seq, function(x){
    
    ## 提取数据
    dataset <- datasets[[names(fr_rst)[x]]]
    object <- fr_rst[[x]]
    name <- names(fr_rst)[x]
    red <- data.frame(sample = rownames(object$red), R = 1 - object$red$U, type = rep(name, nrow(object$red)))
    red <- base::merge(x = red, y = dataset$sample_table, by.x = 'sample', by.y = 'ID', all.x = T, all.y = T)
    colnames(red) <- c('ID', 'r', 'type', 'Elevation', 'SlopeType', 'group')
    rownames(red) <- red$ID
    
    ## 正态性检验
    ns_shapiro.test <- shapiro.test(red[which(red$SlopeType == 'NS'),]$r)
    cs_shapiro.test <- shapiro.test(red[which(red$SlopeType == 'CS'),]$r)
    e2900_shapiro.test <- shapiro.test(red[which(red$Elevation == '2900'),]$r)
    e3102_shapiro.test <- shapiro.test(red[which(red$Elevation == '3102'),]$r)
    e3194_shapiro.test <- shapiro.test(red[which(red$Elevation == '3194'),]$r)
    
    ## 方差齐性检验
    ns_bartlett.test <- bartlett.test(list(
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$r))
    cs_bartlett.test <- bartlett.test(list(
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$r))
    e2900_bartlett.test <- bartlett.test(list(
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 2900 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$r))
    e3102_bartlett.test <- bartlett.test(list(
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3102 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$r))
    e3194_bartlett.test <- bartlett.test(list(
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'NS'),]$ID,]$r,
        red[dataset$sample_table[which(dataset$sample_table$Elevation == 3194 & dataset$sample_table$SlopeType == 'CS'),]$ID,]$r))
    
    ## Kruskal-Wallis检验  
    ns_kruskal.test <- kruskal.test(red[which(red$SlopeType == 'NS'),]$r, red[which(red$SlopeType == 'NS'),]$Elevation)
    cs_kruskal.test <- kruskal.test(red[which(red$SlopeType == 'CS'),]$r, red[which(red$SlopeType == 'CS'),]$Elevation)
    e2900_kruskal.test <- kruskal.test(red[which(red$Elevation == '2900'),]$r, red[which(red$Elevation == '2900'),]$SlopeType)
    e3102_kruskal.test <- kruskal.test(red[which(red$Elevation == '3102'),]$r, red[which(red$Elevation == '3102'),]$SlopeType)
    e3194_kruskal.test <- kruskal.test(red[which(red$Elevation == '3194'),]$r, red[which(red$Elevation == '3194'),]$SlopeType)
    
    ## 整合检验结果
    test_rst <- data.frame(type = rep(name, 5),
                           group = c('NS', 'CS', '2900', '3102', '3194'),
                           shapiro.test.statistic = c(ns_shapiro.test$statistic %>% as.numeric, cs_shapiro.test$statistic %>% as.numeric, 
                                                      e2900_shapiro.test$statistic %>% as.numeric, e3102_shapiro.test$statistic %>% as.numeric, e3194_shapiro.test$statistic %>% as.numeric),
                           shapiro.test.p.value = c(ns_shapiro.test$p.value %>% as.numeric, cs_shapiro.test$p.value %>% as.numeric, 
                                                    e2900_shapiro.test$p.value %>% as.numeric, e3102_shapiro.test$p.value %>% as.numeric, e3194_shapiro.test$p.value %>% as.numeric),
                           bartlett.test.statistic = c(ns_bartlett.test$statistic %>% as.numeric, cs_bartlett.test$statistic %>% as.numeric, 
                                                       e2900_bartlett.test$statistic %>% as.numeric, e3102_bartlett.test$statistic %>% as.numeric, e3194_bartlett.test$statistic %>% as.numeric),
                           bartlett.test.p.value = c(ns_bartlett.test$p.value %>% as.numeric, cs_bartlett.test$p.value %>% as.numeric, 
                                                     e2900_bartlett.test$p.value %>% as.numeric, e3102_bartlett.test$p.value %>% as.numeric, e3194_bartlett.test$p.value %>% as.numeric),
                           kruskal.test.statistic = c(ns_kruskal.test$statistic %>% as.numeric, cs_kruskal.test$statistic %>% as.numeric, 
                                                      e2900_kruskal.test$statistic %>% as.numeric, e3102_kruskal.test$statistic %>% as.numeric, e3194_kruskal.test$statistic %>% as.numeric),
                           kruskal.test.p.value = c(ns_kruskal.test$p.value %>% as.numeric, cs_kruskal.test$p.value %>% as.numeric, 
                                                    e2900_kruskal.test$p.value %>% as.numeric, e3102_kruskal.test$p.value %>% as.numeric, e3194_kruskal.test$p.value %>% as.numeric))
    list(red = red, test_rst = test_rst)
})
fred <- lapply(X = fred_raw %>% length %>% seq, function(x) fred_raw[[x]]$red) %>% do.call('rbind', .)
test_fred <- lapply(X = fred_raw %>% length %>% seq, function(x) fred_raw[[x]]$test_rst); names(test_fred) <- names(fr_rst)
saveMutipleXlsx(objectList = test_fred, path = outdir_07, filename = 'fred多样性统计检验.xlsx', row.names = F, col.names = T)
fred$type <- factor(fred$type, levels = args$taxa_type_full_order)
fred$Elevation <- factor(fred$Elevation, levels = args$elev_order)
fred$SlopeType <- factor(fred$SlopeType, levels = args$slop_order)

## 绘制功能冗余图变化趋势图
compare_list_elev <- list(c('2900', '3102'), c('3102', '3194'), c('2900', '3194'))
p.red.all <- ggplot(fred, aes(x = Elevation, y = r)) + 
    geom_boxplot(aes(fill = SlopeType), outlier.size = 0.5) + 
    xlab("Elevation (m)") + ylab("Functional redundancy") + 
    scale_fill_manual(values = args$slop_color) + 
    ggpubr::stat_compare_means(aes(group = SlopeType), label = "p.signif", method = "wilcox.test", vjust = 1, size = 6, hide.ns = T) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.2, 0.62),
          strip.text = element_text(size = 12)) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
p.red.ns <- ggplot(fred[which(fred$SlopeType == 'NS'),], aes(x = Elevation, y = r)) + 
    geom_boxplot(aes(fill = Elevation), outlier.size = 0.5) + 
    theme_bw() +xlab("Elevation (m)") + ylab("Functional redundancy") +   scale_fill_manual(values = args$elev_color) + 
    ggpubr::stat_compare_means(comparisons = compare_list_elev, label = "p.signif", method = "wilcox.test",
                               hide.ns = FALSE, size = 6, color = 'black', vjust = 1) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.2, 0.62), 
          legend.title = element_blank(), strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
p.red.cs <- ggplot(fred[which(fred$SlopeType == 'CS'),], aes(x = Elevation, y = r)) + 
    geom_boxplot(aes(fill = Elevation), outlier.size = 0.5) + 
    theme_bw() + xlab("Elevation (m)") + ylab("Functional redundancy") +  scale_fill_manual(values = args$elev_color) + 
    ggpubr::stat_compare_means(comparisons = compare_list_elev, label = "p.signif", method = "wilcox.test",
                               hide.ns = FALSE, size = 6, color = 'black', vjust = 1) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = c(0.2, 0.62), 
          legend.title = element_blank(), strip.text = element_text(size = 12), axis.title.x = element_blank()) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
(p.red.slope <- cowplot::plot_grid(p.red.ns, p.red.cs, align = 'hv', labels = c('A', 'B')))
saved <- savePDF(object = p.red.slope, path = outdir_07, filename = 'NS和CS中FR多样性随海拔的变化趋势.pdf', width = 8.24, height = 9.89)

## 计算所有ASV功能冗余和不同ASV功能冗余的相关性
all_red_data <- fred[which(fred$type == 'all'),]
colnames(all_red_data)[2] <- "r_all"
all_red_data <- all_red_data[,-which(colnames(all_red_data) == 'type')]
split_dat_names <- c('hyper_rare', "rare","abundant")
corr_df_data_red <- lapply(split_dat_names, function(name){
    tmp_red_dat <- fred[which(fred$type == name),]
    if (length(unique(all_red_data$ID == tmp_red_dat$ID)) > 1 || !unique(all_red_data$ID == tmp_red_dat$ID)){
        stop('error!')
    }
    res <- data.frame(all_red_data, r_split = tmp_red_dat$r, type = rep(name, nrow(tmp_red_dat)))
}) %>% do.call('rbind', .)
corr_df_data_red$type <- factor(corr_df_data_red$type, levels = args$taxa_type_part_order)
corr_df_data_red$SlopeType <- factor(corr_df_data_red$SlopeType, levels = args$slop_order)
p.red.2 <- ggplot(corr_df_data_red, aes(x = r_all, y = r_split, group = SlopeType)) + 
    geom_point(aes(fill = SlopeType), color = 'gray70', shape = 21) + 
    geom_smooth(method = 'lm', aes(color = SlopeType)) + 
    scale_fill_manual(values = args$slop_color) + 
    scale_color_manual(values = args$slop_color) + 
    ggpubr::stat_cor(aes(color = SlopeType), method = 'pearson') + 
    xlab("Functional redundancy of all taxa") + ylab("Functional redundancy of ranked taxa") + 
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), 
          strip.text = element_text(size = 12), legend.position = c(0.8, 0.82)) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
saved <- savePDF(object = p.red.2, path = outdir_07, filename = "不同类群和整体red的相关分析.pdf", width = 5.84, height = 9.89)
(p.red <- cowplot::plot_grid(p.red.all, p.red.2, nrow = 1, labels = c("A", "B"), label_size = 18, rel_widths = c(3,5)))
saved <- savePDF(object = p.red, path = outdir_07, filename = "功能冗余变化趋势图.pdf", width = 10.35, height = 9.64)

## 功能冗余驱动因素
fred_cor_rst <- lapply(X = unique(fred$type), function(type){
    tmp <- fred[which(fred$type == type),]
    tmp.ns <- tmp[which(tmp$SlopeType == 'NS'),]
    tmp.cs <- tmp[which(tmp$SlopeType == 'CS'),]
    ns_meta <- met_table[tmp.ns$ID,]
    cs_meta <- met_table[tmp.cs$ID,]
    corr_rst_ns <- psych::corr.test(x = tmp.ns$r, y = ns_meta, method = 'spearman', adjust = 'fdr')
    corr_rst_cs <- psych::corr.test(x = tmp.cs$r, y = cs_meta, method = 'spearman', adjust = 'fdr')
    res <- data.frame(r = c(corr_rst_ns$r, corr_rst_cs$r), 
                      p = c(corr_rst_ns$p, corr_rst_cs$p),
                      var = rep(corr_rst_ns$r %>% colnames, 2),
                      slope = c(rep("NS", corr_rst_ns$r %>% length), rep("CS", corr_rst_cs$r %>% length)),
                      type = rep(type, c(corr_rst_ns$r, corr_rst_cs$r) %>% length))    
}) %>% do.call('rbind', .)
fred_cor_rst$signif <- cut(fred_cor_rst$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
fred_cor_rst$var <- factor(fred_cor_rst$var, levels = args$met_order)
fred_cor_rst$type <- factor(fred_cor_rst$type, levels = args$taxa_type_full_order)
fred_cor_rst$slope <- factor(fred_cor_rst$slope, levels = args$slop_order)
(p_fred_drivers <- ggplot(fred_cor_rst, aes(x = var, y = type, fill = r)) + 
        geom_tile() + geom_text(aes(label = paste0(round(r, 2), "\n", signif)), vjust = 0.6, size = 3, fontface = 'bold') + 
        scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Spectral"), na.value = "white") + theme_gray() + 
        theme(legend.position = 'right', axis.title = element_blank(), axis.text = element_text(size = 8),
              strip.text = element_text(size = 10), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
        facet_wrap(~slope, nrow = 2))
saved <- savePDF(object = p_fred_drivers, path = outdir_07, filename = "功能冗余变化驱动因子.pdf", width = 8.74, height = 6.60)