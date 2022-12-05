
library(microbiomeMarker)
library(phyloseq)
library(pheatmap)
library(dplyr)
library(magrittr)

rm(list = ls())
setwd('~/Daily/jindundun/20221201-LEfSe')

#raw_grp <- grp <- c('CS', 'PS')
raw_grp <- c('Control', 'Treatment1', 'Treatment2')
grp <- c('Treatment', 'Control')
colors <- setNames(c('#ff9289', '#4DBBD5'), grp)
pfx <- paste0(grp[1], '-vs-', grp[2], '_', collapse = '')

sam_tab <- read.delim('input/metadata1.txt', row.names = 2) %>%
  subset(Group %in% raw_grp)

if ('Control' %in% grp) {
  sam_tab$Group %<>% sub('Treatment[12]', 'Treatment', .)
}

# build phyloseq object ---------
seq_tab <- read.delim('input/jjh_16S_genus.txt', row.names = 1) %>%
  .[rownames(sam_tab)]

tax_tab <- data.frame(row.names = rownames(seq_tab), 
                      Genus = rownames(seq_tab))  %>%
  as.matrix

ps <- phyloseq(otu_table(seq_tab, taxa_are_rows = TRUE), tax_table(tax_tab),
               sample_data(sam_tab))

# group comparision ------------
ps_pair_test <- run_test_two_groups(
    ps,
    norm = 'CPM',
    group = "Group",
    method = "welch.test"
)# %>% marker_table()

ps_pair_test@marker_table$feature %<>% sub('.*__', '', .) 
marker_tab1 <- ps_pair_test %>% marker_table() %>% as.data.frame()
write.table(marker_tab1, file = paste0(pfx, 'diff-abundance_Genus_welch-test.txt'),
            sep = '\t', row.names = F, quote = F)

## heatmap
marker_tab1_order <- marker_tab1[order(-marker_tab1$ef_diff_mean), ]
markers <- marker_tab1_order$feature
coln_anno <- data.frame(row.names = markers,
                        Enrich_group = marker_tab1_order$enrich_group)
### column-scaled
pdf(paste0(pfx, 'diff-abundance_Genus_welch-test_heatmap_col-scaled.pdf'),
    6, 6, family = 'Times')
dat_scaled <- t(seq_tab[markers, ]) %>% scale(center = F)
pheatmap(dat_scaled,
         annotation_row = sam_tab,
         annotation_col = coln_anno,
         annotation_colors = list(
           Group = colors,
           Enrich_group = colors
         ),
         cluster_rows = F, cluster_cols = F,
         show_rownames = F,
         gaps_row = as.integer(table(sam_tab$Group)[1]))
dev.off()

### raw abundance
pdf(paste0(pfx, 'diff-abundance_Genus_welch-test_heatmap.pdf'), 6, 6,
    family = 'Times')
dat <- t(seq_tab[markers, ])
pheatmap(dat,
         annotation_row = sam_tab,
         annotation_col = coln_anno,
         annotation_colors = list(
           Group = colors,
           Enrich_group = colors
         ),
         cluster_rows = F, cluster_cols = F,
         show_rownames = F,
         gaps_row = as.integer(table(sam_tab$Group)[1]))
dev.off()

## boxplot
pdf(paste0(pfx, 'diff-abundance_Genus_welch-test_boxplot.pdf'), 5, 4,
    family = 'Times')
seq_tab[markers, ] %>%
  mutate(Genus0 = rownames(.)) %>%
  reshape2::melt(id.vars = 'Genus0') %>%
  mutate(
    Group = mapvalues(variable, from = rownames(sam_tab), to = sam_tab$Group),
    cols =  mapvalues(Genus0, from= marker_tab1_order$feature,
                      to=colors[marker_tab1_order$enrich_group]),
    Genus =  paste0("<span style=\"color: ", cols, "\">", Genus0, "</span>")) %>%
  mutate(Genus = factor(Genus, levels = mapvalues(markers,
                                                  from = Genus0, to = Genus)
    )) %>%
  ggplot(aes(value, Genus, fill = Group)) +
    geom_boxplot() +
    labs(x = 'Relative abundance') +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(panel.border = element_rect(size = 1),
          text = element_text(color = 'black'),
          axis.text.y = ggtext::element_markdown(),
          legend.position = c(.8, .15))

dev.off()
system(paste0('open ', pfx, 'diff-abundance_Genus_welch-test_boxplot.pdf'))


# LEfSe analysis -----------
ps_lefse <- run_lefse(
    ps = ps, 
    group = "Group",
    norm = "CPM",
   # transform = 'log10p',
    wilcoxon_cutoff = 0.05,
    kw_cutoff = 0.05,
    lda_cutoff = 2
)

lefse_dat <- ps_lefse %>% marker_table() %>% as.data.frame() 
lefse_dat$feature %<>% sub('.*__', '', .)

write.table(lefse_dat, file = paste0(pfx, 'LEfSe_result.txt'), sep = '\t',
                                     row.names =  F, quote = F)

## visual
lefse_dat2 <- lefse_dat %>%
  mutate(ef_lda_sign = ifelse(enrich_group == grp[2],
                              ef_lda, -ef_lda)) %>%
  .[order(.$ef_lda_sign), ] %>%
  mutate(color = colors[enrich_group],
         pos = ifelse(enrich_group == grp[2], 2, 4))

pdf(paste0(pfx, 'LEfSe_barplot.pdf'), 6, 6, family = 'Times')
barplot(lefse_dat2$ef_lda_sign, xlab="LDA SCORE", horiz=TRUE,
        space=0.25, width=0.8, col=lefse_dat2$color, names=FALSE,xlim = c(-6,6))

abline(v=setdiff(seq(-6,6,2), 0), col="lightgray", lty=2)

for (i in 1:nrow(lefse_dat2)) {
  text(0, i:i-0.4, lefse_dat2[i,]$feature, cex=0.8, offset=0.2, xpd=T, srt=0,
       pos = lefse_dat2[i,]$pos)
}

legend('topleft', legend = grp, pch=22,  col = NA, fill = colors, bty="n",
       cex=1)
dev.off()
system(paste0('open ', pfx, 'LEfSe_barplot.pdf'))
