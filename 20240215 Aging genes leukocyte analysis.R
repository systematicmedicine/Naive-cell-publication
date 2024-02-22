## 20240215 Aging genes leukocyte analysis.R
# Author: Cameron
# Date: 15th February 2024
#
# For the aging genes identified in Peters et al 2015, determine proportion of 
# expression attributable to each of the 18 leukocyte subtypes from the Human Blood 
# Atlas.
# 
# Peters: https://doi.org/10.1038/ncomms9570
# Human Blood Atlas: https://v19.proteinatlas.org/humanproteome/blood
#
# Script requires functions from 20240215 aging genes leukocyte analysis - functions.R

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load libraries
library(tidyverse)
library(openxlsx)
library(biomaRt)
library(cumstats)
library(ggpubr)
source("20240215 aging genes leukocyte analysis - functions.R")

## Define hard coded variables
PETERS_SUP1_PATH <- "41467_2015_BFncomms9570_MOESM436_ESM.xlsx"
BLOOD_ATLAS_PATH <- "rna_immune_cell.tsv"

## Load Peters supplementary data 1 (aging genes)
ptr <- read.xlsx(PETERS_SUP1_PATH, startRow=3) %>%
  dplyr::select(c("NEW-Gene-ID", "NEW-Entrez-ID", "RANK", "META:Direction")) %>%
  mutate(RANK = gsub("n", NA, RANK) %>% as.numeric(),
         `NEW-Entrez-ID` = as.character(`NEW-Entrez-ID`))

# Determine Ensembl ID using biomaRt
bm <- getBM(
  filters = "entrezgene_id",
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  values = unique(ptr$`NEW-Entrez-ID`),
  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl")))

ptr <- ptr %>%
  mutate(ensembl.id = bm$ensembl_gene_id[match(`NEW-Entrez-ID`, bm$entrezgene_id)], .after = `NEW-Entrez-ID`)  

## Load Human Blood Atlas (HBA)
hba <- read.table(BLOOD_ATLAS_PATH, sep = "\t", header = TRUE) %>%
  filter(Immune.cell != "total PBMC") %>%
  group_by(Gene) %>%
  mutate(prop.exp = nTPM/sum(nTPM)) %>%
  ungroup()

## Combine Peters and HBA datasets
comb <- hba %>%
  filter(Gene %in% ptr$ensembl.id | Gene.name %in% ptr$`NEW-Gene-ID`) %>%
  
  # Attempt to match by Ensembl ID
  mutate(
    rank = ptr$RANK[match(Gene, ptr$ensembl.id)],
    direction = ptr$`META:Direction`[match(Gene, ptr$ensembl.id)]
  ) %>%
  
  # For failed matches, attempt to match by Gene Symbol
  mutate(
    rank = ifelse(is.na(rank), ptr$RANK[match(Gene.name, ptr$`NEW-Gene-ID`)], rank),
    direction = ifelse(is.na(direction), ptr$`META:Direction`[match(Gene.name, ptr$`NEW-Gene-ID`)], direction)
  )

# Determine how many Peters aging genes successfully mapped
map.ptr.sig <- comb$rank %>% na.omit() %>% unique() %>% length()

## Loop through all 18 leukocyte types
fig1.lst <- list()
tbl1.lst <- list()

for (leuk in unique(comb$Immune.cell)){
  
  # All genes reported in Peters
  all.genes <- comb %>%
    filter(Immune.cell == leuk)
  all.prop <- median(all.genes$prop.exp, na.rm = TRUE)
  
  # All significant ageing genes reported in Peters
  sig.genes <- all.genes %>%
    filter(! is.na(rank))
  sig.prop <- median(sig.genes$prop.exp, na.rm = TRUE)
  sig.pval <- wilcox.test(x = na.omit(sig.genes$prop.exp), y = na.omit(all.genes$prop.exp))$p.value
  
  # Top 20 aging genes reported in Peters
  t20.genes <- all.genes %>%
    filter(rank <= 20) %>%
    arrange(rank)
  t20.prop <- median(t20.genes$prop.exp, na.rm = TRUE)
  t20.pval <- wilcox.test(x = na.omit(t20.genes$prop.exp), y = na.omit(all.genes$prop.exp))$p.value
  
  t20.plt <- ggplot(data = t20.genes, aes(x = rev(rank), y = prop.exp * 100, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("dodgerblue3", "darkgoldenrod")) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_discrete(limits = rev(factor(1:20)), breaks = c(1, 5, 10, 15, 20)) +
    labs(x = "Rank", y = paste("Expression attributable (%)")) +
    geom_abline(slope = 0, intercept = all.prop * 100, linetype = "dashed", colour = "red", linewidth=0.3) +
    ggtitle(leuk) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 8), 
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 12)
    ) +
    guides(fill = guide_legend(title = "Direction of expression change"))
 
  # Update table 1 list
  tbl1.lst[[length(tbl1.lst) + 1]] <- c(leuk, all.prop, sig.prop, sig.pval, t20.prop, t20.pval)
  
  # Update figure 1 list
  fig1.lst[[length(fig1.lst) + 1]] <- t20.plt
}

## Calculate stats for Naive-T group (CD4 + CD8)
naiveT.stats <- list()

# All genes reported in Peters
naiveT.all.genes <- comb %>%
  filter(Immune.cell == "naive CD8 T-cell" | Immune.cell == "naive CD4 T-cell") %>%
  group_by(Gene) %>%
  summarise(prop.exp = sum(prop.exp),
            rank = rank[1]) %>%
  ungroup()

naiveT.stats$all.prop <- median(naiveT.all.genes$prop.exp, na.rm = TRUE)

# All significant ageing genes reported in Peters
naiveT.sig.genes <- naiveT.all.genes %>%
  filter(! is.na(rank))

naiveT.stats$sig.prop <- median(naiveT.sig.genes$prop.exp, na.rm = TRUE)
naiveT.stats$sig.pval <- wilcox.test(x = na.omit(naiveT.sig.genes$prop.exp), y = na.omit(naiveT.all.genes$prop.exp))$p.value

# Top 20 aging genes reported in Peters
naiveT.t20.genes <- naiveT.all.genes %>%
  filter(rank <= 20)

naiveT.stats$t20.prop <- median(naiveT.t20.genes$prop.exp, na.rm = TRUE)
naiveT.stats$t20.pval <- wilcox.test(x = na.omit(naiveT.t20.genes$prop.exp), y = na.omit(naiveT.all.genes$prop.exp))$p.value
naiveT.stats$t20.fc <- naiveT.stats$t20.prop / naiveT.stats$all.prop

# Rank 1 and Rank 2
naiveT.stats$r1.pct <- naiveT.sig.genes$prop.exp[naiveT.sig.genes$rank == 1] * 100
naiveT.stats$r2.pct <- naiveT.sig.genes$prop.exp[naiveT.sig.genes$rank == 2] * 100

## Perform stats on direction of expression
direc.stats <- list()

direc.stats$sig.pos <- (ptr$`META:Direction`[! is.na(ptr$RANK)] == "+") %>% sum()
direc.stats$sig.neg <- (ptr$`META:Direction`[! is.na(ptr$RANK)] == "-") %>% sum()
direc.stats$sig.pval <- binom.test(direc.stats$sig.neg, direc.stats$sig.pos + direc.stats$sig.neg, p = 0.5)$p.value

direc.stats$t20.pos <- (ptr$`META:Direction`[ptr$RANK <= 20] == "+") %>% na.omit() %>% sum()
direc.stats$t20.neg <- (ptr$`META:Direction`[ptr$RANK <= 20] == "-") %>% na.omit %>% sum()
direc.stats$t20.pval <- binom.test(direc.stats$t20.neg, direc.stats$t20.pos + direc.stats$t20.neg, p = 0.5)$p.value

## Calculate cumulative median expression (plotted in fig 2)
cum.med <- comb %>%
  filter(! is.na(rank)) %>%
  arrange(rank) %>%
  dplyr::select(c(rank, Immune.cell, prop.exp)) %>%
  tidyr::spread(key = Immune.cell, value = prop.exp) %>%
  mutate(naiveT.group = `naive CD4 T-cell` + `naive CD8 T-cell`) %>%
  mutate_at(2:ncol(.), cummedian)

## Create table 1
tbl1 <- do.call(rbind.data.frame, tbl1.lst) %>%
  `colnames<-`(c("Immune.cell", "all.prop", "sig.prop", "sig.pval", "t20.prop", "t20.pval")) %>%
  mutate(Immune.cell = sapply(Immune.cell, fix_cell_name)) %>%
  mutate_at(2:ncol(.), as.numeric) %>%
  mutate(sig.logfc = log2(sig.prop / all.prop) %>% round(digits = 2), .after = sig.prop) %>%
  mutate(sig.padj = p.adjust(sig.pval, method="BH") %>% formatC(format = "e", digits = 2), .after = sig.pval) %>%
  mutate(t20.logfc = log2(t20.prop / all.prop) %>% round(digits = 2), .after = t20.prop) %>%
  mutate(t20.padj = p.adjust(t20.pval, method="BH") %>% formatC(format = "e", digits = 2), .after = t20.pval) %>%
  arrange(desc(t20.logfc))

# Write table 1 to disk
write.csv(tbl1, "manuscript/Table 1.csv", row.names = FALSE)

## Create figure 1

# Order by highest fold change
fig1.order <- match(tbl1$Immune.cell, comb$Immune.cell %>% unique() %>% fix_cell_name())
fig1.lst <- fig1.lst[fig1.order]

# Add stars to title
for (k in 1:length(fig1.lst)){
  pval <- tbl1$t20.padj[k] %>% as.numeric()
  title <- paste(tbl1$Immune.cell[k], assign_stars(pval))
  fig1.lst[[k]] <- fig1.lst[[k]] + ggtitle(title)
}

# Combine subplots
fig1 <- ggarrange(plotlist = fig1.lst,
                   ncol = 3, nrow = 6,
                   common.legend = TRUE,
                   legend = "bottom")

# Write figure 1 to disk
ggsave("manuscript/Figure 1.tiff", device = "tiff", plot = fig1, dpi = 300, width = 180, height = 270, unit = "mm")

## Create figure 2
fig2 <- ggplot(data = cum.med, aes(x = rank, y = naiveT.group * 100)) +
  geom_path(colour = "gray50", linewidth=1) +
  geom_abline(slope = 0, intercept = naiveT.stats$all.prop * 100, colour = "red", linetype = "dashed") +
  labs(y = "Median expression attributable to Naive T cells (%)", x = "Top n ranked genes") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_classic()

ggsave("manuscript/Figure 2.tiff", device = "tiff", fig2, dpi=300, width=180, height=120, unit = "mm")







