### Creating FAB Operon Conservation Image 20251217 ###

# ============================================================
# 1. Setup- set directory and load necessary packages
# ============================================================
setwd("C:/Users/ceeja/OneDrive/Desktop/Lactobacilli_family_tree_20250505/fab_conservation_trees")

library(ggplot2)
library(gggenes)
library(rtracklayer)
library(dplyr)
library(readr)

# ============================================================
# 2. Read in the GFF file downloaded from geneious
# ============================================================
gff <- import("36_fab_operons.gff")

# convert to dataframe
gff_df <- as.data.frame(gff)

# ============================================================
# 3. Define FAB gene colors
# ============================================================
gene_colors <- c(
  fabI = "#0D1164",  accA = "#8C1007",  accD = "#B07AA1",
  accC = "#FFC0CB",  fabZ = "#FFFF00",  accB = "#06D001",
  fabF = "#E15759",  fabG = "#640D5F",  fabD = "#EA2264",
  fabK = "#3ABEF9",  acpP = "#F28E2B",  fabH = "#59A14F",
  fabT = "#4E79A7",  other = "grey70")

# make lowercase
names(gene_colors) <- tolower(names(gene_colors))

# ============================================================
# 4. Build gene table from GFF- extract info and clean names
# ============================================================
genes_df <- gff_df %>%
  filter(type == "CDS") %>%
  transmute(
    molecule = seqnames,
    Name_raw = Name,
    start    = start,
    end      = end,
    strand   = strand,
    forward  = strand == "+"
  )

# Clean molecule (species) names (for plotting)
genes_df$molecule <- genes_df$molecule |>
  gsub("%20", " ", x = _) |>
  gsub(" \\(reversed\\)", "", x = _) |>
  gsub("_", " ", x = _)

# ============================================================
# 5. Read CURATED mapping table 
# ============================================================
mapping_final <- read_tsv(
  "fab_gene_annotation_mapping_suggested.tsv",
  col_types = cols(.default = "c")
)

# clean up 
normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>  # removes hidden unicode
    trimws() |>
    gsub("\\s+", " ", x = _) |>
    tolower()
}
# Apply Mapping Table
genes_df <- genes_df %>%
  left_join(mapping_final, by = c("Name_raw")) %>%
  mutate(
    gene_clean = tolower(trimws(gene_final))
  ) %>%
  filter(!is.na(gene_clean))

# SAFETY CHECK 
stopifnot(all(unique(genes_df$gene_clean) %in% names(gene_colors)))

# ============================================================
# 6. Define Operon Order based on ANI tree
# ============================================================
operon_order <- c(
  "Lacticaseibacillus paracasei subsp tolerans",
  "Lacticaseibacillus rhamnosus",
  "Lacticaseibacillus chiayiensis",
  "Lacticaseibacillus zeae",
  "Lacticaseibacillus casei",
  "Lacticaseibacillus huelsenbergensis",
  "Lacticaseibacillus parahuelsenbergensis",
  "Lacticaseibacillus styriensis",
  "Lactobacillus delbrueckii subsp.bulgaricus",
  "Lactobacillus delbrueckii subsp.indicus",
  "Lactobacillus delbrueckii subsp lactis",
  "Lactobacillus delbrueckii subsp.delbrueckii",
  "Lactobacillus delbrueckii subsp.jakobsenii",
  "Lactobacillus acetotolerans",
  "Lactobacillus amylolyticus",
  "Lactobacillus crispatus" ,
  "Lactobacillus amylovorus",
  "Lactobacillus helveticus",
  "Lactobacillus kefiranofaciens",
  "Lactiplantibacillus pentosus",
  "Lactiplantibacillus paraplantarum",
  "Lactiplantibacillus argentoratensis",
  "Lactiplantibacillus plantarum",
  "Limosilactobacillus reuteri",
  "Limosilactobacillus portuensis",
  "Limosilactobacillus vaginalis",
  "Ligilactobacillus animalis",
  "Ligilactobacillus murinus",
  "Limosilactobacillus mucosae",
  "Limosilactobacillus fermentum",
  "Ligilactobacillus salvarius",
  "Ligilactobacillus ruminis",
  "Ligilactobacillus agilis",
  "Ligilactobacillus acidipiscis",
  "Lacticaseibacillus pabuli",
  "Lacticaseibacillus pantheris"
)
# apply order
genes_df$molecule <- factor(
  genes_df$molecule,
  levels = rev(operon_order)  # rev() puts first entry at the TOP
)

# ============================================================
# 7. Anchor on fabF
# ============================================================
fabF_anchor <- genes_df %>%
  filter(gene_clean == "fabf") %>%
  group_by(molecule) %>%
  summarise(
    fabF_mid = (start + end) / 2,
    .groups = "drop"
  )

genes_df <- genes_df %>%
  left_join(fabF_anchor, by = "molecule") %>%
  mutate(
    start_aligned = start - fabF_mid,
    end_aligned   = end   - fabF_mid
  )

# ============================================================
# 8. Plot 
# ============================================================
p <-ggplot(genes_df,
    aes(
    xmin = start_aligned,
    xmax = end_aligned,
    y = molecule,
    fill = gene_clean,
    forward = forward)) +
  geom_gene_arrow() +
  scale_fill_manual(values = gene_colors,
                    breaks = c(
                      "fabz", "fabt", "fabh", "acpp", "fabk",
                      "fabd", "fabg", "fabf", "accb", "accc", 
                      "accd", "acca", "fabi", "other"),
                    labels = c(
                      "FabZ", "FabT", "FabH", "AcpP", "FabK", "FabD",
                      "FabG", "FabF", "AccB", "AccC", "AccD",
                      "AccA", "FabI", "Other")) +
  theme_genes() +
  theme(
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.title.y = element_blank(),
    legend.position = "right") +
  labs(
    x = "Position relative to fabF (bp)",
    fill = "FAB operon gene"
  )

print(p)

# ============================================================
# 9. Save Plot as PDF; choose size
# ============================================================
ggsave("36_fab_operon_10x10plot20251217.pdf", plot = p, width = 10, height = 10, units = "in")
ggsave("36_fab_operon_10x8plot20251217.pdf", plot = p, width = 8, height = 10, units = "in")
ggsave("36_fab_operon_10x7plot20251217.pdf", plot = p, width = 7, height = 10, units = "in")
ggsave("36_fab_operon_10x6plot20251217.pdf", plot = p, width = 6, height = 10, units = "in")
