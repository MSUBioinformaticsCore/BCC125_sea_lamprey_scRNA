# Stage composition along the germ cell pseudotime trajectory
#
# Purpose: quantify whether larval and transformer cells span overlapping ranges
# of pseudotime, or segregate to opposite ends. If they overlap, the trajectory
# is not reducible to the July/November batch axis, and this is the single most
# useful figure for reviewer comment 3.
#
# Adjust the object and column names in the CONFIG block to match your session.

library(SingleCellExperiment)
library(TSCAN)
library(ggplot2)
library(dplyr)
library(tidyr)

## ---- CONFIG -----------------------------------------------------------------

sce      <- sce.germ          # germ cell subclustered SCE
pt_out   <- out               # object returned by quickPseudotime()
sample_col <- "Sample"        # colData column holding the library ID
larval_ids <- c("C00004336", "C00004362")   # NB: resolve C00003262 vs C00004362

outdir <- "."

## ---- ASSEMBLE ---------------------------------------------------------------

pt <- as.data.frame(pt_out$ordering)
colnames(pt) <- paste0("path", seq_len(ncol(pt)))

df <- data.frame(
  cell    = colnames(sce),
  sample  = as.character(colData(sce)[[sample_col]]),
  cluster = as.character(colLabels(sce)),
  pt,
  check.names = FALSE
) %>%
  mutate(
    stage = ifelse(sample %in% larval_ids, "larva", "transformer"),
    library = case_when(
      sample == "C00004337" ~ "female transformer",
      sample == "C00004374" ~ "male transformer",
      TRUE ~ "larva"
    )
  )

long <- df %>%
  pivot_longer(starts_with("path"), names_to = "path", values_to = "pseudotime") %>%
  filter(!is.na(pseudotime))

## ---- 1. CLUSTER COMPOSITION -------------------------------------------------
## Simplest and most convincing evidence. If germ cell clusters contain cells from
## both stages, fastMNN is aligning genuinely shared populations rather than
## inventing correspondence.

comp <- df %>%
  count(cluster, stage) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

comp_wide <- comp %>%
  select(-n) %>%
  pivot_wider(names_from = stage, values_from = prop, values_fill = 0) %>%
  arrange(desc(larva))

print(comp_wide, n = Inf)
write.csv(comp_wide, file.path(outdir, "cluster_stage_composition.csv"), row.names = FALSE)

# how many clusters are genuinely mixed (both stages >= 10%)?
mixed <- comp_wide %>% filter(larva >= 0.10, transformer >= 0.10)
cat("\nClusters with >=10% of cells from both stages:", nrow(mixed),
    "of", nrow(comp_wide), "\n")

## ---- 2. PSEUDOTIME RANGE BY STAGE, PER BRANCH -------------------------------
## The numbers to quote in the response.

rng <- long %>%
  group_by(path, stage) %>%
  summarise(
    n      = n(),
    min    = min(pseudotime),
    q10    = quantile(pseudotime, 0.10),
    median = median(pseudotime),
    q90    = quantile(pseudotime, 0.90),
    max    = max(pseudotime),
    .groups = "drop"
  )
print(rng, n = Inf)
write.csv(rng, file.path(outdir, "pseudotime_range_by_stage.csv"), row.names = FALSE)

# Overlap of the central 80% of each stage's distribution, per branch.
# A large overlap is the result that weakens the reviewer's objection.
ov <- rng %>%
  select(path, stage, q10, q90) %>%
  pivot_wider(names_from = stage, values_from = c(q10, q90)) %>%
  mutate(
    overlap_lo = pmax(q10_larva, q10_transformer),
    overlap_hi = pmin(q90_larva, q90_transformer),
    overlap_width = pmax(0, overlap_hi - overlap_lo),
    union_width   = pmax(q90_larva, q90_transformer) - pmin(q10_larva, q10_transformer),
    overlap_frac  = overlap_width / union_width
  )
print(ov, n = Inf)

## ---- 3. STAGE PROPORTION IN PSEUDOTIME DECILES ------------------------------

dec <- long %>%
  group_by(path) %>%
  mutate(decile = ntile(pseudotime, 10)) %>%
  count(path, decile, stage) %>%
  group_by(path, decile) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write.csv(dec, file.path(outdir, "stage_prop_by_pseudotime_decile.csv"), row.names = FALSE)

## ---- 4. FIGURES -------------------------------------------------------------

p1 <- ggplot(long, aes(pseudotime, fill = stage)) +
  geom_density(alpha = 0.55, colour = NA) +
  facet_wrap(~ path, scales = "free", ncol = 1) +
  labs(x = "Pseudotime", y = "Density",
       title = "Stage composition along each trajectory branch") +
  theme_bw(base_size = 11)

p2 <- ggplot(dec, aes(factor(decile), prop, fill = stage)) +
  geom_col(width = 0.9) +
  facet_wrap(~ path, ncol = 1) +
  labs(x = "Pseudotime decile", y = "Proportion of cells",
       title = "Proportion of cells from each stage, by pseudotime decile") +
  theme_bw(base_size = 11)

p3 <- ggplot(comp, aes(reorder(cluster, -prop * (stage == "larva")), prop, fill = stage)) +
  geom_col(width = 0.9) +
  labs(x = "Germ cell subcluster", y = "Proportion of cells",
       title = "Stage composition of each germ cell subcluster") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(outdir, "pseudotime_stage_density.pdf"), p1, width = 7, height = 7)
ggsave(file.path(outdir, "pseudotime_stage_deciles.pdf"), p2, width = 7, height = 7)
ggsave(file.path(outdir, "cluster_stage_composition.pdf"), p3, width = 7, height = 4)

## ---- 5. SENTENCE FOR THE RESPONSE -------------------------------------------

for (pp in unique(dec$path)) {
  d <- dec %>% filter(path == pp, stage == "larva")
  mixed_dec <- sum(d$prop > 0.10 & d$prop < 0.90)
  cat(sprintf(
    "\n%s: %d of 10 pseudotime deciles contain between 10%% and 90%% larval cells.\n",
    pp, mixed_dec))
}
cat("\nIf most deciles are mixed, state that directly with these numbers.\n",
    "If the early deciles are ~100% larval and the late ~100% transformer,\n",
    "the ordering cannot be distinguished from the batch axis and the claim\n",
    "must be softened rather than defended.\n")
