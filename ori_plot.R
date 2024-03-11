##### Analysing origin of Bcell clones V1
rm(list = ls())

library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

trees <- readRDS("./fake_treedata.rds")
#function to filter for only trees with >= 1 synovium + Week16 clone
filter_family <- function(tree_list) {
  sel <- function(sq) {
    keep <- sum(str_detect(ifelse(sq$tissue == "synovium" & sq$time == "Week16",
                                  "TRUE", "FALSE"),
                           "TRUE"))
    keep <- ifelse(keep > 0, TRUE, FALSE)
    return(keep)
  }
  patient <- tree_list
  to_keep <- map_lgl(patient$Sequences, sel)
  keep <- patient[to_keep, ]
  return(keep)
}
tree_syn_16w <- map(trees, filter_family)
remove_after_asterisk <- function(df) {
  df$Sequences <- map(df$Sequences,
                      ~mutate(., J.allele = sub("^(.*?)\\*.*", "\\1", J.allele),
                      V.allele = sub("^(.*?)\\*.*", "\\1", V.allele)
  ))
  return(df)
}
tree_syn_16w <- map(tree_syn_16w, remove_after_asterisk)
#extract the information what the most ancestral(root) node is
roots <- function(tree, sequ) {
  tree_tibble <- as_tibble(tree)
  tree <- drop.tip(tree, tip = "Germline")
  root_number <- find_root(tree)
  tree_tibble[which(tree_tibble$node == root_number), "label"] <-
    tree_tibble[which(tree_tibble$node == 1), "label"]
  rootID <- as.character(tree_tibble[which(tree_tibble$node == 1), "label"])
  rootID <- gsub("ID_", "", rootID)
  sequ <- sequ %>% filter(Clone.ID == rootID)
  sequ$V.allele <- gsub("IGHV", "V", sequ$V.allele)
return(sequ)
}
origins <- map2_df(tree_syn_16w, names(tree_syn_16w),
                   function(patient, name) {
                     tree <- patient$Tree
                     Seq <- patient$Sequences
                     ori <- map2_df(tree, Seq, roots)
                    ori$sample <- name
return(ori)
})
## count of origin of tissue by V_allele
v_origins <- origins %>%
  group_by(V.allele, tissue) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
# Order V.allele by the sum of count grouped by V.allele
ordered_alleles <- v_origins %>%
  group_by(V.allele) %>%
  summarise(sum_count = sum(count)) %>%
  arrange(desc(sum_count)) %>%
  pull(V.allele)
# Reorder V.allele factor levels
v_origins$V.allele <- factor(v_origins$V.allele, levels = ordered_alleles)
# Create stacked bar plot
ggplot(v_origins, aes(x = V.allele, y = count, fill = factor(tissue))) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = count, fill = factor(tissue)), stat = "identity") +
  labs(x = "V gene", y = "Total trees", fill = "Percentage") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
## separate for response/non-response
## Responders: Patient 1,2 Non-Responders Patient 3,4
origins_res <- map_df(tree_syn_16w[c(1, 2)], function(patient) {
  tree <- patient$Tree
  Seq <- patient$Sequences
  ori <- map2_df(tree, Seq, roots)
  return(ori)
})
origins_nonres <- map_df(tree_syn_16w[c(3, 4)], function(patient) {
  tree <- patient$Tree
  Seq <- patient$Sequences
  ori <- map2_df(tree, Seq, roots)
  return(ori)
})
# which V.alleles are only in either group 
v_genes_res <- unique(origins_res$V.allele)
v_genes_nonres <- unique(origins_nonres$V.allele)
setdiff(v_genes_nonres, v_genes_res)
setdiff(v_genes_res, v_genes_nonres)
## count of origin of tissue by V_allele in responders
v_origins_res <- origins_res %>%
  group_by(V.allele, tissue) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
## count of origin of tissue by V_allele in non responders
v_origins_nonres <- origins_nonres %>%
  group_by(V.allele, tissue) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
# Order V.allele by the sum of count grouped by V.allele
ordered_alleles_res <- v_origins_res %>%
  group_by(V.allele) %>%
  summarise(sum_count = sum(count)) %>%
  arrange(desc(sum_count)) %>%
  pull(V.allele)
# Order V.allele by the sum of count grouped by V.allele
ordered_alleles_nonres <- v_origins_nonres %>%
  group_by(V.allele) %>%
  summarise(sum_count = sum(count)) %>%
  arrange(desc(sum_count)) %>%
  pull(V.allele)
# Reorder V.allele factor levels
v_origins_res$V.allele <- factor(v_origins_res$V.allele,
                                 levels = ordered_alleles_res)
v_origins_nonres$V.allele <- factor(v_origins_nonres$V.allele,
                                    levels = ordered_alleles_nonres)
v_origins_res$Response <- "R"
v_origins_nonres$Response <- "NR"
v_origins_combined <- rbind(v_origins_res, v_origins_nonres)
# Create stacked bar plot response
ggplot(v_origins_res, aes(x = V.allele, y = count, fill = factor(tissue))) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = count, fill = factor(tissue)), stat = "identity") +
  labs(x = "V gene", y = "Total trees", fill = "Percentage") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
# Create stacked bar plot non response
ggplot(v_origins_nonres, aes(x = V.allele, y = count, fill = factor(tissue))) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = count, fill = factor(tissue)), stat = "identity") +
  labs(x = "V gene", y = "Total trees", fill = "Percentage") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
barwidth <- 0.35
v_origins_nonres$V.allele <- factor(v_origins_nonres$V.allele,
                                    levels = ordered_alleles_res)
## This was an attempt to plot it stacked, dodged
## but the V.allele not being numeric, I couldn"t find a solution
# ggplot(v_origins_res,
#        aes(x = V.allele, y = count, fill = factor(tissue))) +
#   geom_bar(stat = "identity") +
#   geom_bar(aes(y = count,
#                fill = factor(tissue)), stat = "identity") +
#   labs(x = "V gene", y = "Total trees", fill = "Percentage") +
#   scale_fill_manual(values = c("red", "blue")) +
#   theme_minimal() +
#   theme(legend.position = c(0.85, 0.8),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_bar(data = v_origins_nonres, stat = "identity",
#            mapping = aes(x = (V.allele + barwidth + 0.01),
#                          y = count, fill = as.factor(tissue))) +
#   geom_bar(aes(y = count,
#                fill = factor(tissue)), stat = "identity") +
#   labs(x = "V gene", y = "Total trees", fill = "Percentage") +
#   scale_fill_manual(values = c("red", "blue")) +
#   theme_minimal() +
#   theme(legend.position = c(0.85, 0.8),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
v_origins_combined$V.allele <- factor(v_origins_combined$V.allele,
                                    levels = ordered_alleles_res)
df <- v_origins_combined %>%
  mutate(x_label = factor(paste0(V.allele, " ", Response)))
# Create stacked bar plot
ggplot(df, aes(x = x_label, y = count, fill = factor(tissue))) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "V gene", y = "Total trees", fill = "Percentage") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.2))
