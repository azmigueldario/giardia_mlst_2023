# ===================================
# Dependencies

library(tidyverse)
library(janitor)
library(ggbreak)
library(svglite)

# ===================================
# Data cleaning

quast_report <- 
  read_tsv("./processed_data/quast/transposed_report.tsv", 
           na = c(" ", "", "-")) |> 
  clean_names()

# ===================================
# Exploratory plots 


# three extreme outliers with high N50
ggplot(data = quast_report,
       mapping = aes (x = "", y = n50)) +
  geom_violin(width = 0.7, color = "red")  +
  geom_boxplot(width = 0.1, color = "red4") +
  scale_y_break(breaks = c(600000, 1900000), scales = 0.2) +
  scale_y_break(breaks = c(125000, 500000), scales = 0.2) +
  labs(y = "N50 (basepairs)",
       x = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold",
                                    size = 15)) 

# distribution skewed to the right
ggplot(data = quast_report,
       mapping = aes (x = "", y = number_contigs)) +
  geom_violin(width = 0.7, color = "blue") +
  geom_boxplot(width = 0.1, color = "blue4") +
  scale_y_cut(breaks = 2000, space = 0.5, which = 1, scales = 0.5) +
  labs(y = "Number of Contigs",
       x = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold",
                                    size = 15))

# ===================================
# Filter poor quality genomes 

classified_data <- 
  quast_report |> 
  mutate(data_quality = case_when(n50 <= 30000 ~ "Poor",
                                  number_contigs >= 1300 ~ "Poor",
                                  total_length > 14500000 | total_length < 9660000 ~ "Poor",
                                  .default = "High")) |> 
  select(assembly, , data_quality, number_contigs, total_length, n50, 
         number_misassembled_contigs, genome_fraction_percent,
         number_ns_per_100_kbp) |> 
  pivot_longer(cols = c(-1, -2), names_to = "var", values_to = "value")

classified_data |> 
  filter(data_quality == "high") |> 
  select(assembly) |> 
  unique() |> 
  write_tsv("./output/selected_genomes.txt",
            col_names = F,
            append = F)

##########################
# Comparative plotting


# labels are uninformative and unclear
labels_clean <-   c("Ref. genome coverage", "N50", "Contig number", 
                    "Misassembled contigs (n)",
                    "N's per 100 kbp", "Genome length")
names(labels_clean) <- c("genome_fraction_percent", "n50", "number_contigs", "number_misassembled_contigs", 
                  "number_ns_per_100_kbp", "total_length")

# text and labels must be readable from afar
plot_quality <- 
  classified_data |> 
  ggplot(mapping = aes(x = "", y = value, color = data_quality)) +
  geom_jitter(size = 8) +
  facet_wrap(~var, 
             scales = "free_y", 
             nrow = 3,
             labeller = labeller(var = labels_clean)) +
  scale_color_manual(values = c("#6497bf", "#d8031c")) +
  theme_light() +
  labs(color = "Data quality \ncategory\n") +
  theme(strip.background=element_rect(colour="black",fill="grey30"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 22),
        legend.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 16)); plot_quality

ggsave("./output/bactopia_qc.svg",
       plot=plot_quality, 
       device = "svg",
       scale = 2,
       dpi = 600)

######################
# Palette colors

#6497bf	(100,151,191)
#9fcbee	(159,203,238)
#01016f	(1,1,111) Dark blue
#d8031c	(216,3,28)
#5a5a5a
#4BBF64
