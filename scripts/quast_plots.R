# ===================================
# Dependencies

library(tidyverse)
library(janitor)
library(ggbreak)

# ===================================
# Data cleaning

quast_report <- 
  read_tsv("./quast/transposed_report.tsv", 
           na = c(" ", "", "-")) |> 
  clean_names()


quast_report<- 
  tibble(quast_report,
         .name_repair = make_clean_names)

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

quast_report |> 
  filter (n50 > 30000 & number_contigs < 1500) |> 
  select(assembly) |> 
  write_tsv("./output/selected_genomes.txt",
            col_names = F,
            append = T)

