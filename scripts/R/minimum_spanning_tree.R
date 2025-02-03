# Minimum spanning tree

## Required modules
library(tidyverse)
library(igraph)
library(ape)
library(ggnetwork)

###################################
# Load data -----------------------

cgmlst_alleles <- read_tsv("./processed_data/chewBBACA/cgMLST_call.tsv")

cgmlst_metadata <- 
  read_tsv("./processed_data/metadata/giardia_metadata.tsv") |> 
  filter(study_accession == "PRJNA280606") |> 
  mutate(label = paste0(location, "_", collection_date_year),
         run_alias = str_replace_all(run_alias, "/", "_")) |> 
  select(id, location) 
  


trial_mst <- left_join(cgmlst_alleles, cgmlst_metadata, join_by(FILE == id)) 

trial_mst |> 
  relocate(location, .after = 1) |>
  mutate(location =  replace_na(location, "other")) |> 
  select(FILE, location) |> 
  write_tsv("./output/trial_MST_metadata.tsv")


cgmlst_hamming <- 
  read_tsv("./scripts/cgmlst_hamming.tsv")

###################################
# Calculate distances ------------


modded_alleles <- cgmlst_alleles |> 
  column_to_rownames(var = "FILE")  
cgmlst_dist  <- ape::dist.gene(modded_alleles) |> as.matrix()

<- ape::dist.gene(cgmlst_alleles) 

# left_join(cgmlst_data)

cgmlst_metadata$run_alias 
cgmlst_alleles$FILE

  # Data must be transformed to reflect distances among observations
cgmlst_dist <- dist(cgmlst_data, )



# Plot data -----------------------

cgmlst_graph <- 
  ape::mst(cgmlst_dist) |> 
  graph_from_adjacency_matrix(adjmatrix = _, mode = "undirected")
 
ggraph(cgmlst_graph, layout)
attr(cgmlst_dist, "Label") <- cgmlst_metadata$label


ggnetwork(cgmlst_graph,
          arrow.gap = 0,
          layout = layout_with_fr(cgmlst_graph)) |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = name), size = 2) +  theme_blank() +
  geom_nodetext(aes(label = name), color = "black", size = 2.5) +
  theme(plot.margin = unit(c(0, 1, 1, 6), "cm"))+
  theme(legend.position = c(0, 0.14),
        legend.background = element_blank(),
        legend.text = element_text(size = 7))

# Example MST -----------------------
  
load(url("https://web.stanford.edu/class/bios221/data/dist2009c.RData"))
dist2009c["Label"]
country09 = attr(dist2009c, "Label")
mstree2009 = ape::mst(dist2009c)
gr09 = graph.adjacency(mstree2009, mode = "undirected")
gg = ggnetwork(gr09, arrow.gap = 0, layout = layout_with_fr(gr09))
ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = name), size = 2) +  theme_blank() +
  geom_nodetext(aes(label = name), color = "black", size = 2.5) +
  theme(plot.margin = unit(c(0, 1, 1, 6), "cm"))+
  theme(legend.position = c(0, 0.14),
        legend.background = element_blank(),
        legend.text = element_text(size = 7))