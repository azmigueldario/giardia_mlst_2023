# Scripts to download and pre-process dataset

1. `download_data_repositories.sh`: Uses nf-core fetchngs and curl to download reference genome and all available NGS data for Giardia duodenalis assemblages A and B
2. `assembly_and_qc.sh`:  Runs bactopia/bactopia pipeline for all accessions with short-read data. As _Giardia spp._ genome is mostly redundant, we can use bacterial sequencing assembly tools
2. `hybrid_assembly.sh`: For accessions with short and long sequencing reads, we perform hybrid assembly in Bactopia. 
