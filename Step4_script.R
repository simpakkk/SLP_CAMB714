# the essentials ----
library(tidyverse)
library(rhdf5)
library(edgeR)

archs4.mouse <- "mouse_matrix_v10.h5" # if you placed the hdf5 file in your working directory, just use "human_matrix_v8.h5" as the path
all.samples.mouse <- h5read(archs4.mouse, name="meta/samples/geo_accession")
mySamples <- c("GSM2310941", # WT_unstim_rep1
               "GSM2310942", # WT_unstim_rep2
               "GSM2310943", # Ripk3_unstim_rep1
               "GSM2310944", # Ripk3_unstim_rep2
               "GSM2310945", # Ripk3Casp8_unstim_rep1
               "GSM2310946", # Ripk3Casp8_unstim_rep2
               "GSM2310947", # WT_LPS.6hr_rep1
               "GSM2310948", # WT_LPS.6hr_rep2
               "GSM2310949", # Ripk3_LPS.6hr_rep1
               "GSM2310950", # Ripk3_LPS.6hr_rep2
               "GSM2310951", # Ripk3Casp8_LPS.6hr_rep1
               "GSM2310952") # Ripk3Casp8_LPS.6hr_rep2

my.sample.locations <- which(all.samples.mouse %in% mySamples)
genes <- h5read(archs4.mouse, "meta/genes/gene_symbol")
expression <- h5read(archs4.mouse, "data/expression", 
                     index=list(my.sample.locations, 1:length(genes)))

expression <- t(expression)
rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)

keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")
archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

sample_source_name <- h5read(archs4.mouse, "meta/samples/source_name_ch1")
sample_title <- h5read(archs4.mouse, name="meta/samples/title")
sample_characteristics<- h5read(archs4.mouse, name="meta/samples/characteristics_ch1")

studyDesign <- tibble(Sample_title = sample_title[my.sample.locations], 
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations], 
                      genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                      treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=treatment, shape=genotype) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

