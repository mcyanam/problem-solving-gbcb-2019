#Normalization and differential expression using NOISeq

library(ggplot2)
library(gplots)
library(NOISeq)

#Fibrobacter succinogenes
fscount = read.csv(“Fs.csv”, row.names = 1)
fsfactors = read.csv(“Fs_metadata.csv”)
fsdata = readData(data = fscount, factors = fsfactors)
fsnoiseq = noiseq(fsdata, k = NULL, factor = “condition”, norm = “tmm”, replicates = “no”, pnr = 0.2, nss = 5, v = 0.02, lc = 0)
fsnoiseq@results[[1]]
write.csv(fsnoiseq@results, file = “fs_normalized_data.csv”)


###Upregulated and downregulated genes:
fsdata.high = degenes(fsnoiseq, q = 0.8, M = “up”)
write.csv(fsdata.up, file = “DE_fs.highmethane.csv”)
fsdata.low = degenes(fsnoiseq, q = 0.8, M = “down”)
write.csv(fsdata.down, file = “DE_fs.lowmethane.csv”)

#Prevotella ruminicola
prcount = read.csv(“Pr_count.csv”, row.names = 1)
prfactors = read.csv(“Pr_metadata.csv”)
prdata = readData(data = prcount, factors = prfactors)
prnoiseq = noiseq(prdata, k = NULL, factor = “condition”, norm = “tmm”, replicates = “no”, pnr = 0.2, nss = 5, v = 0.02, lc = 0)
prnoiseq@results[[1]]
write.csv(prnoiseq@results, file = “pr_normalized_data.csv”)


###Upregulated and downregulated genes:
prdata.high = degenes(prnoiseq, q = 0.8, M = “up”)
write.csv(prdata.up, file = “DE_pr.highmethane.csv”)
prdata.low = degenes(prnoiseq, q = 0.8, M = “down”)
write.csv(prdata.down, file = “DE_pr.lowmethane.csv”)

#Treponema succinifaciens
tscount = read.csv(“Ts_count.csv”, row.names = 1)
tsfactors = read.csv(“Ts_metadata.csv”)
tsdata = readData(data = tscount, factors = tsfactors)
tsnoiseq = noiseq(tsdata, k = NULL, factor = “condition”, norm = “tmm”, replicates = “no”, pnr = 0.2, nss = 5, v = 0.02, lc = 0)
tsnoiseq@results[[1]]
write.csv(tsnoiseq@results, file = “ts_normalized_data.csv”)


###Upregulated and downregulated genes:
fsdata.high = degenes(tsnoiseq, q = 0.8, M = “up”)
write.csv(tsdata.up, file = “DE_ts.highmethane.csv”)
tsdata.low = degenes(tsnoiseq, q = 0.8, M = “down”)
write.csv(tsdata.down, file = “DE_ts.lowmethane.csv”)

#Slackia heliotrinireducens
shcount = read.csv(“Sh_count.csv”, row.names = 1)
shfactors = read.csv(“Sh_metadata.csv”)
shdata = readData(data = shcount, factors = shfactors)
shnoiseq = noiseq(shdata, k = NULL, factor = “condition”, norm = “tmm”, replicates = “no”, pnr = 0.2, nss = 5, v = 0.02, lc = 0)
shnoiseq@results[[1]]
write.csv(shnoiseq@results, file = “sh_normalized_data.csv”)


###Upregulated and downregulated genes:
shdata.high = degenes(shnoiseq, q = 0.8, M = “up”)
write.csv(shdata.up, file = “DE_sh.highmethane.csv”)
shdata.low = degenes(shnoiseq, q = 0.8, M = “down”)
write.csv(shdata.down, file = “DE_sh.lowmethane.csv”)


#Methanobrevibacter ruminantium
mrcount = read.csv(“Mr_count.csv”, row.names = 1)
mrfactors = read.csv(“Mr_metadata.csv”)
mrdata = readData(data = mrcount, factors = mrfactors)
mrnoiseq = noiseq(mrdata, k = NULL, factor = “condition”, norm = “tmm”, replicates = “no”, pnr = 0.2, nss = 5, v = 0.02, lc = 0)
mrnoiseq@results[[1]]
write.csv(mrnoiseq@results, file = “mr_normalized_data.csv”)


###Upregulated and downregulated genes:
mrdata.high = degenes(mroiseq, q = 0.8, M = “up”)
write.csv(mrdata.up, file = “DE_mr.highmethane.csv”)
mrdata.low = degenes(mrnoiseq, q = 0.8, M = “down”)
write.csv(mrdata.down, file = “DE_mr.lowmethane.csv”)



#Taxonomic profile based on OTU count matrix and taxa tables

library(phyloseq)
library(ape)
library(vegan)
otumat = read.csv(“OTU_table.csv”, row.names = 1)
otumat = as.matrix(otumat)
OTU = otu_table(otumat, taxa_as_rows = TRUE)

taxmat = read.csv(“TAX_table.csv”, row.names = 1)
taxmat = as.matrix(taxmat)
TAX = tax_table(taxmat)

sample1 = read.csv(“all_samples.csv”, row.names = 1)
sampledata = sample_data(sample1)

physeq = phyloseq(OTU, TAX, sampledata)

#Making bar plot with alpha diversity:
plot_richness(physeq, color = “Samples”) + geom_point(size = 3) 

#Making bar plot with compartmentalized:
plot_bar(physeq, “Kingdom”, fill = “Phylum”, facet_grid = ~Sample)

#Making heatmap:
plot_heatmap(physeq)


#Heatmap using ggplot or ggplot2
library(gplofs)
library(ggplot2)
library(vegan)
library(Heatplus)
library(RColorBrewer)

#Following the data that we already have as otumat (otu table that we used before):
#Make it based on abundance  
data.prop = otumat/rowSums(otumat)
data.dist = vegdist(data.prop, method = “bray” or “jaccard”)
row.clus = hclust(data.dist, “aver”)
data.dist.g = vegdist(t(data.prop), method = “bray”)
col.clus = hclust(data.dist.g, “aver”)
heatmap.2(as.matrix(data.prop), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = brewer.pal(6,”Accent”), scale = “row”, trace = “none”, margins = c(7,10), colsep = 1:3, sepwidth = c(0.001, 0.001))
