################# INSTALL REQUIRED PACKAGES ############################

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("KEGGgraph")
BiocManager::install("RBGL")

################# LOAD PACKAGES #######################################
library(KEGGgraph)
library(RBGL)

################ SPECIES SPECIFIC REACTIONS ###########################

############################### MST #####################################
mst_files <- c("mst00010.xml","mst00020.xml","mst00030.xml",
               "mst00040.xml","mst00051.xml","mst00052.xml",
               "mst00500.xml","mst000511.xml","mst00620.xml",
               "mst00680.xml","mst01200.xml")
files_mst <- system.file(paste("/extdata/", mst_files, sep=""), package="KEGGgraph")
# parse KGML files 
kgml_mst <- sapply(files_mst, parseKGML)
# make individual reaction graphs
reaction_mst <- sapply(kgml_mst, KEGGpathway2reactionGraph)
# merge reaction graphs
pan_mst <- mergeGraphs(reaction_mst)
# plot merged reaction graph
mst_plot <- plot(pan_mst,attrs=list(node=list(label="foo", fillcolor="lightgreen"),edge=list(color="black")))

############################### MRU #####################################
mru_files <- c("mru00010.xml","mru00020.xml","mru00030.xml",
               "mru00040.xml","mru00051.xml","mru00052.xml",
               "mru00061.xml","mru00061.xml","mru00500.xml",
               "mru000511.xml","mru00620.xml",
               "mru00680.xml","mru01200.xml")
files_mru <- system.file(paste("/extdata/", mru_files, sep=""), package="KEGGgraph")
kgml_mru <- sapply(files_mru, parseKGML)
reaction_mru <- sapply(kgml_mru, KEGGpathway2reactionGraph)
pan_mru <- mergeGraphs(reaction_mru)
mru_plot <- plot(pan_mru,attrs=list(node=list(label="foo", fillcolor="lightgreen"),edge=list(color="black")))

############################### FSU ###################################

fsu_files <- c("fsu00010.xml","fsu00020.xml","fsu00030.xml",
           "fsu00040.xml","fsu00051.xml","fsu00052.xml",
           "fsu00061.xml","fsu00061.xml","fsu00500.xml",
           "fsu000511.xml","fsu00620.xml",
           "fsu00680.xml","fsu01200.xml")
files_fsu <- system.file(paste("/extdata/", fsu_files, sep=""), package="KEGGgraph")
kgml_fsu <- sapply(files_fsu, parseKGML)
reaction_fsu <- sapply(kgml_fsu, KEGGpathway2reactionGraph)
pan_fsu <- mergeGraphs(reaction_fsu)
plot(pan_fsu,attrs=list(node=list(label="foo", fillcolor="lightgreen"),edge=list(color="black")))

################################# PRU ###################################

pru_files <- c("pru00010.xml","pru00020.xml","pru00030.xml",
               "pru00040.xml","pru00051.xml","pru00052.xml",
               "pru00061.xml","pru00061.xml","pru00500.xml",
               "pru000511.xml","pru00620.xml",
               "pru00680.xml","pru01200.xml")

files_pru <- system.file(paste("/extdata/", pru_files, sep=""), package="KEGGgraph")
kgml_pru <- sapply(files_pru, parseKGML)
reaction_pru <- sapply(kgml_pru, KEGGpathway2reactionGraph)
pan_pru <- mergeGraphs(reaction_pru)
pru_plot <- plot(pan_pru,attrs=list(node=list(label="foo", fillcolor="lightgreen"),edge=list(color="black")))

################################## TSU  #################################

tsu_files <- c("tsu00010.xml","tsu00020.xml","tsu00030.xml",
               "tsu00040.xml","tsu00051.xml","tsu00052.xml",
               "tsu00061.xml","tsu00061.xml","tsu00500.xml",
               "tsu000511.xml","tsu00620.xml",
               "tsu00680.xml","tsu01200.xml")
files_tsu <- system.file(paste("/extdata/", tsu_files, sep=""), package="KEGGgraph")
kgml_tsu <- sapply(files_tsu, parseKGML)
reaction_tsu <- sapply(kgml_tsu, KEGGpathway2reactionGraph)
pan_tsu <- mergeGraphs(reaction_tsu)
tsu_plot <- plot(pan_tsu,attrs=list(node=list(label="foo", fillcolor="lightgreen"),edge=list(color="black")))

################################### SHI #################################

shi_files <- c("shi00010.xml","shi00020.xml","shi00030.xml",
               "shi00040.xml","shi00051.xml","shi00052.xml",
               "shi00061.xml","shi00061.xml","shi00500.xml",
               "shi000511.xml","shi00620.xml",
               "shi00680.xml","shi01200.xml")

files_shi <- system.file(paste("/extdata/", shi_files, sep=""), package="KEGGgraph")

kgml_shi <- sapply(files_shi, parseKGML)
reaction_shi <- sapply(kgml_shi, KEGGpathway2reactionGraph)
pan_shi <- mergeGraphs(reaction_shi)
shi_plot <- plot(pan_shi,attrs=list(node=list(label="foo", fillcolor="lightgreen"),edge=list(color="black")))

#########################################################################


############## REACTION FOR ALL SIX SPECIES ##############################
files <- c("fsu00010.xml","fsu00020.xml","fsu00030.xml",
           "fsu00040.xml","fsu00051.xml","fsu00052.xml",
           "fsu00061.xml","fsu00500.xml","fsu000511.xml",
           "fsu00620.xml","fsu00680.xml","fsu01200.xml",
           "mst00010.xml","mst00020.xml","mst00030.xml",
           "mst00040.xml","mst00051.xml","mst00052.xml",
           "mst00500.xml","mst000511.xml","mst00620.xml",
           "mst00680.xml","mst01200.xml",
           "mru00010.xml","mru00020.xml","mru00030.xml",
           "mru00040.xml","mru00051.xml","mru00052.xml",
           "mru00061.xml","mru00500.xml","mru000511.xml",
           "mru00620.xml","mru00680.xml","mru01200.xml",
           "pru00010.xml","pru00020.xml","pru00030.xml",
           "pru00040.xml","pru00051.xml","pru00052.xml",
           "pru00061.xml","pru00500.xml","pru000511.xml",
           "pru00620.xml","pru00680.xml","pru01200.xml",
           "tsu00010.xml","tsu00020.xml","tsu00030.xml",
           "tsu00040.xml","tsu00051.xml","tsu00052.xml",
           "tsu00061.xml","tsu00500.xml","tsu000511.xml",
           "tsu00620.xml","tsu00680.xml","tsu01200.xml",
           "shi00010.xml","shi00020.xml","shi00030.xml",
           "shi00040.xml","shi00051.xml","shi00052.xml",
           "shi00061.xml","shi00500.xml","shi000511.xml",
           "shi00620.xml","shi00680.xml","shi01200.xml")

files <- system.file(paste("/extdata/", files, sep=""), package="KEGGgraph")
# parse KGML files 
kgml <- sapply(files, parseKGML)
# make individual reaction graphs
reaction <- sapply(kgml, KEGGpathway2reactionGraph)
# merge reaction graphs
pan <- mergeGraphs(reaction)
# plot merged reaction graph
all_reac <- plot(pan,attrs=list(node=list(label="foo", fillcolor="lightgreen", fontsize="200"),edge=list(color="black")))
