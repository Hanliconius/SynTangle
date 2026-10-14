library(GENESPACE)
library(stringr)

wd <- "/netscratch/dep_mercier/grp_marques/marques/Rhync_tenuis_pangenome_project/GENESPACE/"

path2mcscanx = "/opt/share/software/scs/appStore/bookwormApps/synteny/MCScanX/vb1ca533/src"

path2orthofinder = "orthofinder"

path2diamond = "diamond"

genomes2run <- c("Rbreviuscula", "Raustrobrasiliensis",
                "Rhync_tenuis_ref.hap1", "Rhync_tenuis_ref.hap2",
                "Rhync_tenuis_6523A.hap1", "Rhync_tenuis_6523A.hap2",
                "Rhync_tenuis_6524A.hap1", "Rhync_tenuis_6524A.hap2",
                "Rhync_tenuis_6344A.JGV.hap1", "Rhync_tenuis_6344A.JGV.hap2",
                "Rhync_tenuis_6344B.PECP2.hap1", "Rhync_tenuis_6344B.PECP2.hap2",
                "Rhync_tenuis_6365A.035_5.hap1", "Rhync_tenuis_6365A.035_5.hap2",
                "Rhync_tenuis_6344C.PECP3.hap1", "Rhync_tenuis_6344C.PECP3.hap2",
                "Rhync_tenuis_6228D.036.hap1", "Rhync_tenuis_6228D.036.hap2",
                "Rhync_tenuis_6365B.036_7.hap1", "Rhync_tenuis_6365B.036_7.hap2")

gpar <- init_genespace(
  wd = wd,
  genomeIDs = genomes2run,
  ploidy = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder,
  path2diamond = path2diamond)


#customize plot
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

#customPal <- colorRampPalette(c("#D61F27","#F9AC60","#FCF7BF","#ADD8E7","#307BB5"))
  
customPal <- colorRampPalette(c("#F9AC60", "#307BB5", "#D61F27", "#ADD8E7", "#FCF7BF"))


pdf(paste(wd, "genespace_all_acc.pdf", sep="/"),
    family="Helvetica", height=11.7, width=8.3)
par(mai = c(0.4, 1, 0.1, 0.7)); # margin: bottom, left, top, right

ripd <- plot_riparian(
  gsParam = gpar,
  refGenome = "Rbreviuscula",
  genomeIDs = genomes2run,
  labelTheseGenomes = c("R. breviuscula", "R. austro-brasiliensis",
                "REF h1", "REF h2",
                "JGV 016 h1", "JGV 016 h2",
                "JGV 017 h1", "JGV 017 h2",
                "JGV h1","JGV h2",
                "PECP2 h1","PECP2 h2",
                "035-6 h1","035-6 h2",
                "PECP3 h1","PECP3 h2",
                "036 h1","036 h2",
                "036-7 h1","036-7 h2"),
  palette = customPal,
  addThemes = ggthemes,
  braidAlpha = .75,
  chrExpand = 0.75,
  chrLabFontSize = 8,
  chrBorderLwd = 0.3,
  chrBorderCol = "black",
  useOrder = FALSE,
  useRegions = TRUE,
  invertTheseChrs = data.frame(genome=c("Rbreviuscula"), chr=c("Chr3_h1")),
  #chrLabFun = function(x) str_replace_all(x, chr_dict_char),
  customRefChrOrder=c("Chr2_h1", "Chr5_h1", "Chr1_h1", "Chr4_h1", "Chr3_h1"))

dev.off()
