library(ggtree)
library(ggplot2)
setwd("/home/ubuntu/ICMC/steno-babiker/bactopia_gtdbtk/stenotrophomonas-maltophilia/panOut/bactopia-tools/pangenome/pangenome/iqtree")

tree<-ape::read.tree('core-genome.aln.treefile')

ggtree::ggtree(tree)+
  geom_tiplab()+
  geom_treescale()+
  theme(text = element_text(size = 13))
  
  

ggsave('simpleTree.jpeg', width = 23, height = 16, units = 'in', dpi=300)



ggtree::ggtree(tree, layout="circular")+
  geom_tiplab()+
  geom_treescale()+
  theme(text = element_text(size = 12))



ggsave('simpleTreeCirc.jpeg', width = 14, height = 14, units = 'in', dpi=300)