library(phangorn)
library(phytools)
library(ggtree)
library(ggplot2)
require(scales)
require(gridExtra)

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

tree1_file = args[1] # Tree from original (pre-shuffle) alignment
#tree2_file = args[2] # Tree from cleaned (post-shuffle) alignment
#population = args[3] # Population
#tree_summary_file = args[4]

#tree1_file = "/Volumes/TimMcInerney/1000_Genomes_Test/raijin/analyses/shuffle/chrMT/ACB/1kGP_P3_ACB_chrMT_MPboot_branchLengths.treefile"
#tree2_file = "/Volumes/TimMcInerney/1000_Genomes_Test/raijin/analyses/shuffle/chrMT/ACB/1kGP_P3_ACB_chrMT_SHUFFLE.iterative.1_MPboot_branchLengths.treefile"
#tree_summary_file = "/Volumes/TimMcInerney/1000_Genomes_Test/raijin/analyses/shuffle/chrMT/ACB/1kGP_P3_ACB_chrMT_MPboot_branchLengths_treeSummary.csv"

#species = basename(tree1_file)
#species = gsub("\\..*","", species)
#species = sub("_", " ", species)

tree1 = read.tree(tree1_file)
#tree2 = read.tree(tree2_file)

tree1 = midpoint(tree1)
#tree2 = midpoint(tree2)

tree1 = force.ultrametric(tree1)
#tree2 = force.ultrametric(tree2)

#RF.dist(tree1, tree2, check.labels = T)

out_plot = sub(".treefile", "_cladogram.png", tree1_file)

tree1_plot = ggtree(tree1, ladderize = T) +
  geom_tiplab(size = rel(2)) +
  geom_label(aes(x=branch, label = edge.length))

plot_multiplier = 2
wd = 297 * plot_multiplier
ht = 210 * plot_multiplier
ggsave(filename = out_plot, plot = tree1_plot, width = wd, height = ht, units = "mm", dpi = 300, limitsize = F)
