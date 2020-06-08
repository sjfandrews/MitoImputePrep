library(ggtree)
library(ggplot2)
require(scales)
require(gridExtra)

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

tree1_file = args[1] # Tree from original (pre-shuffle) alignment
tree2_file = args[2] # Tree from cleaned (post-shuffle) alignment
#population = args[3] # Population
#tree_summary_file = args[4]

#tree1_file = "/Volumes/TimMcInerney/1000_Genomes_Test/raijin/analyses/shuffle/chrMT/ACB/1kGP_P3_ACB_chrMT_MPboot_branchLengths.treefile"
#tree2_file = "/Volumes/TimMcInerney/1000_Genomes_Test/raijin/analyses/shuffle/chrMT/ACB/1kGP_P3_ACB_chrMT_SHUFFLE.iterative.1_MPboot_branchLengths.treefile"
#tree_summary_file = "/Volumes/TimMcInerney/1000_Genomes_Test/raijin/analyses/shuffle/chrMT/ACB/1kGP_P3_ACB_chrMT_MPboot_branchLengths_treeSummary.csv"

species = basename(tree1_file)
species = gsub("\\..*","", species)
species = sub("_", " ", species)

tree1 = read.tree(tree1_file)
tree2 = read.tree(tree2_file)

#tree_summary = read.csv(tree_summary_file, header = T)

out_plot = sub(".treefile", ".png", tree1_file)

tree1_plot = ggtree(tree1, layout="equal_angle") +
  labs(title = paste0("Unrooted Maximum Parsimony tree from original alignment\n",
                      species))
#tree1_plot

tree2_plot = ggtree(tree2, layout="equal_angle") +
  labs(title = paste0("Unrooted Maximum Parsimony tree from post-shuffle alignment\n",
                      species))
#tree2_plot

both_trees = grid.arrange(arrangeGrob(tree1_plot, tree2_plot,
                                      ncol = 2))


message(message(paste0("SAVING EBSP PLOTS TO ", out_plot)))
plot_multiplier = 2
wd = 25.4 * plot_multiplier
ht = 11.95 * plot_multiplier
ggsave(filename = out_plot, plot = both_trees, width = wd, height = ht, units = "cm", dpi = 300, limitsize = F)
