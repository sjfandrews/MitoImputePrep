library(phangorn)
library(phytools)

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

tree1_file = args[1] # Tree from original (pre-shuffle) alignment

tree1 = read.tree(tree1_file)
tree1 = midpoint(tree1)
tree1 = force.ultrametric(tree1)

write.tree(tree1, sub(x = tree1_file, pattern = ".treefile", "_ultrametric.treefile"))
