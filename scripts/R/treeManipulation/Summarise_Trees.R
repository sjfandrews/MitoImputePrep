library(phangorn)
library(compiler)

phylo.str = function(x) {
  internal.nodes = x$Nnode
  terminal.nodes = length(x$tip.label)
  total.nodes = internal.nodes + terminal.nodes
  star = terminal.nodes + 1
  bif = 2*terminal.nodes - 2
  str.score = (total.nodes - star) / (bif - star)
  return(str.score)
}
phylo.str = cmpfun(phylo.str)

message("=============================")
message("======== USAGE HINTS ========")
message("ARGUEMENT 1:   TREE FILE 1")
message("ARGUEMENT 2:   TREE FILE 2")
message("ARGUEMENT 3:   SHUFFLE FILE 1")
message("ARGUEMENT 4:   SHUFFLE FILE 2")
message("=============================")

args = commandArgs(trailingOnly = TRUE)
tree_file1 = args[1] # Tree inferred from the original/pre-shuffle alignment
tree_file2 = args[2] # Tree inferred from the post-shuffle alignment
shuffle_file1 = args [3] # Shuffle summary file for the original/pre-shuffle alignment
shuffle_file2 = args [4] # Shuffle summary file for the post-shuffle alignment

out_file = paste0(sub(".treefile", "", tree_file1), "_treeSummary.csv")

#tree_file1 = "~/Desktop/SANDBOX/mtVertebrates/Ovis_aries.9940.MT.muscle.refined_MPboot_branchLengths.treefile"
#tree_file2 = "~/Desktop/SANDBOX/mtVertebrates/Ovis_aries.9940.MT.muscle.refined_SHUFFLE.iterative.1_MPboot_branchLengths.treefile"
#shuffle_file1 = "~/Desktop/SANDBOX/mtVertebrates/Ovis_aries.9940.MT.muscle.refined_SHUFFLE.sites.csv"
#shuffle_file2 = "~/Desktop/SANDBOX/mtVertebrates/Ovis_aries.9940.MT.muscle.refined_SHUFFLE.iterative.1.csv"

tree1 = di2multi(read.tree(file = tree_file1), tol = 1e-20)
tree2 = di2multi(read.tree(file = tree_file2), tol = 1e-20)

shuffle1 = read.csv(shuffle_file1, header = T)
shuffle2 = read.csv(shuffle_file2, header = T)

tree_names = c(basename(tree_file1), basename(tree_file2))
topology_structure = c(phylo.str(tree1), phylo.str(tree2))
n_node = c(tree1$Nnode, tree2$Nnode)
mean_bootstrap = c(mean(as.numeric(tree1$node.label), na.rm = T), mean(as.numeric(tree2$node.label), na.rm = T))
mean_branchLength = c(mean(tree1$edge.length, na.rm = T), mean(tree2$edge.length, na.rm = T))
Co_scores = c(mean(shuffle1$Co, na.rm = T), mean(shuffle2$Co, na.rm = T))


outDF = data.frame(tree_names, topology_structure, n_node, mean_bootstrap, mean_branchLength, Co_scores)

message(paste0("WRITING SUMMARY FILE TO: ", out_file))
write.csv(x = outDF, file = out_file, quote = F, row.names = F)
