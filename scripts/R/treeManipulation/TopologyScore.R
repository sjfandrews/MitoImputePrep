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

args = commandArgs(trailingOnly = TRUE)
tree_file1 = args[1]
tree_file2 = args[2]

out_file = paste0(sub(".treefile", "", basename(tree_file1)), "_topologyScore.csv")

#tree_file1 = "~/Desktop/SANDBOX/mtVertebrates/Ovis_aries.9940.MT.muscle.refined_MPboot_branchLengths.treefile"
#tree_file2 = "~/Desktop/SANDBOX/mtVertebrates/Ovis_aries.9940.MT.muscle.refined_SHUFFLE.iterative.1_MPboot_branchLengths.treefile"

tree1 = read.tree(file = tree_file1)
tree2 = read.tree(file = tree_file2)

tree_names = c(basename(tree_file1), basename(tree_file2))
topology_structure = c(phylo.str(tree1), phylo.str(tree2))

outDF = data.frame(tree_names, topology_structure)

write.csv(x = outDF, file = out_file, quote = F, row.names = F)