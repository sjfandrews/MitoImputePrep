library(phangorn)

args <- commandArgs(trailingOnly = TRUE)
fasta_file = args[1]
tree_file = args[2]

message("===========================")
message("======= USAGE HINTS =======")
message("ARGUEMENT 1:   FASTA FILE")
message("ARGUEMENT 2:   TREE FILE")
message("===========================")
#fasta_file = "~/Desktop/SANDBOX/mtvertebrates/Ursus_arctos.9644.MT.muscle.refined.fa"
#tree_file = "~/Desktop/SANDBOX/mtvertebrates/Ursus_arctos.9644.MT.muscle.refined_MPboot.treefile"

out_tree1 = paste0(sub(".treefile", "", tree_file), "_branchLengths.treefile")
out_tree2 = paste0(sub(".treefile", "", tree_file), "_branchLengths_mprooted.treefile")

fasta = read.phyDat(fasta_file, format = "fasta", type = "DNA")
tree = read.tree(tree_file)

#fid = names(fasta)
#
#for (tip in 1:length(tree$tip.label)) {
#  tid = gsub("\\..*","", tree$tip.label[tip])
#  loc = match(tid, gsub("\\..*", "", fid))
#  fid[loc] = tree$tip.label[tip]
#}
#
#names(fasta) = fid

new.tree = acctran(tree, fasta)
new.tree = di2multi(new.tree)
#new.tree = midpoint(new.tree) # DONT MIDPOINT ROOT!

message(paste0("WRITING UNROOTED TREE TO:  ", out_tree1))
write.tree(phy = new.tree, file = out_tree1)

new.tree = midpoint(new.tree)

message(paste0("WRITING MID-POINT ROOTED TREE TO:  ", out_tree2))
write.tree(phy = new.tree, file = out_tree2)
