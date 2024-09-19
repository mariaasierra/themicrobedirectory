import ete3
from ete3 import NCBITaxa
ncbi = NCBITaxa() #download taxa from ncbi
taxonomy = [line.rstrip() for line in open("metasub_hosts_taxa.txt", "r")] #Include query taxonomy
name2taxid = ncbi.get_name_translator(taxonomy)
taxid=name2taxid.values()
taxnames=name2taxid.items()

taxid=sum(taxid, [])
tree = ncbi.get_topology(taxid, intermediate_nodes=True)
ncbi.annotate_tree(tree, taxid_attr="name")

for node in tree.traverse():
    node.name = node.sci_name
tree.write(outfile="tree_hosts_bac_metasub.txt", format=1) #Must change spaces with _ before importing in R