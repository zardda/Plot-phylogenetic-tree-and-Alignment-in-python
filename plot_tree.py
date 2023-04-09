from Bio import Phylo
import matplotlib.pyplot as plt
def draw_tree(file,ax):
    tree = Phylo.read(file,"newick")
    tree.as_phyloxml()
    tree.rooted = True
    for pos in ("left", "right", "top", "bottom"):
        ax.spines[pos].set_visible(False)
    ax.tick_params(left=False, labelleft=False,bottom=False,labelbottom=False)   
    Phylo.draw(tree,axes=ax)

