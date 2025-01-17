##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
## If using Espalier or this code, please cite:
##
##      Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.
##
############################

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import re
from Espalier.viz import balticmod as bt

def plot(tree_files,fig_name,**kwargs):
    
    """
        
        Plot trees in tree_files as a tanglegram
        
        Parameters:     
            tree_files (list): list of newick tree files containing local trees
            fig_name (str): figure name for output png file
              
        Optional keyword arguments:
            tree_labels (list) : Tree labels corresponding to trees files.
            numerical_taxa_names (boolean): set True if taxa names are numeric
            cmap (matplotlib colormap): colormap
            dispalce_scalar (float): scales displacement between neighboring trees based on max tree height 
            branch_width (float): branch thinkness/width in plotted trees
            xax_border (float): size of border on x-axis
            yax_border (float): size of border on y-axis
            tip_font_size (float): tip label font size
            
        New keyword arguments for Mark:
            tip_marker_size (float): Size of tip node markers (circles). Setting to zero removes tip markers
            tip_color_map (dict): Dictionary mapping each tip label to a color used for lines connecting nodes. Overrides colors set by cmap.
            scale_bars (boolean): adds scale bars showing tree heights if True
    
    """

    tree_labels = kwargs.get('tree_labels', None)
    numerical_taxa_names = kwargs.get('numerical_taxa_names', False) 

    # Can fiddle with these params to improve fig appearance    
    cmap = kwargs.get('cmap', mpl.cm.viridis) # mpl.cm.Spectral is nice too
    dispalce_scalar = kwargs.get('dispalce_scalar', 0.5)
    branch_width = kwargs.get('branch_width', 2) 
    xax_border = kwargs.get('xax_border', 0.1) # should be relative to cumulative-displace
    yax_border = kwargs.get('yax_border', 2) # should be relative to ySpan of trees
    tip_font_size = kwargs.get('tip_font_size', 12)
    
    # Set new args passed from keywords
    tip_marker_size =  kwargs.get('tip_marker_size', 30)
    tip_color_map = kwargs.get('tip_color_map', None)
    scale_bars = kwargs.get('scale_bars', False)
    
    # Load trees into tree dict
    trees={}
    segments = list(range(len(tree_files)))
    for idx,tr in enumerate(tree_files):
        output_tree = tr.replace('.tre','.nexus')
        convert2nexus(tr,output_tree,numerical_taxa_names)
        ll=bt.loadNexus(output_tree,absoluteTime=False)
        #ll.setAbsoluteTime(2020.0)
        trees[idx]=ll
    print('\nDone!')
    
    # Rescale tree heights so they are all equal
    tree_heights = []
    for t,tr in enumerate(trees.keys()): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        tree_heights.append(cur_tree.treeHeight)
    max_height_cap = max(tree_heights)
    for t,tr in enumerate(trees.keys()): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        for k in cur_tree.Objects: ## iterate over a flat list of branches
            k.length = k.length * (max_height_cap/cur_tree.treeHeight)
        cur_tree.traverse_tree() ## required to set heights
        cur_tree.treeStats() ## report stats about tree
    
    # Compute displaceAmount on the same scale as the tree heights
    displaceAmount= max_height_cap * dispalce_scalar
    
    # Extract tip positions
    tip_positions={x:{} for x in segments} ## remember the position of each tip in each tree
    for t,tr in enumerate(trees.keys()): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        for k in cur_tree.Objects:
            if k.branchType=='leaf':
                tip_positions[tr][k.name]=(k.height,k.y) ## remember (X, Y) position of tip
    
    for X in range(10): ## 10 untangling iterations
        print('iteration %d'%(X+1))
        for t,tr in enumerate(segments): ## iterate over each tree
            print(tr)
            ptr=segments[t-1] ## previous tree
            ntr=segments[t] ## next tree
            seg=trees[ptr] ## fetch appropriate tree
            nex_seg=trees[ntr]
            for k in sorted(nex_seg.Objects,key=lambda q:q.height): ## iterate over branches from most recent to oldest
                if k.branchType=='node': ## can only sort nodes
                    leaves=[[seg.tipMap[tip] for tip in w.leaves if tip in seg.tipMap] if w.branchType=='node' else [w.name] for w in k.children] ## descendent tips in current order
                    
                    for c in range(len(leaves)):
                        leaves[c]=sorted(leaves[c],key=lambda x:tip_positions[ntr][x][1] if x in tip_positions[ntr] else 0.0) ## sort leaves according to their positions in the next tree
                    
                    ys=[sorted([tip_positions[ntr][w][1] for w in cl if w in tip_positions[ntr]]) for cl in leaves] ## extract y positions of descendents
                    merge_ys=sum(ys,[]) ## flatten list of tip y coordinates
                    ypos=range(min(merge_ys),max(merge_ys)+1) ## get y positions of tips in current order
                    order={i:x for i,x in enumerate(leaves)} ## dict of tip order: tip name
                    
                    new_order=sorted(order.keys(),key=lambda x:-np.mean([(tip_positions[ptr][order[x][w]][1]-ypos[w]) for w in range(min([len(order[x]),len(ypos)])) if order[x][w] in tip_positions[ptr]])) ## get new order by sorting existing order based on y position differences
                    if new_order!=range(len(leaves)): ## if new order is not current order
                        k.children=[k.children[i] for i in new_order] ## assign new order of child branches
                        nex_seg.drawTree() ## update y positions
    
                        for w in nex_seg.Objects: ## iterate over objects in next tree
                            if w.branchType=='leaf':
                                tip_positions[ntr][w.name]=(w.height,w.y) ## remember new positions
                    
            if t==0: ## if first tree
                trees[segments[t]].drawTree() ## update positions
                lvs=sorted([w for w in trees[segments[t]].Objects if w.branchType=='leaf'],key=lambda x:x.y) ## get leaves in y position order
                
                norm=mpl.colors.Normalize(0,len(lvs))
                pos_colours={w.name:cmap(norm(w.y)) for w in lvs} ## assign colour
                
    
    # Plotting all trees
    fig,ax = plt.subplots(figsize=(48,10),facecolor='w')
    cumulative_displace=0 ## this tracks the "current" x position, so trees are plotted one after another
    tree_names = segments
    ref_tree = segments[0] # first tree
    
    for t,tr in enumerate(tree_names): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        
        x_attr=lambda k: k.height+cumulative_displace
        #x_attr=lambda k: (k.height*(max_height/cur_tree.treeHeight))+cumulative_displace
        
        b_func=lambda k: branch_width 
        s_func=lambda k: tip_marker_size # 30
        su_func=lambda k: 2 * tip_marker_size # 60
        if tip_color_map:
            #k_color = tip_color_map[k.name]
            ct_func=lambda k: tip_color_map[k.name]
        else:
            ct_func=lambda k: cmap(tip_positions[ref_tree][k.name][1]/float(cur_tree.ySpan))
        cu_func=lambda k: 'k'
        z_func=lambda k: 100
        zu_func=lambda k: 99
        
        # For tip naming
        text_func = lambda k: k.name.replace('_',' ')
        target_func = lambda k: k.is_leaf()
        position_func = lambda k: (k.height+cumulative_displace+0.2, k.y)
        
        def colour_func(node):
            #if traitName in node.traits:
            #    return 'indianred' if node.traits[traitName]=='V' else 'blue'
            #else:
                return 'k'
            
        cn_func=colour_func
        
        cur_tree.plotTree(ax,x_attr=x_attr,branchWidth=b_func,colour_function=cn_func)
        cur_tree.plotPoints(ax,x_attr=x_attr,size_function=s_func,colour_function=ct_func,zorder_function=z_func)
        cur_tree.plotPoints(ax,x_attr=x_attr,size_function=su_func,colour_function=cu_func,zorder_function=zu_func)
        
        # Add tip label if at last tree
        if t == len(tree_names) - 1: # last_tree
            cur_tree.addText(ax, text=text_func, position=position_func,fontsize=tip_font_size)
        
        for k in cur_tree.Objects: ## iterate over branches
            if isinstance(k,bt.leaf): ## if leaf...
                y=k.y
                
                if tip_color_map:
                    k_color = tip_color_map[k.name]
                else:
                    pos_in_first_tree=tip_positions[ref_tree][k.name][1] ## fetch y coordinate of same tip in the first tree
                    frac_pos=pos_in_first_tree/float(cur_tree.ySpan) ## normalize coordinate to be within interval [0.0,1.0]
                    k_color = cmap(frac_pos)
    
                if t!=len(tree_names)-1: ## as long as we're not at the last tree - connect tips with coloured lines
                    next_x,next_y=tip_positions[tree_names[t+1]][k.name] ## fetch coordinates of same tip in next tree
                    next_x+=cumulative_displace+cur_tree.treeHeight+displaceAmount ## adjust x coordinate by current displacement and future displacement
                    nextIncrement=cumulative_displace+cur_tree.treeHeight
                    ax.plot([x_attr(k),nextIncrement+0.05*displaceAmount,nextIncrement+0.95*displaceAmount,next_x],[y,y,next_y,next_y],lw=1,ls='-',color=k_color,zorder=0) ## connect current tip with same tip in the next tree
        
        if tree_labels:
            add_tree_label(ax,cur_tree,tree_labels[tr],cumulative_displace)
        
        if scale_bars:
            label_str = f"Height = {tree_heights[t]:.2f}"
            add_tree_scale(ax,cur_tree,label_str,cumulative_displace)
        
        cumulative_displace+=cur_tree.treeHeight+displaceAmount ## increment displacement by the height of the tree
    
    [ax.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']]
    
    ax.tick_params(axis='x',size=0)
    ax.tick_params(axis='y',size=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    ax.set_ylim(-yax_border,cur_tree.ySpan+yax_border) ## set y limits
    ax.set_xlim(-xax_border,cumulative_displace+xax_border)
    
    plt.savefig(fig_name, dpi=300)
    
def convert2nexus(in_tree,out_tree,numerical_taxa_names=False):
    
    """
        Convert newick file to nexus for tanglegram plotting
    """
    
    myTree=bt.loadNewick(in_tree, absoluteTime=False)
    myTree.traverse_tree() ## required to set heights
    myTree.treeStats() ## report stats about tree
    names = []
    for idx,k in enumerate(myTree.Objects): ## iterate over a flat list of branches
        if k.branchType=='leaf':
            curr_name = k.numName
            names.append(curr_name)
    
    date_str = '' #'_2020.00'
    
    # Write taxa names
    nex=open(out_tree,"w")
    nex.write("#NEXUS\n")
    nex.write("Begin taxa;\n")
    nex.write("\tDimensions ntax=" + str(len(names)) + ";\n")	
    nex.write("\t\tTaxlabels\n")
    for n in names:
        nex.write("\t\t\t" + n + date_str + "\n")    
    nex.write("\t\t\t;\n")
    nex.write("End;\n")		
    
    # Write translation 	
    nex.write("Begin trees;\n")	
    nex.write("\tTranslate\n")	
    for idx,n in enumerate(names):
        if numerical_taxa_names:
            nex.write("\t\t" + n + ' ' + n + date_str + "\n") # if taxa names are numbers
        else:
            nex.write("\t\t" + str(idx+1) + ' ' + n + date_str + "\n") # if taxa names are non-numerical strings
       
    nex.write(";\n")
    
    # Write tree
    with open(in_tree, 'r') as file:
        tree_str = file.read().replace('\n', '')
    if not numerical_taxa_names:
        for idx,n in enumerate(names):
            tree_str = re.sub(n, str(idx+1), tree_str) # if taxa names are non-numerical strings    
    nex.write("tree TREE1 = " + tree_str + "\n")
    nex.write("End;\n")

def add_tree_label(ax,tree,label_str,cumulative_displace):
    
    """
        Add a label to tree
    """
    
    curr_min_x = np.Inf
    curr_max_x = -np.Inf
    curr_min_y = np.Inf
    curr_max_y = -np.Inf
    for k in tree.Objects:
        if k.x > curr_max_x:
            curr_max_x = k.x
        if k.x < curr_min_x:
            curr_min_x = k.x
        if k.y > curr_max_y:
            curr_max_y = k.y
        if k.y < curr_min_y:
            curr_min_y = k.y
    x_text_pos = cumulative_displace + (curr_max_x - curr_min_x) / 2
    y_text_pos = curr_max_y + (curr_max_y - curr_min_y) * 0.05
    ax.text(x_text_pos,y_text_pos,label_str,horizontalalignment='center',fontsize=12)
    
def add_tree_scale(ax,tree,label_str,cumulative_displace):
    
    """
        Add scale bar to show true (unscaled) tree height
    """
    
    curr_min_x = np.Inf
    curr_max_x = -np.Inf
    curr_min_y = np.Inf
    curr_max_y = -np.Inf
    for k in tree.Objects:
        if k.x > curr_max_x:
            curr_max_x = k.x
        if k.x < curr_min_x:
            curr_min_x = k.x
        if k.y > curr_max_y:
            curr_max_y = k.y
        if k.y < curr_min_y:
            curr_min_y = k.y
    
    x1 = cumulative_displace + curr_min_x
    x2 = cumulative_displace + curr_max_x
    y1 = curr_min_y - (curr_max_y - curr_min_y) * 0.05
    y2 = y1
    ax.plot([x1, x2], [y1, y2], 'ko-', linewidth=1.5, markersize=4)
    
    x_text_pos = cumulative_displace + (curr_max_x - curr_min_x) / 2
    y_text_pos = curr_min_y - (curr_max_y - curr_min_y) * 0.1
    ax.text(x_text_pos,y_text_pos,label_str,horizontalalignment='center',fontsize=12)


if  __name__ == '__main__':
    
    "Test with multiple trees"
    segments = 20
    path = "/Users/mfarman/MCCT5/"
    ML_tree_files = [path + "MCCT" + str(i) + ".tre" for i in range(segments)]
    tanglegram_fig_name = 'tanglegram-test.png' 
    
    "List of strings/labels for each individual tree"
    tree_labels = ['Tree ' + str(s) for s in range(1,segments+1)]
    
    """
        Map each tip label to a color used for connecting tip nodes across trees
        Can be any named color in matplotlib: https://matplotlib.org/stable/gallery/color/named_colors.html
    """
    """
    tip2color = {'1':'indianred',
                 '2':'indianred',
                 '3':'indianred',
                 '4':'indianred',
                 '5':'indianred',
                 '6':'blue',
                 '7':'blue',
                 '8':'blue',
                 '9':'blue',
                 '10':'blue'}
    """
    tip2color = {'FY25':'blue',
                 'FY72':'blue',
                 'TO43':'blue',
                 'NL35':'blue',
                 'GR001':'blue',
                 'GR03':'blue',
                 'Fg197':'deeppink',
                 'Fg199':'deeppink',
                 'Fg205':'deeppink',
                 'Fg207':'deeppink',
                 'Fg233423':'blue',
                 'Fg241165':'deeppink',
                 'Fg15':'deeppink',
                 'Fg7':'deeppink',
                 'CML3064':'deeppink',
                 'CML3065':'deeppink',
                 'CML3066':'deeppink',
                 'CML3067':'deeppink',
                 'CML3068':'deeppink',
                 'CML3069':'deeppink',
                 'CML3070':'deeppink',
                 'CML3071':'deeppink',
                 'CML3402':'deeppink',
                 'CML3403':'deeppink',
                 'CML3404':'deeppink',
                 'CML3405':'deeppink',
                 'CML3406':'deeppink',
                 'CML3407':'deeppink',
                 'CML3409':'deeppink',
                 'CS3005':'blue',
                 'DAOM180378':'blue',
                 'Fg233423':'blue',
                 'Fg241165':'deeppink',
                 'Gz3639':'blue',
                 'ITEM124':'deeppink',
                 'NRRL28336':'deeppink',
                 'F335':'deeppink',
                 '#F329':'deeppink',
                 'F332':'deeppink',
                 '66042':'lime',
                 '37425':'deeppink',
                 '37401':'deeppink',
                 '#46422':'deeppink',
                 '37525':'deeppink',
                 '47588':'deeppink',
                 '46434':'deeppink',
                 '52195':'deeppink',
                 '52008':'deeppink',
                 '#52847':'deeppink',
                 '#52429':'deeppink',
                 '38380':'indianred',
                 '36905':'indianred',
                 '#38405':'indianred',
                 '47605':'blue',
                 '45373':'blue',
                 '45156':'blue',
                 '#44211':'lime',
                 '44078':'blue',
                 '#44070':'blue',
                 '43884':'lime',
                 '43161':'lime',
                 '53173':'blue',
                 '47659':'lime',
                 '37408':'blue',
                 '37450':'blue',
                 '37509':'blue',
                 '47571':'blue',
                 '52005':'blue',
                 '52129':'blue',
                 '52512':'blue',
                 '52840':'blue',
                 '52850':'blue',
                 '52955':'blue',
                 'F333':'blue',
                 'F336':'blue',
                 'F337':'blue',
                 'F338':'blue',
                 'F325':'blue',
                 'F327':'blue',
                 'F330':'blue',
                 'F331':'blue',
                 'F341':'deeppink',
                 'F343':'deeppink',
                 'F342':'blue',
                 'F344':'blue',
                 '66030':'lime',
                 '66031':'lime',
                 '66040':'lime',
                 '66041':'blue',
                 'F326':'deeppink',
                 '66043':'lime',
                 '66044':'lime',
                 '66045':'lime',
                 '66047':'lime',
                 '66049':'blue',
                 'F328':'lime',
                 '53046':'deeppink',
                 '53184':'deeppink',
                 'F334':'deeppink',
                 'KY410':'deeppink'}
    
    cmap=mpl.cm.Spectral
    plot(ML_tree_files, tanglegram_fig_name,
         tree_labels=tree_labels,
         numerical_taxa_names=False,
         cmap=cmap,
         tip_font_size=10,
         tip_marker_size=30,
         tip_color_map=tip2color,
         scale_bars=True)

