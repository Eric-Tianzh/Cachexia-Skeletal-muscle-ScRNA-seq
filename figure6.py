####################################################
####------------------figure 6------------------####
####################################################

###code in Python###
###--scRNA-seq data analysis--###
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 
from matplotlib.colors import LinearSegmentedColormap

sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
adata_all=sc.read_h5ad('gM.after_batch.marker.h5ad')
color_dict = {'B': "#1f77b4",'EC-1': "#38917E",'EC-2': "#66c2a4",'EC-3': "#56b4b0",'EC-4': "#238b45",
                'EC-5': "#a8ddb5",'FAPs-1': "#dd7c4f",'FAPs-2': "#CC5200",'FAPs-3': "#FF9933",'FAPs-4': "#FFD699",
                'Glial_cell': "#af8dc3",'Mono/Macro': "#c51b7d",'MUSC': "#d0dd97",'Myonuclei': "#4D261A",'Neutrophil': "#17becf",
                'Schwann cell': "#963B79",'Smooth M-1': "#f1a0c4",'Smooth M-2': "#fbb4b9",'Smooth M-3': "#F08080",'T': "#FFD700",
                'Tenocytes-1': "#7895C1",'Tenocytes-2': "#2171b5",'Tenocytes-3': "#4292c6",'unknown': "#808080"}
sc.pl.umap(adata_all,color=["labels"],legend_loc="on data",legend_fontsize=6,frameon=False,palette=color_dict,save='figure6a.pdf')

data=sc.read_h5ad('gM-MUSC.after_batch.bbknn.marker.h5ad')
color_dict = {"Activated-1": "#E4A031", "Activated-2" : "#F6E2C1", "Activated-3" : "#D75B4E", "Activated-4" : "#DBE4FB",
              "Transitional" : "#A1D4A2", "Myoblast-1" : "#7C4D77", "Myoblast-2" : "#B6B3D6","Myocyte" : "#26445E",
              "Quiescent" : "#1f77b4"}

sc.pl.umap(data,color=["labels"],legend_loc="on data",legend_fontsize=6,frameon=False,palette=color_dict,save='figure6b.pdf')


import loompy
import os
import re
import pandas as pd
import numpy as np #<2.0
import scanpy as sc
import velocyto as vlc
import anndata
from matplotlib import pyplot as plt
import matplotlib
class mplDeprecation(UserWarning): pass
matplotlib.cbook.mplDeprecation = mplDeprecation
import scvelo as scv
sc.logging.print_header()
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 

###merge samples###
os.chdir("/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/Trajectory/RNAvelocity")
gMsample=pd.read_table('gM/gMsample.txt',header=None)
gM=gMsample.iloc[:,0].values.tolist()
loompy.combine(gM,"gM/gM-merged.loom",key="Accession")

os.chdir("/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/Trajectory/RNAvelocity/gM")
velogM = scv.read("gM-merged.loom", cache=True)
datagM = sc.read_h5ad('/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/gM/MUSC/gM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad')
#unify the cell barcode and merge
velogM.obs = velogM.obs.rename(index = lambda x: re.sub('x','',re.sub('.....:', '',x)))
datagM.obs = datagM.obs.rename(index = lambda x: re.sub('-.*-.*','',x))
datagM = scv.utils.merge(datagM,velogM)
datagM.write('gM_MUSC_merged_loom.h5ad')

scv.pp.filter_and_normalize(datagM)
scv.pp.moments(datagM, n_pcs=30, n_neighbors=30) #Computes moments for velocity estimation.
scv.tl.recover_dynamics(datagM,n_jobs=30)
scv.tl.velocity(datagM, mode='dynamical')
#Visualization
mpl.rcParams['figure.figsize']=(8,12)
scv.tl.velocity_graph(datagM)
scv.pl.velocity_embedding_stream(datagM, basis='X_umap',color = 'labels',save='figure6d.pdf')