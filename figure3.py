####################################################
####------------------figure 3------------------####
####################################################

###code in Python###
###--MuSCs velocity analysis--###
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
qMsample=pd.read_table('qM/qMsample.txt',header=None)
qM=qMsample.iloc[:,0].values.tolist()
loompy.combine(qM,"qM/qM-merged.loom",key="Accession")

os.chdir("/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/Trajectory/RNAvelocity/qM")
veloqM = scv.read("qM-merged.loom", cache=True)
dataqM = sc.read_h5ad('/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/qM/MUSC/qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad')
#unify the cell barcode and merge
veloqM.obs = veloqM.obs.rename(index = lambda x: re.sub('x','',re.sub('.....:', '',x))) #231080B_32_qM_1_H8LR1:AAACCCAGTGTCCAATx
dataqM.obs = dataqM.obs.rename(index = lambda x: re.sub('-.*-.*','',x)) #231080B_32_qM_1_AAACCCAAGGCTTAAA-1-21
dataqM = scv.utils.merge(dataqM,veloqM)
dataqM.write('qM_MUSC_merged_loom.h5ad')

scv.pp.filter_and_normalize(dataqM)
scv.pp.moments(dataqM, n_pcs=30, n_neighbors=30) #Computes moments for velocity estimation.
scv.tl.recover_dynamics(dataqM,n_jobs=30)
scv.tl.velocity(dataqM, mode='dynamical')
#Visualization
color_dict = {"Activated-1": "#E4A031", "Activated-2" : "#F6E2C1", "Activated-3" : "#D75B4E", 
              "Transitional" : "#A1D4A2", "Myoblast-1" : "#7C4D77", "Myoblast-2" : "#B6B3D6","Myocyte" : "#26445E",
              "Quiescent-1" : "#DBE4FB","Quiescent-2" : "#1f77b4"}
mpl.rcParams['figure.figsize']=(8,12)
scv.tl.velocity_graph(dataqM)
scv.pl.velocity_embedding_stream(dataqM, basis='X_umap',color = 'labels',save='figure3a.pdf')
scv.tl.latent_time(dataqM)
scv.pl.scatter(dataqM, color='latent_time', color_map='gnuplot', size=80,save='figureS3a.pdf')
top_genes = dataqM.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(dataqM, var_names=top_genes, sortby='latent_time', col_color='labels', n_convolve=100,save='figureS3b.pdf')

