####################################################
####------------------figure 2------------------####
####################################################

###code in Python###
###--MuSCs scRNA-seq data analysis--###
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib.colors import LinearSegmentedColormap

import math
import itertools
import warnings
from collections import Counter
from SCCAF import *

sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

data=sc.read_h5ad('qM.after_batch.cluster.h5ad')
label_names = pd.read_csv('cell-type1.csv')
label_names['ClusterID'] =label_names['ClusterID'] .astype(str)
labels = label_names.set_index(label_names.columns[0]).loc[data.obs['leiden_res0.5']].iloc[:, 0].values
data.obs['labels'] = labels
cluster_obs=data.obs['labels'][data.obs['labels']=='MuSCs'] #17350
cluster_cell=cluster_obs.index.values.tolist()
adata_all=sc.read_h5ad('qM.after_batch.rm_blacklist_counts.h5ad')
adata_all=adata_all[cluster_cell,:]

##after QC
adata_all.obs.drop(['n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_hb', 'pct_counts_hb'],axis=1,inplace=True)
adata_all.var.drop(['n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'],axis=1,inplace=True)
sc.pp.calculate_qc_metrics(adata_all, qc_vars=['mt','hb'], percent_top=None, log1p=False, inplace=True)
adata_all.write('qM-MUSC.after_batch.rm_blacklist_counts.h5ad')

##normalize
sc.pp.normalize_total(adata_all, target_sum=1e4) #normalize counts per cell
sc.pp.log1p(adata_all)
adata_all.write('qM-MUSC.after_batch.rm_blacklist_logcounts.h5ad')

##Identify highly-variable genes
sc.pp.highly_variable_genes(adata_all, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_all,save="qM-MUSC.highly_variable_genes.pdf")

adata_all.raw = adata_all
adata_all.raw.to_adata().write('qM-MUSC.after_batch_raw.h5ad')

###filtering highly-variable genes
adata_all = adata_all[:, adata_all.var.highly_variable] #15731 × 2767
sc.pp.regress_out(adata_all, ['total_counts', 'pct_counts_mt'])#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
sc.pp.scale(adata_all, max_value=10)

##Principal component analysis
sc.tl.pca(adata_all, svd_solver='arpack',n_comps=100)
sc.pl.pca_variance_ratio(adata_all, log=True, n_pcs=100,save="qM-MUSC.pca_variance_ratio.pdf")
adata_all.write('qM-MUSC.after_batch.pca.h5ad')

##Computing the neighborhood graph
adata_all = sc.read_h5ad("qM-MUSC.after_batch.pca.h5ad")
sc.pp.neighbors(adata_all, n_pcs=50)

##BBKNN batch correction
adata_all = sc.read_h5ad("qM-MUSC.after_batch.pca.h5ad")
sce.pp.bbknn(adata_all,batch_key='donor',n_pcs=50)
sc.tl.umap(adata_all)
for i in np.arange(0.1, 2.1, 0.1):
    i_name='leiden_res'+str(i)[:3]
    i_name_save='qM-MUSC.bbknn.umap.leiden_res'+str(i)[:3]+'.pdf'
    print(i_name)
    print(i_name_save)
    sc.tl.leiden(adata_all,resolution=i,key_added=i_name)
    sc.pl.umap(adata_all, color=i_name, save=i_name_save,legend_fontsize='medium')
adata_all.write('qM-MUSC.after_batch.bbknn.cluster.group.h5ad')

##Finding marker genes
adata_all.uns['log1p']["base"] = None
for i in np.arange(0.1, 2.1, 0.1):
    i_name='leiden_res'+str(i)[:3]
    i_name_add='leiden_res'+str(i)[:3]+'_marker'
    i_name_marker='qM-MUSC.bbknn.rank_genes_groups.leiden_res'+str(i)[:3]+'_marker.pdf'
    dedf_out=i_name_add+"bbknn.csv"
    print(i_name_add)
    sc.tl.rank_genes_groups(adata_all, i_name, method='wilcoxon',key_added=i_name_add, pts=True)
    sc.pl.rank_genes_groups(adata_all, n_genes=25, sharey=False,key=i_name_add,save=i_name_marker)
    dedf = sc.get.rank_genes_groups_df(adata_all,key=i_name_add,group=None)
    dedf.to_csv(dedf_out)
adata_all.write('../qM-MUSC.after_batch.bbknn.marker.h5ad')

##after BBKNN, delete cluster which cell number less than 200,integrate cluster with similar marker genes
data=sc.read_h5ad('qM-MUSC.after_batch.bbknn.marker.h5ad')
cluster_obs = data.obs['leiden_res1.0'][data.obs['leiden_res1.0'].isin(['10','12', '13', '14'])] #294
cluster_cell=cluster_obs.index.values.tolist()
data = data[~data.obs.index.isin(cluster_cell)] 

label_names = pd.read_csv('cell-type-rename_res1.0.csv')
label_names['ClusterID'] =label_names['ClusterID'] .astype(str)
labels = label_names.set_index(label_names.columns[0]).loc[data.obs['leiden_res1.0']].iloc[:, 0].values
data.obs['labels'] = labels
data.obs['labels'] = data.obs['labels'].astype('category')

group=["Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3","Transitional","Myoblast-1","Myoblast-2","Myocyte"]
##calculate the marker of new cluster labels
data.uns['log1p']["base"] = None
##Finding marker genes
sc.tl.rank_genes_groups(data, 'labels', method='wilcoxon',key_added='labels_marker', pts=True)
dedf = sc.get.rank_genes_groups_df(data,key='labels_marker',group=None)
dedf.to_csv('qm_musc_res1.0_labels_marker.csv')
data.write('qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad')

color_dict = {"Activated-1": "#E4A031", "Activated-2" : "#F6E2C1", "Activated-3" : "#D75B4E", 
              "Transitional" : "#A1D4A2", "Myoblast-1" : "#7C4D77", "Myoblast-2" : "#B6B3D6","Myocyte" : "#26445E",
              "Quiescent-1" : "#DBE4FB","Quiescent-2" : "#1f77b4"}

###--UMAP of MuSCs clusters--###
sc.pl.umap(data,color=["labels"],legend_loc="on data",legend_fontsize=6,frameon=False,size=8,palette=color_dict,save='figure2a.pdf') 

###--Dotplot of marker genes--###
group=["Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3","Transitional","Myoblast-1","Myoblast-2","Myocyte"]  
sc.pl.dotplot(data,['Pax7','Myf5','Ttc5','Pde4c','Nr4a1','Ier5l','Igfbp5','Pvr','Zc3h12a','Suco',
                    'Sphk1','Cbx4','Myod1','Srxn1','Neu2','Gcnt2','Myog'], groupby="labels",cmap=LinearSegmentedColormap.from_list("custom_gradient", ["#f0ead2","#A1D4A2","#1f77b4"], N=256),
                          standard_scale="var",figsize=(7,4),categories_order=group,
                          save='figure2b.pdf')

###--UMAP plot of characteristic genes in MuSCs--###
cmapu=sns.color_palette("vlag", as_cmap=True)
sc.pl.umap(data,color=['Pax7','Myod1','Myf5','Myog'],legend_loc="on data",color_map=cmapu,legend_fontsize='x-small',size=25,frameon=False,save='figure2c.pdf')

###--UMAP distribution of donor,group--###
sc.pl.umap(adata_all, color=['donor'], save='figureS2a.pdf',legend_fontsize='medium')
dataqM_LLC=data[data.obs['paired_group_name'] == 'LLC',]
dataqM_PF=data[data.obs['paired_group_name'] == 'PF',]
sc.pl.umap(dataqM_LLC,color=["labels"],legend_fontsize=6,frameon=False,size=20,palette=color_dict,save='figureS2c_left.pdf')
sc.pl.umap(dataqM_PF,color=["labels"],legend_fontsize=6,frameon=False,size=20,palette=color_dict,save='figureS2c_right.pdf')

###--ROC curve of MuSCs clusters--###
y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(data.X, data.obs['labels'], n=1000)
plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc,plot='roc',colors=color_dict,save='figure2d.pdf')


