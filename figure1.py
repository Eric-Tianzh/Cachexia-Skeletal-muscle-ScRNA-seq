####################################################
####------------------figure 1------------------####
####################################################

###code in Python###
###--Quadriceps scRNA-seq data analysis--###
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

adata_all_count = sc.read_h5ad("../merge_all.after_batch.rm_blacklist_counts.h5ad")

cluster_obs=adata_all_count.obs['tissue'][adata_all_count.obs['tissue']=='qM']
cluster_cell=cluster_obs.index.values.tolist()
adata_all=adata_all_count[cluster_cell,:]

adata_all.obs.drop(['n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_hb', 'pct_counts_hb', 'n_genes', 'n_counts'],axis=1,inplace=True)
sc.pp.calculate_qc_metrics(adata_all, qc_vars=['mt','hb'], percent_top=None, log1p=False, inplace=True)
adata_all.write('qM.after_batch.rm_blacklist_counts.h5ad')

##after QC
adata_all=sc.read_h5ad('qM.after_batch.rm_blacklist_counts.h5ad')

##normalize
sc.pp.normalize_total(adata_all, target_sum=1e4)
sc.pp.log1p(adata_all)
adata_all.write('qM.after_batch.rm_blacklist_logcounts.h5ad')

##Identify highly-variable genes
sc.pp.highly_variable_genes(adata_all, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_all,save="qM.highly_variable_genes.pdf")

adata_all.raw = adata_all
adata_all.raw.to_adata().write('qM.after_batch_raw.h5ad')

##filtering highly-variable genes
adata_all = adata_all[:, adata_all.var.highly_variable]
sc.pp.regress_out(adata_all, ['total_counts', 'pct_counts_mt'])#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
sc.pp.scale(adata_all, max_value=10)

##Principal component analysis
sc.tl.pca(adata_all, svd_solver='arpack',n_comps=100)
sc.pl.pca_variance_ratio(adata_all, log=True, n_pcs=100,save="qM.pca_variance_ratio.pdf")
adata_all.write('qM.after_batch.pca.h5ad')

##Computing the neighborhood graph
adata_all = sc.read_h5ad("qM.after_batch.pca.h5ad")
sc.pp.neighbors(adata_all, n_pcs=75)

##Embedding the neighborhood graph
sc.tl.umap(adata_all)
adata_all.write('qM.after_batch.umap.h5ad')

##Clustering the neighborhood graph
for i in np.arange(0.1, 2.1, 0.1):
    i_name='leiden_res'+str(i)[:3]
    i_name_save='qM.umap.leiden_res'+str(i)[:3]+'.pdf'
    print(i_name)
    print(i_name_save)
    sc.tl.leiden(adata_all,resolution=i,key_added=i_name)
    sc.pl.umap(adata_all, color=i_name, save=i_name_save,legend_fontsize='medium')
adata_all.write('qM.after_batch.cluster.h5ad')
#output for Seurat format
data=sc.read_h5ad('qM.after_batch.cluster.h5ad')
data.write_csvs("qM_after_cluster",skip_data=False)#Write the anndata meta data into csv formate

##Finding marker genes
adata_all.uns['log1p']["base"] = None
for i in np.arange(0.1, 2.1, 0.1):
    i_name='leiden_res'+str(i)[:3]
    i_name_add='leiden_res'+str(i)[:3]+'_marker'
    i_name_marker='qM.rank_genes_groups.leiden_res'+str(i)[:3]+'_marker.pdf'
    dedf_out=i_name_add+".csv"
    print(i_name_add)
    sc.tl.rank_genes_groups(adata_all, i_name, method='wilcoxon',key_added=i_name_add, pts=True)
    sc.pl.rank_genes_groups(adata_all, n_genes=25, sharey=False,key=i_name_add,save=i_name_marker)
    dedf = sc.get.rank_genes_groups_df(adata_all,key=i_name_add,group=None)
    dedf.to_csv(dedf_out)
adata_all.write('qM.after_batch.marker.h5ad')

###--UMAP of Quadriceps clusters--###
color_dict = {'B': "#1f77b4",'EC-1': "#38917E",'EC-2': "#66c2a4",'EC-3': "#56b4b0",'EC-4': "#238b45",
                'FAPs-1': "#dd7c4f",'FAPs-2': "#CC5200",'FAPs-3': "#FF9933",'FAPs-4': "#FFD699",
                'Glial_cell': "#af8dc3",'Mono/Macro': "#c51b7d",'MuSCs': "#d0dd97",'Myonuclei': "#4D261A",'Neutrophil': "#17becf",
                'Schwann cell': "#963B79",'Smooth M-1': "#f1a0c4",'Smooth M-2': "#fbb4b9",'Smooth M-3': "#F08080",'T': "#FFD700",
                'Tenocytes-1': "#7895C1",'Tenocytes-2': "#2171b5",'Tenocytes-3': "#4292c6"}
sc.pl.umap(adata_all,color=["labels"],legend_loc="on data",legend_fontsize=6,frameon=False,palette=color_dict,save='figure1c.pdf')

###--Dotplot of marker genes--###
group = ['B', 'Mono/Macro', 'Neutrophil', 'T', 'EC-1', 'EC-2', 'EC-3', 'EC-4', 'FAPs-1', 'FAPs-2', 'FAPs-3', 'FAPs-4', 'Tenocytes-1', 'Tenocytes-2', 'Tenocytes-3', 'Glial_cell', 'Schwann cell', 'MuSCs', 'Myonuclei', 'Smooth M-1', 'Smooth M-2', 'Smooth M-3']
sc.pl.dotplot(adata_all, ['Cd79a', 'Cd19','Cd68','Cd14','Csf3r','S100a9','Cd3d','Nkg7',
                     'Kdr','Cdh5', 'Pdgfra','Dcn','Tnmd','Mkx','Ptn','Fabp7',
                     'Mpz','Mbp','Pax7', 'Myf5','Acta1','Myh4','Myh11','Myl9'], 
                     groupby="labels",figsize=(9,7),swap_axes=True,categories_order=group, standard_scale="var", save='figure1d.pdf')

###--UMAP plot of characteristic genes in MuSCs--###
sc.pl.umap(data,color=['Pax7','Myod1'],legend_loc="on data",legend_fontsize='x-small',frameon=False,save='figure1f.pdf')

###--UMAP distribution of donor, group and specific markers--###
sc.pl.umap(adata_all, color=['donor'], save='figureS1a.pdf',legend_fontsize='medium')
sc.pl.umap(adata_all, color=['paired_group_name'], save='figureS1b.pdf',legend_fontsize='medium')
cmapu = sns.light_palette("red", as_cmap=True)
sc.pl.umap(data,color=['Cd79a','Cd68','Csf3r','Cd3d','Cdh5', 'Pdgfra','Tnmd','Fabp7',
                     'Mpz','Pax7','Acta1','Myh11'],color_map=cmapu,frameon=False,colorbar_loc=None,save='figureS1c.pdf')
