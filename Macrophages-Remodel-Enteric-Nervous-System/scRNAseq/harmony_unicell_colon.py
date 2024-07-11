
import ucdeconvolve as ucd
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt

TOKEN = "uc_D8x68tmYlRmwuBJrwx1A40Q3uPfvvDKaVwHSpSXe8RD30Ldf"
ucd.api.authenticate(TOKEN)

sc.settings.verbosity = 1
sc.settings.figdir = './figures'
obj_save_path = '.'
mat_path = '.'

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400, facecolor = 'white', figsize=(6,6), format='png')

doublet_table = pd.read_csv("doublet_consensus.csv", header=None)
doublets = doublet_table[18]

root_path = "/Users/liuy39/Research/H026/scRNAseq/"
meta_df0 = pd.read_csv(root_path + "/metadata/scRNAseq_meta_summary.txt", sep="\t")
meta_df  = meta_df0.loc[meta_df0['Genoty'] == "wt" ]
samplenames =  meta_df['Animal_ID']
path2 = '/filtered_feature_bc_matrix.h5'
filenames = [root_path + "cellranger_out/"+ str(sample) + path2 for sample in samplenames]
adata_list = []

for i in range(len(filenames)):
  filename = filenames[i]
  sample   = samplenames[i]
  adata0 = sc.read_10x_h5(filename, gex_only = True)
  adata0.var_names_make_unique()
  adata0.obs['doublets'] = adata0.obs_names.isin(doublets)
  adata0 = adata0[adata0.obs.doublets == False, :]
  adata0.obs["Sex"] = meta_df.loc[meta_df.Animal_ID == sample, 'Sex'].values[0]
  adata0.obs["Genoty"] = meta_df.loc[meta_df.Animal_ID == sample, 'Genoty'].values[0]
  adata0.obs["batch"] = meta_df.loc[meta_df.Animal_ID == sample, 'batch'].values[0]
  adata0.obs["tissue"] = meta_df.loc[meta_df.Animal_ID == sample, 'tissue'].values[0]
  sc.pp.filter_cells(adata0, min_genes=200)
  sc.pp.filter_genes(adata0, min_cells=3)
  adata0.var['mt'] = adata0.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
  sc.pp.calculate_qc_metrics(adata0, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
  non_mito_rp_genes_list = [name for name in adata0.var_names if not name.startswith('mt-')] 
  non_mito_rp_genes_list2 = [name for name in non_mito_rp_genes_list if not name.startswith('rp') ] 
  adata0 = adata0[:, non_mito_rp_genes_list2]  
  adata0 = adata0[adata0.obs.n_genes_by_counts > 300 , :]
  adata0 = adata0[adata0.obs.n_genes_by_counts < 4487, :]
  adata0 = adata0[adata0.obs.pct_counts_mt < 20, :]
  adata0.layers["counts"] = adata0.X.copy()
  adata0.raw = adata0
  adata_list.append(adata0)

adata_dict = { k:v for (k,v) in zip(samplenames, adata_list)}
adatas = ad.concat(adata_dict, label="Animal_ID",  join="outer")
adatas.obs_names_make_unique()

cluster_counts = adatas.obs['Animal_ID'].value_counts() 
sc.pl.violin(adatas, 'n_genes', groupby='Animal_ID', stripplot=False, inner='box')  
sc.pl.violin(adatas, 'pct_counts_mt', groupby='Animal_ID', stripplot=False, inner='box') 

#keep = cluster_counts.index[cluster_counts >= 100] 
#adatas = adatas[adatas.obs['sample_name'].isin(keep)].copy()

sc.pl.violin(adatas, ['n_genes_by_counts'], save='_n_genes', jitter=0.4)
sc.pl.violin(adatas, ['total_counts'], save='_total_counts', jitter=0.4)
sc.pl.violin(adatas, ['pct_counts_mt'], save='_mito_pct', jitter=0.4)

adata_colon       = adatas[adatas.obs.tissue == "colon", :]
#adata_small_bowel = adatas[adatas.obs.tissue == "small_bowel", :]

adatas = adata_colon 

sc.pp.normalize_total(adatas, target_sum=1e4)
sc.pp.log1p(adatas) 
sc.pp.highly_variable_genes(adatas, min_mean=0.0125, max_mean=3, min_disp=0.25)
adatas.raw = adatas
adatas = adatas[:,adatas.var.highly_variable]
sc.pp.scale(adatas, max_value=10)
sc.tl.pca(adatas, svd_solver='arpack')
sc.pl.pca_variance_ratio(adatas, log=True, n_pcs=50, save='') # scanpy generates the filename automatically

##### without integration by harmony
#sc.pp.neighbors(adatas, n_neighbors=10, n_pcs=50) 
#sc.tl.umap(adatas)
#sc.tl.leiden(adatas, resolution=0.5)
#####

def one_col_lgd(umap):
  legend = umap.legend(bbox_to_anchor=[1.00, 0.5],
  loc='center left', ncol=1, prop={'size': 6})
  legend.get_frame().set_linewidth(0.0)
  for handle in legend.legendHandles:
    handle.set_sizes([25.0])
  return legend

#donor_umap = sc.pl.umap(adatas, color=['sample_name'],show=False, palette=sns.color_palette("husl", 24), legend_fontsize=6, frameon=True, title='Donor')
#lgd = one_col_lgd(donor_umap)
#fig = donor_umap.get_figure()
#fig.set_size_inches(5, 5)
#fig.savefig(str(sc.settings.figdir) + '/umap_lgd_sample', dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')

# by cluster
#leiden_umap = sc.pl.umap(adatas, color=['leiden'], show=False, palette=sns.color_palette("husl", 24),legend_fontsize=6, frameon=True, title='Leiden')
#lgd = one_col_lgd(leiden_umap)
#fig = leiden_umap.get_figure()
#fig.set_size_inches(5, 5)
#fig.savefig(str(sc.settings.figdir) + '/umap_lgd_leiden', dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')

#harmony
sce.pp.harmony_integrate(adatas, 'Animal_ID') 
adatas.obsm['X_pca'] = adatas.obsm['X_pca_harmony']
sc.pp.neighbors(adatas, n_neighbors=20, n_pcs=50)  # final n_neighbors ==>  
sc.tl.umap(adatas)
sc.tl.leiden(adatas, resolution=0.3)              # final resolution   ==>  

# by sample
donor_umap = sc.pl.umap(adatas, color=['Animal_ID'],show=False, palette=sns.color_palette("husl", 24),legend_fontsize=6, frameon=True, title='Donor')
lgd = one_col_lgd(donor_umap)
fig = donor_umap.get_figure()
fig.set_size_inches(5, 5)
fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_sample2',dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')

# by cluster
leiden_umap = sc.pl.umap(adatas, color=['leiden'],show=False, palette=sns.color_palette("husl", 24),legend_fontsize=6, frameon=True, title='Leiden')
lgd = one_col_lgd(leiden_umap)
fig = leiden_umap.get_figure()
fig.set_size_inches(5, 5)
fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_leiden2', dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')

adatas.write_h5ad("colon_harmony_allcells.h5ad")





#ucd run
ucd.tl.base(adatas)
#ucd.pl.base_clustermap(adatas, groupby = 'leiden', n_top_celltypes=75)
ucd.pl.base_clustermap(adatas, groupby = 'leiden', category = 'raw', n_top_celltypes = 75)
plt.show()

#####
adata_colon_microphage = adatas[adatas.obs.leiden.isin(("1","14","11")), :]
adata_colon_microphage.write_h5ad("colon_macrophages.h5ad")
adata_colon_neuron = adatas[adatas.obs.leiden.isin(("13","17")), :]
adata_colon_neuron.write_h5ad("colon_neuron.h5ad")

adatas = adata_colon_microphage

#harmony
sce.pp.harmony_integrate(adatas, 'Animal_ID')
adatas.obsm['X_pca'] = adatas.obsm['X_pca_harmony']
sc.pp.neighbors(adatas, n_neighbors=20, n_pcs=50)  # final n_neighbors ==>
sc.tl.umap(adatas)
sc.tl.leiden(adatas, resolution=0.3)              # final resolution   ==>

# by sample
donor_umap = sc.pl.umap(adatas, color=['Animal_ID'],show=False, palette=sns.color_palette("husl", 24),legend_fontsize=6, frameon=True, title='Donor')
lgd = one_col_lgd(donor_umap)
fig = donor_umap.get_figure()
fig.set_size_inches(5, 5)
fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_sample2',dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')

# by cluster
leiden_umap = sc.pl.umap(adatas, color=['leiden'],show=False, palette=sns.color_palette("husl", 24),legend_fontsize=6, frameon=True, title='Leiden')
lgd = one_col_lgd(leiden_umap)
fig = leiden_umap.get_figure()
fig.set_size_inches(5, 5)
fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_leiden2', dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')

# marker genes for each cluster
sc.tl.rank_genes_groups(adatas, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adatas, n_genes=25, sharey=False)
markers = adatas.uns['rank_genes_groups']
groups = markers['names'].dtype.names
df = pd.DataFrame( { group + "_" +key[:15]: markers[key][group] for group in groups for key in ['names','scores','pvals','pvals_adj','logfoldchanges']})
df.head()
df.to_csv('marker_gene_pval_colon.csv')


#ucd run
ucd.tl.base(adatas)
#ucd.pl.base_clustermap(adatas, groupby = 'leiden', n_top_celltypes=75)
ucd.pl.base_clustermap(adatas, groupby = 'leiden', category = 'raw', n_top_celltypes = 75)
plt.show()


# worked

celltypes = ucd.utils.assign_top_celltypes(adatas, category = "raw", groupby = "leiden", inplace = False)
#ucd.tl.explain(adatas, celltypes = celltypes, groupby = "leiden", group_n = 64)


# add cell type from unicell
label = "celltype_ucdbase_propagated"
adatas.obs[label] = 'unknown'
adatas.obs.loc[adatas.obs.leiden.isin(("0",)), label] = "navie b cell"
adatas.obs.loc[adatas.obs.leiden.isin(("1",)), label] = "t cell"
adatas.obs.loc[adatas.obs.leiden.isin(("2",)), label] = "t cell"
adatas.obs.loc[adatas.obs.leiden.isin(("3",)), label] = "transit amplifying cell of small intestine"
adatas.obs.loc[adatas.obs.leiden.isin(("4",)), label] = "enterocyte of epithelium of large intestine"
adatas.obs.loc[adatas.obs.leiden.isin(("5",)), label] = "t cell"
adatas.obs.loc[adatas.obs.leiden.isin(("6",)), label] = "transit amplifying cell of small intestine"
adatas.obs.loc[adatas.obs.leiden.isin(("7",)), label] = "b cell"
adatas.obs.loc[adatas.obs.leiden.isin(("8",)), label] = "macrophage"
adatas.obs.loc[adatas.obs.leiden.isin(("9",)), label] = "intestinal tuft cell"
adatas.obs.loc[adatas.obs.leiden.isin(("10",)), label] = "plasma cell"
adatas.obs.loc[adatas.obs.leiden.isin(("11",)), label] = "b cell"
adatas.obs.loc[adatas.obs.leiden.isin(("12",)), label] = "lymphoid cell"
adatas.obs.loc[adatas.obs.leiden.isin(("13",)), label] = "enterocyte of epithelium of large intestine"
adatas.obs.loc[adatas.obs.leiden.isin(("14",)), label] = "mast cell"
adatas.obs.loc[adatas.obs.leiden.isin(("15",)), label] = "unknown"

sc.pl.umap(adatas, color = 'celltype_ucdbase_propagated', legend_loc = 'on data', legend_fontsize = "small", frameon = False)

adatas.write_h5ad("IEC_harmony_unicell.h5ad") 



