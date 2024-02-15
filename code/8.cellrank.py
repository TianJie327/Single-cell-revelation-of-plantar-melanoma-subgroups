# import sys
# 
# # if "google.colab" in sys.modules:
# #     !pip install -q git+https://github.com/theislab/cellrank@dev
    
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

adata = sc.read_h5ad("/gdata01/user/tianjie/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/7.cellrank/MMscvelo.h5ad")
adata.obs[celltype] = adata.obs[celltype].astype('str')
sc.pl.umap(adata,color=celltype, legend_loc='on data',legend_fontsize=8, save=celltype+'.pdf')

adata.obs['celltype'] = adata.obs[celltype].copy()
adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'],ordered=True)
adata.obs['celltype']
celltype='celltype'
sc.pl.umap(adata, color=['celltype'],legend_loc='on data',legend_fontsize=8, save=celltype+'.pdf')



adata.X = adata.raw.to_adata().X.copy()
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata, enforce=True)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata, n_jobs=20)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(
    adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, save=celltype+'embedding_stream_1.pdf'
)
scv.pl.velocity_embedding_stream(adata, basis="umap", color=celltype, save=celltype+'embedding_stream_2.pdf'
)
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color=celltype,legend_fontsize=12, title="", smooth= 0.5, min_mass=0.5, save=celltype+'embedding_stream_3.pdf'
)

cr.tl.terminal_states(adata, cluster_key=celltype, weight_connectivities=0.2)
cr.pl.terminal_states(adata, save=celltype+'terminal_states1.pdf')
cr.pl.terminal_states(adata, same_plot=False, save=celltype+'terminal_states2.pdf')
cr.pl.terminal_states(adata, discrete=True, save=celltype+'terminal_states3.pdf')
cr.tl.initial_states(adata, cluster_key=celltype)
cr.pl.initial_states(adata, discrete=True, save=celltype+'initial_states.pdf')
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False, save=celltype+'lineages.pdf')
cr.pl.lineages(adata, same_plot=True, save=celltype+'absorption.pdf')
scv.tl.recover_latent_time(
    adata, root_key="initial_states_probs", end_key="terminal_states_probs"
)
scv.tl.paga(
    adata,
    groups=celltype,
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime",
)
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key=celltype,
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA",
    save=celltype+'paga_pie.pdf'
)
