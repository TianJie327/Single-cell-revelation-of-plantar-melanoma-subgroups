
#https://pyscenic.readthedocs.io/en/latest/tutorial.html
#分析环境 source activate pyscenic_py3.1

import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 20

RESOURCES_FOLDERNAME = "/gdata01/AnnotationPackage/pyscenic/resources/"
AUXILLIARIES_FOLDERNAME = "/gdata01/AnnotationPackage/pyscenic"
RESULTS_FOLDERNAME = "results/"
FIGURES_FOLDERNAME = "figures/"

sc.settings.figdir = 'results/'

BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"


def savesvg(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:
    """
    Save figure as vector-based SVG image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format='svg')
    
    
def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
    
    
    
################################### 三、二进制热图 ###################################
# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scanpy.plotting._tools.scatterplots import plot_scatter
import seaborn as sns


def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f
    
import matplotlib.pyplot as plt


ad=sc.read_h5ad("MM.h5ad")
sc.pl.umap(ad,color=celltype, legend_loc='on data',legend_fontsize=18, save=celltype+'.pdf')

ad.obs['celltype'] = ad.obs[celltype].copy()
ad.obs['celltype'] = pd.Categorical(ad.obs['celltype'],categories=["MM0", "MM1", "MM2", "MM3", "MM4"],ordered=True)
ad.obs['celltype']
celltype='celltype'
sc.pl.umap(ad, color=['celltype'],legend_loc='on data',legend_fontsize=18, save=celltype+'.pdf')


ad.obs.celltype.value_counts()
cells = ad.obs.celltype.sort_values()
cells.index

# scenic output
lf = lp.connect( "MM_SCENIC.loom", mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

lf.close()

thresholds = pd.read_csv('/gpfs/share/home/2011110275/cellranger.7.0.1/pyscenic/500_5000_10/2.MMsub_copykat/thresholds.csv', index_col=0).threshold
thresholds.adjust = thresholds.copy()
thresholds.adjust = thresholds.adjust.astype('float64')
thresholds.adjust = thresholds.adjust.to_dict()
thresholds.adjust

bin_mtx2, thresholds2 = binarize(auc_mtx, threshold_overides = thresholds.adjust, num_workers=30) 
ad.obs.celltype.value_counts()
cells = ad.obs.celltype.sort_values()
cells.index

#COLORS = sns.color_palette("hls", 5) 
COLORS = ["#92140C",
  "#000000",
  "#EE964B",
  "#CFF27E",
  "#A09BE7"]

COLORS
cell_type_color_lut = dict(zip(sorted(list(set(ad.obs.celltype))), COLORS)) 
cell_id2cell_type_lut = ad.obs.celltype.to_dict() 

#黑白色盘
bw_palette = sns.xkcd_palette(["white", "black"]) 
sns.set()
sns.set_style("whitegrid")
fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])
savesvg('legend - GSE115978 - on_off.svg', fig)

#热图表头
sns.set()
sns.set(font_scale=0.8)
fig = palplot(sns.color_palette(COLORS), sorted(list(set(ad.obs.celltype))), size=1.0)
savesvg('legend - cell_type_colors.svg', fig)

#二进制热图
sns.set()
sns.set(font_scale=1.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})

g = sns.clustermap(bin_mtx2.T[cells.index], col_cluster=False,
               col_colors=cells.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
               cmap=bw_palette, figsize=(8,8))
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cell Type'])
g.cax.set_visible(False)
g.fig.savefig('fig2h.pdf', format='pdf')

TF = [
 'NFIB(+)',
 'NFKB1(+)', 
 'MYC(+)',
 'EGR3(+)',
 'TWIST1(+)',
 'TCF4(+)',
 'HMGA2(+)'
        ]
sc.pl.umap(ad_v1,color=TF, legend_loc='on data',legend_fontsize=16, save='TFs.pdf', ncols= 1)

fig = plt.figure(figsize=(6,30))
ax1 = plt.subplot(711) 
plot_binarization(auc_mtx, 'NFIB(+)', thresholds.adjust['NFIB(+)'])
ax1.set_title('NFIB(+)')

ax2 = plt.subplot(712) 
plot_binarization(auc_mtx, 'NFKB1(+)', thresholds.adjust['NFKB1(+)'])
ax2.set_title('NFKB1(+)')

ax3 = plt.subplot(713) 
plot_binarization(auc_mtx, 'MYC(+)', thresholds.adjust['MYC(+)'])
ax3.set_title('MYC(+)')

ax4 = plt.subplot(714) 
plot_binarization(auc_mtx, 'EGR3(+)', thresholds.adjust['EGR3(+)'])
ax4.set_title('EGR3(+)')

ax5 = plt.subplot(715) 
plot_binarization(auc_mtx, 'TWIST1(+)', thresholds.adjust['TWIST1(+)'])
ax5.set_title('TWIST1(+)')

ax6 = plt.subplot(716) 
plot_binarization(auc_mtx, 'TCF4(+)', thresholds.adjust['TCF4(+)'])
ax6.set_title('TCF4(+)')

ax7 = plt.subplot(717) 
plot_binarization(auc_mtx, 'HMGA2(+)', thresholds.adjust['HMGA2(+)'])
ax7.set_title('HMGA2(+)')

plt.tight_layout()
plt.savefig("my_figure_tfs.pdf")
plt.show()