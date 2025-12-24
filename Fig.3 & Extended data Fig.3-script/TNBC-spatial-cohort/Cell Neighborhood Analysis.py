import scanpy as sc
import pandas as pd
import squidpy as sq
import matplotlib.pyplot as plt
adata = sc.read_visium("/sEA6/")

adata_obs_add = pd.read_csv(r"sEA6_obs.csv", header=None)
adata_obs_add.columns = adata_obs_add.iloc[0]
adata_obs_add = adata_obs_add.set_index(adata_obs_add.columns[0])
adata_obs_add = adata_obs_add.iloc[1:]

common_rows = adata.obs.index.intersection(adata_obs_add.index)

adata = adata[common_rows, :].copy()

for col in adata_obs_add.columns:
    adata.obs[col] = adata_obs_add.loc[adata.obs.index, col].values

cell_subtype_colors = {
    'Fibroblast': '#FF7F0EFF', 'IM.MES': '#2CA02CFF', 'OXPHOS': '#D62728FF', 'B_cell': '#9467BDFF',
    'C5_CCL18': '#8C564BFF', 'C6_MMP9': '#E377C2FF', 'T_cell': '#7F7F7FFF', 'Endothelial_cell': '#BCBD22FF',
    'LAR': '#17BECFFF', 'ASDC': '#5050FFFF', 'C1_CD14': '#CE3D32FF', 'C10_LYVE1': '#749B58FF',
    'C3_CCR2': '#F0E685FF', 'C4_CCL2': '#466983FF', 'C7_CX3CR1': '#BA6338FF', 'C8_MT1G': '#5DB1DDFF',
    'C9_SLC2A1': '#802268FF', 'CD4_EM': '#6BD76BFF', 'CD4_N': '#D595A7FF', 'CD4_REG': "#F5E8E2FF",
    'CD4_REG_Proliferating': '#837B8DFF', 'CD8_EM': '#C75127FF', 'CD8_EMRA': '#D58F5CFF', 'LAR_OXPHOS': '#7A65A5FF',
    'CD8_N': '#E4AF69FF', 'CD8_RM': '#3B1B53FF', 'cDC1': '#CDDEB7FF', 'cDC2': '#612A79FF',
    'IM.regulatory': '#AE1F63FF', 'LanghDC': '#E7C76FFF', 'mregDC_migDC': '#5A655EFF', 'Myeloid_cell': '#CC9900FF',
    'NK_CYTO': '#99CC00FF', 'pDC': '#A9A9A9FF', 'IM.MES_LAR': '#FFB347FF', 'IM.MES_OXPHOS': '#FFFF00',
    'IM.MES_IM.regulatory': '#FF6B6BFF', 'CD8_EX': '#C4A0D3FF', 'IM.regulatory_LAR': '#FFE66DFF',
    'IM.regulatory_OXPHOS': '#7DDDFAFF'
}

adata_filtered.obs["fist_cell"] = adata_filtered.obs["fist_cell"].astype("category")
categories = sorted(adata_filtered.obs["fist_cell"].cat.categories.tolist())

colors_list = []
missing_cells = []

for cell_type in categories:
    if cell_type in cell_subtype_colors:
        colors_list.append(cell_subtype_colors[cell_type])
    else:
        # 如果不在字典中，使用灰色
        colors_list.append('#CCCCCC')  # 灰色
        missing_cells.append(cell_type)

adata_filtered.uns['fist_cell_colors'] = colors_list


sq.gr.spatial_neighbors(adata_filtered)
adata_filtered.obs['fist_cell'] = adata_filtered.obs['fist_cell'].astype('category')
sq.gr.nhood_enrichment(
    adata_filtered,
    cluster_key="fist_cell"
)
result = sq.gr.nhood_enrichment(
    adata_filtered,
    cluster_key="fist_cell",
    copy=True
)

zscore_matrix = result.zscore
counts_matrix = result.counts

sq.pl.nhood_enrichment(adata_filtered, cluster_key="fist_cell")


fig, ax = plt.subplots(figsize=(5, 5))
sq.pl.nhood_enrichment(adata_filtered,
                       cluster_key="fist_cell",
                       method="average",
                       ax=ax)
plt.savefig(r'nhood_enrichment.pdf',
            dpi=300,
            bbox_inches='tight',
            format='pdf')
plt.close()

sq.gr.interaction_matrix(adata_filtered, cluster_key="fist_cell")

adata_filtered.uns["fist_cell_interactions"]
sq.pl.interaction_matrix(adata_filtered, cluster_key="fist_cell", method="average", figsize=(5, 5))

sq.gr.co_occurrence(adata_filtered, cluster_key="fist_cell")
sq.pl.co_occurrence(adata_filtered, cluster_key="fist_cell", clusters="OXPHOS", figsize=(8, 5))