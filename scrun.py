import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# scanpy settings
sc.settings.verbosity = 3
sc.set_figure_params(dpi_save=300, color_map='plasma', vector_friendly=True, fontsize=14)
sc.logging.print_versions()

# set rcParams for matplotlib
plt.rcParams['figure.figsize'] = (8,8)
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['savefig.dpi'] = 300

# setting folder path to save analysis results
results = './results/'

# reading the count matrix
nestorowa_data = sc.read('./data/nestorowa_corrected_log2_transformed_counts.txt', cache=True)
nestorowa_data.raw = nestorowa_data
print(nestorowa_data)

# apply the QC recipe
sc.pp.recipe_weinreb17(nestorowa_data, log=False)

# examine the QC annotations.
print(nestorowa_data)

# Calculate the principal components
sc.tl.pca(nestorowa_data, svd_solver='arpack')

# Visualize the data using a PCA plot
sc.pl.pca(nestorowa_data, na_color='#ffcc00',add_outline=True, outline_color=('black','#b85323'), outline_width=(0.08,0.05), title='Nestorowa PCA Plot', save='Figure1.png', show=False)

# Elbow plot to visualize the variance ratio for each principal component
sc.pl.pca_variance_ratio(nestorowa_data, n_pcs=18, save="Figure2.png", show=False)

# Compute a neighborhood graph of observations using pcs
sc.pp.neighbors(nestorowa_data, n_neighbors=4, n_pcs=20)

# Apply diffusion map algorithm to denoise data
sc.tl.diffmap(nestorowa_data)

# Compute a neighborhood graph of observations using diffmap representation
sc.pp.neighbors(nestorowa_data, n_neighbors=10, use_rep='X_diffmap')

# Compare UMAP, TSNE, Force-directed graph data visualizations.
sc.tl.umap(nestorowa_data)
sc.tl.tsne(nestorowa_data, n_pcs=20)
sc.tl.draw_graph(nestorowa_data, layout='fa', random_state=1)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,12))
fig.tight_layout()
sc.pl.pca(nestorowa_data, na_color='#ffcc00',add_outline=True, outline_color=('black','#b85323'), outline_width=(0.08,0.05), ax=fig.add_subplot(axes[0, 0]), show=False)
sc.pl.tsne(nestorowa_data, na_color='#ffcc00',add_outline=True, outline_color=('black','#b85323'), outline_width=(0.08,0.05), ax=fig.add_subplot(axes[1, 0]), show=False)
sc.pl.umap(nestorowa_data, na_color='#ffcc00',add_outline=True, outline_color=('black','#b85323'), outline_width=(0.08,0.05), ax=fig.add_subplot(axes[0, 1]), show=False)
sc.pl.draw_graph(nestorowa_data, na_color='#ffcc00',add_outline=True, outline_color=('black','#b85323'), outline_width=(0.08,0.05), ax=fig.add_subplot(axes[1, 1]), show=False)
# plt.show()
fig.savefig('./figures/Figure3.png', dpi=300)

# Save the processed data
# nestorowa_data.write(results + "nestorowa_data_processed.h5ad")

# Apply leiden algorithm for community detection
sc.tl.leiden(nestorowa_data, resolution=1, key_added='Leiden Clusters (r=1)')
sc.tl.leiden(nestorowa_data, resolution=0.5, key_added='Leiden Clusters (r=0.5)')

# Investigate clusters with known marker gene expression in at two leiden resolutions
sc.tl.umap(nestorowa_data)
sc.pl.umap(nestorowa_data, color=['Leiden Clusters (r=1)', 'Leiden Clusters (r=0.5)'], legend_loc='on data', palette='tab20', save='Figure4-a.png', add_outline=True, outline_color=('black','#b85323'), outline_width=(0.08,0.05), show=False)
sc.pl.umap(nestorowa_data, color=['Leiden Clusters (r=0.5)', 'Procr', 'Hba-a2', 'Elane'], legend_loc='on data', palette='tab20', save='Figure4-b.png', show=False)
sc.pl.umap(nestorowa_data, color=['Leiden Clusters (r=0.5)', 'Irf8', 'Itga2b', 'Prss34'], legend_loc='on data', palette='tab20', save='Figure4-c.png', show=False)
sc.pl.umap(nestorowa_data, color=['Leiden Clusters (r=0.5)', 'Cd19', 'Ms4a2', 'Dntt'], legend_loc='on data', palette='tab20', save='Figure4-d.png', show=False)


# Lookup marker gene expression in each cluster based on table1
fig, axes = plt.subplots(ncols=3,figsize=(15,5))
fig.tight_layout()
sc.pl.dotplot(nestorowa_data, ['Procr', 'Hba-a2', 'Elane', 'Irf8', 'Itga2b', 'Prss34', 'Cd19', 'Cma1', 'Ms4a2', 'Dntt'], groupby='Leiden Clusters (r=0.5)', cmap='plasma', show=False, ax=fig.add_subplot(axes[0]))
sc.pl.matrixplot(nestorowa_data, ['Procr', 'Hba-a2', 'Elane', 'Irf8', 'Itga2b', 'Prss34', 'Cd19', 'Cma1', 'Ms4a2','Dntt'], groupby='Leiden Clusters (r=0.5)', cmap='plasma', standard_scale='var', show=False, ax=fig.add_subplot(axes[1]))
sc.pl.stacked_violin(nestorowa_data, ['Procr', 'Hba-a2', 'Elane', 'Irf8', 'Itga2b', 'Prss34', 'Cd19','Cma1', 'Ms4a2', 'Dntt'], groupby='Leiden Clusters (r=0.5)', cmap='plasma', show=False, ax=fig.add_subplot(axes[2]))
# plt.show()
fig.savefig('./figures/figure5.png', dpi=300)

# Heatmap to visualize global marker gene expression
sc.pl.heatmap(nestorowa_data, ['Procr', 'Hba-a2', 'Elane', 'Irf8', 'Itga2b', 'Prss34', 'Cd19', 'Cma1', 'Ms4a2', 'Dntt'], groupby='Leiden Clusters (r=0.5)', cmap='plasma', swap_axes=True, dendrogram=True, save="Figure6.png", show=False)

# Force-directed graph visualization of marker genes at Leiden r=0.5 resolution
sc.tl.draw_graph(nestorowa_data, layout='fa', random_state=1)

sc.pl.draw_graph(nestorowa_data, color=['Leiden Clusters (r=0.5)', 'Procr', 'Hba-a2', 'Elane'], legend_loc='on data', save="Figure7-a.png", show=False)
sc.pl.draw_graph(nestorowa_data, color=['Leiden Clusters (r=0.5)', 'Irf8', 'Itga2b','Prss34'], legend_loc='on data', save="Figure7-b.png", show=False)
sc.pl.draw_graph(nestorowa_data, color=['Leiden Clusters (r=0.5)', 'Cd19', 'Ms4a2', 'Dntt'], legend_loc='on data', save="Figure7-c.png", show=False)

# Annotate clusters based on Leiden r=0.5 and marker gene expression
nestorowa_data.obs['Clusters'] = nestorowa_data.obs['Leiden Clusters (r=0.5)']
nestorowa_data.obs['Clusters'].cat.categories = ['0', '1', '2', '3/Erythroid', '4/HSCs', '5', '6', '7', '8', '9', '10', '11', '12/Neutrophils','13/Monocytes', '14/Megakaryocytes', '15/Basophils', '16/Bcell', '17', '18']

# Read in the published cell types to compare
reference_cell_types = pd.read_csv('./data/nestorowa_corrected_population_annotation.txt', delimiter=' ')

# Extract published clusters annotations for each cell
reference_cell_types = [reference_cell_types.loc[cell, 'celltype']
                        if cell in reference_cell_types.index else 'other' for cell in nestorowa_data.obs_names]

# Apply annotations to cells in a new column
nestorowa_data.obs['Reference Cell Types'] = reference_cell_types

# compare our annotations to the published using force-directed graph
sc.pl.draw_graph(nestorowa_data, color=['Clusters', 'Reference Cell Types'],na_color='#ffcc00', add_outline=True, outline_color=('black', '#6c7373'),outline_width=(0.08,0.05),
                legend_loc="on data", legend_fontsize='small', palette='Set1', save="Figure8.png", show=False)

# compare our annotations to the published using umap
sc.pl.umap(nestorowa_data, color=['Clusters', 'Reference Cell Types'],na_color='#ffcc00', add_outline=True, outline_color=('black', '#6c7373'),outline_width=(0.08,0.05), legend_loc='on data', legend_fontsize='x-small', palette='tab20', save="Figure9.png", show=False)

# Choose a starting cell for the trajectory analysis
root_cell = np.where(nestorowa_data.obs['Reference Cell Types'] == 'ESLAM')[0][0]
nestorowa_data.uns['iroot'] = root_cell

# Compute diffusion pseudotime
sc.tl.dpt(nestorowa_data)

# using a different name for the dpt_pseudotime key so it shows up clearly on the plot.
nestorowa_data.obs['Diffusion Pseudotime'] = nestorowa_data.obs['dpt_pseudotime']

# visualize the diffusion pseudotime ordering.
sc.pl.umap(nestorowa_data, color=['Diffusion Pseudotime', 'Reference Cell Types'], palette='Set1',legend_loc='on data', save="Figure10.png",legend_fontsize='small', show=False)
sc.pl.draw_graph(nestorowa_data, color=['Diffusion Pseudotime', 'Reference Cell Types'], palette='Set1', legend_loc='on data',legend_fontsize='x-small', save="Figure11.png", show=False)

# Generate paga representation
sc.tl.paga(nestorowa_data, groups='Clusters')

# Look at marker genes along some trajectories - genes for different stages.
# The first 5 genes mark different stages of the erythroid trajectory - going from early to late. The remaining genes neutrophil (Elane, Cebpe).
trajectory_genes = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
                    'Elane', 'Cebpe', 'Gfi1']                    # neutrophil

# Lets create paths of clusters to visualize the expression of trajectory_genes
paths = [('Erythrocytes', ['4/HSCs', '5', '8', '0', '2', '1', '11', '10', '3/Erythroid']),
         ('Neutrophils', ['4/HSCs', '5', '8','0', '2', '7', '12/Neutrophils'])]

# using a different name for the dpt_pseudotime key so it shows up pretty on the plot.
nestorowa_data.obs['Pseudotime'] = nestorowa_data.obs['dpt_pseudotime']

# Plot the heatmap for gene expression changes over pseudotime
fig, axes = plt.subplots(ncols=2, figsize=(12,6))
plt.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for i, (celltype, pseudopath) in enumerate(paths):
    fig = sc.pl.paga_path(
        nestorowa_data, pseudopath, trajectory_genes,
        show_node_names=False,
        ax=axes[i],
        left_margin=0.15,
        n_avg=50,
        ytick_fontsize=12,
        annotations=['Pseudotime'],
        show_yticks=True if i == 0 else False,
        show_colorbar=False,
        color_map='YlGn',
        color_maps_annotations={'Pseudotime': 'plasma'},
        title='{} path'.format(celltype),
        show=False)
plt.savefig('./figures/Figure12.png', dpi=300)
# plt.show()

# save the annotated/analazyed data.
nestorowa_data.write(results + "nestorowa_final.h5ad")
print(nestorowa_data)
