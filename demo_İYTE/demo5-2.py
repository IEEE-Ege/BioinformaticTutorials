# -*- coding: utf-8 -*-
"""Scanpy Analysis App (Purple and Violet Themed)"""

from shiny import App, ui, render, reactive
from shiny.types import FileInfo
import matplotlib.pyplot as plt
import scanpy as sc

# UI
app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.panel_title("SCA-Web"),
            ui.input_file("file", "Upload .h5ad file", accept=[".h5ad"]),
            ui.input_slider("min_genes", "Minimum genes per cell", min=0, max=500, value=200),
            ui.input_slider("max_genes", "Maximum genes per cell", min=500, max=10000, value=5000),
            ui.input_select("hvg_method", "Highly Variable Gene Selection Method",
                            choices=["seurat_v3", "seurat", "cell_ranger"],
                            selected="seurat_v3"),
            ui.input_slider("n_top_genes", "Top HVGs to keep", min=100, max=5000, value=2000),
            ui.input_select("reduction_method", "Dimensionality Reduction Method", choices=["pca"], selected="pca"),
            ui.input_slider("perplexity", "t-SNE Perplexity", min=5, max=50, value=30),
            ui.input_action_button("run", "Run Analysis", class_="btn-violet"),  # Violet button class
            ui.panel_title("Clustering"),
            ui.input_select("clustering_method", "Clustering Algorithm",
                            choices=["leiden", "louvain"], selected="leiden"),
            ui.input_slider("resolution", "Resolution", min=0.1, max=2.0, value=0.8),
            ui.input_action_button("find_clusters", "Find Clusters", class_="btn-violet")  # Violet button class
        ),
        [
            ui.card(
                ui.card_header("Dataset Summary"),
                ui.output_text_verbatim("summary"),
                id="summary_card"
            ),
            ui.card(
                ui.card_header("Highly Variable Genes"),
                ui.output_plot("hvg_plot")
            ),
            ui.card(
                ui.card_header("UMAP Plot"),
                ui.output_plot("umap_plot")
            ),
            ui.card(
                ui.card_header("t-SNE Plot"),
                ui.output_plot("tsne_plot")
            ),
            ui.card(
                ui.card_header("Clustered UMAP Plot"),
                ui.output_plot("clustered_umap_plot")
            ),
            ui.card(
                ui.card_header("Clustered t-SNE Plot"),
                ui.output_plot("clustered_tsne_plot")
            ),
            ui.card(
                ui.card_header("Marker Gene MatrixPlot"),
                ui.output_plot("marker_matrix_plot")
            )
        ]
    ),
    ui.TagList(
        ui.tags.style("""
            /* General Card Styling */
            .card {
                background-color: white;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                border-radius: 8px;
                margin-bottom: 20px;
                padding: 20px;
            }

            /* Header Styling */
            .card-header {
                font-size: 1.2rem;
                font-weight: 600;
                background-color: #4B0082; /* Dark Purple */
                color: white;
                padding: 10px;
                border-top-left-radius: 8px;
                border-top-right-radius: 8px;
            }

            /* Sidebar styling */
            .shiny-input-slider {
                margin-bottom: 15px;
            }
            .shiny-input-select, .shiny-input-file {
                margin-bottom: 15px;
                background-color: #f5f5f5 !important; /* Light Grey with !important to override defaults */
                border-radius: 5px;
            }
            .shiny-output-container {
                margin: 0 auto;
            }
            .shiny-image-output {
                max-width: 700px;
                max-height: 600px;
            }

            /* Button Styling */
            .btn-violet {
                background-color: #6a0dad; /* Deep Violet */
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
            }

            .btn-violet:hover {
                background-color: #9b30ff; /* Lighter Violet for hover */
            }

            /* Text and Summary Styling */
            #summary {
                font-size: 1rem;
                color: #660000; /* Dark Red */
                font-family: 'Arial', sans-serif;
                margin-top: 10px;
            }
            
            /* Image and Plot Styling */
            .shiny-output-container {
                display: flex;
                justify-content: center;
                align-items: center;
            }
        """)
    )
)

# Server logic
def server(input, output, session):

    # NEW Marker Gene Dictionary
    marker_genes_dict = {
        'B-cell': ['CD79A', 'MS4A1'],
        'T-cell': ['CD3D'],
        'T-cell CD8+': ['CD8A', 'CD8B'],
        'NK': ['GNLY', 'NKG7'],
        'Myeloid': ['CST3', 'LYZ'],
        'Monocytes': ['FCGR3A'],
        'Dendritic': ['FCER1A']
    }

    @reactive.calc
    def process_file():
        file = input.file()
        if not file:
            return None

        path = file[0]["datapath"]
        adata = sc.read_h5ad(path)

        # Compute n_genes before filtering
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

        # Filter cells
        sc.pp.filter_cells(adata, min_genes=input.min_genes())
        adata = adata[adata.obs.n_genes < input.max_genes(), :]

        # Normalize and log transform
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.var_names = adata.var_names.astype(str)

        # Set raw data before HVG filtering
        adata.raw = adata

        # HVG selection
        sc.pp.highly_variable_genes(
            adata,
            flavor=input.hvg_method(),
            n_top_genes=input.n_top_genes()
        )
        adata = adata[:, adata.var.highly_variable]

        # Scaling and PCA
        sc.pp.scale(adata, max_value=10)
        sc.pp.pca(adata)

        # Neighbors and embeddings
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.tsne(adata, perplexity=input.perplexity())

        return adata

    @reactive.calc
    def find_clusters():
        adata = process_file()
        if adata is None:
            return None

        method = input.clustering_method()
        resolution = input.resolution()

        if method == "louvain":
            sc.tl.louvain(adata, resolution=resolution, random_state=42)
        elif method == "leiden":
            sc.tl.leiden(adata, resolution=resolution, random_state=42)

        return adata

    @output
    @render.text
    def summary():
        adata = process_file()
        if adata is None:
            return "Upload a file and click Run."
        return f"{adata.shape[0]} cells, {adata.shape[1]} genes after filtering."

    @output
    @render.plot
    def hvg_plot():
        adata = process_file()
        if adata is None:
            return
        plt.figure()
        sc.pl.highly_variable_genes(adata, show=False)
        return plt.gcf()

    @output
    @render.plot
    def umap_plot():
        adata = process_file()
        if adata is None:
            return
        plt.figure()
        if "cell_type" in adata.obs.columns:
            sc.pl.umap(adata, color="cell_type", show=False)
        else:
            sc.pl.umap(adata, color="n_genes", show=False)
        return plt.gcf()

    @output
    @render.plot
    def tsne_plot():
        adata = process_file()
        if adata is None:
            return
        plt.figure()
        if "cell_type" in adata.obs.columns:
            sc.pl.tsne(adata, color="cell_type", show=False)
        else:
            sc.pl.tsne(adata, color="n_genes", show=False)
        return plt.gcf()

    @output
    @render.plot
    def clustered_umap_plot():
        adata = find_clusters()
        if adata is None:
            return
        plt.figure()
        clustering_key = "leiden" if "leiden" in adata.obs.columns else "louvain"
        sc.pl.umap(adata, color=clustering_key, show=False)
        return plt.gcf()

    @output
    @render.plot
    def clustered_tsne_plot():
        adata = find_clusters()
        if adata is None:
            return
        plt.figure()
        clustering_key = "leiden" if "leiden" in adata.obs.columns else "louvain"
        sc.pl.tsne(adata, color=clustering_key, show=False)
        return plt.gcf()

    @output
    @render.plot
    def marker_matrix_plot():
        adata = find_clusters()
        if adata is None:
            plt.figure()
            plt.text(0.5, 0.5, 'No data loaded.', ha='center', va='center')
            plt.axis('off')
            return plt.gcf()

        if adata.raw is None:
            plt.figure()
            plt.text(0.5, 0.5, 'No raw data available.', ha='center', va='center')
            plt.axis('off')
            return plt.gcf()

        plt.figure()
        clustering_key = "leiden" if "leiden" in adata.obs.columns else "louvain"

        # Ensure all genes are strings
        fixed_marker_genes_dict = {
            group: [str(gene) for gene in genes]
            for group, genes in marker_genes_dict.items()
        }

        available_genes = set(adata.var_names)

        filtered_marker_genes_dict = {
            cluster: [gene for gene in genes if gene in available_genes]
            for cluster, genes in fixed_marker_genes_dict.items()
        }

        sc.pl.matrixplot(
            adata,
            filtered_marker_genes_dict,
            groupby=clustering_key,
            use_raw=False,
            show=False
        )
        return plt.gcf()

app = App(app_ui, server)
