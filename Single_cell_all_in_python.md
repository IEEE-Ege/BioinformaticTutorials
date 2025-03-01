**Contents**

-   [Contents](https://github.com/seandavi/awesome-single-cell#contents)

-   [Software
    packages](https://github.com/seandavi/awesome-single-cell#software-packages)

    -   [RNA-seq](https://github.com/seandavi/awesome-single-cell#rna-seq)

    -   [Quality
        control](https://github.com/seandavi/awesome-single-cell#quality-control)

    -   [Marker and differential gene expression
        identification](https://github.com/seandavi/awesome-single-cell#marker-and-differential-gene-expression-identification)

    -   [Cell
        clustering](https://github.com/seandavi/awesome-single-cell#cell-clustering)

    -   [Dimension
        reduction](https://github.com/seandavi/awesome-single-cell#dimension-reduction)

    -   [Archetypal
        analysis](https://github.com/seandavi/awesome-single-cell#archetypal-analysis)

    -   [Batch-effect
        removal](https://github.com/seandavi/awesome-single-cell#batch-effect-removal)

    -   [Cell projection and unimodal
        integration](https://github.com/seandavi/awesome-single-cell#cell-projection-and-unimodal-integration)

    -   [Simulation](https://github.com/seandavi/awesome-single-cell#simulation)

    -   [Pseudotime and trajectory
        inference](https://github.com/seandavi/awesome-single-cell#pseudotime-and-trajectory-inference)

    -   [Cell type identification and
        classification](https://github.com/seandavi/awesome-single-cell#cell-type-identification-and-classification)

    -   [Doublet
        Identification](https://github.com/seandavi/awesome-single-cell#doublet-identification)

    -   [Cell
        subsampling](https://github.com/seandavi/awesome-single-cell#cell-subsampling)

    -   [Feature (Gene)
        imputation](https://github.com/seandavi/awesome-single-cell#feature-gene-imputation)

    -   [Variant
        calling](https://github.com/seandavi/awesome-single-cell#variant-calling)

    -   [Epigenomics](https://github.com/seandavi/awesome-single-cell#epigenomics)

    -   [Multi-assay data
        integration](https://github.com/seandavi/awesome-single-cell#multi-assay-data-integration)

    -   [Rare cell
        detection](https://github.com/seandavi/awesome-single-cell#rare-cell-detection)

    -   [Single cell large
        model](https://github.com/seandavi/awesome-single-cell#single-cell-large-model)

    -   [Other
        applications](https://github.com/seandavi/awesome-single-cell#other-applications)

    -   [Spatial
        transcriptomics](https://github.com/seandavi/awesome-single-cell#spatial-transcriptomics)

-   [Tutorials and
    workflows](https://github.com/seandavi/awesome-single-cell#tutorials-and-workflows)

-   [Web portals, apps, and
    databases](https://github.com/seandavi/awesome-single-cell#web-portals-apps-and-databases)

    -   [Web portals and
        databases](https://github.com/seandavi/awesome-single-cell#web-portals-and-databases)

    -   [Interactive visualization and
        analysis](https://github.com/seandavi/awesome-single-cell#interactive-visualization-and-analysis)

    -   [Big data approach
        overview](https://github.com/seandavi/awesome-single-cell#big-data-approach-overview)

    -   [Experimental
        design](https://github.com/seandavi/awesome-single-cell#experimental-design)

-   [Similar lists and
    collections](https://github.com/seandavi/awesome-single-cell#similar-lists-and-collections)

-   [People](https://github.com/seandavi/awesome-single-cell#people)

    -   [Female](https://github.com/seandavi/awesome-single-cell#female)

    -   [Male](https://github.com/seandavi/awesome-single-cell#male)

**Software packages**

**RNA-seq**

-   [anchor](https://github.com/yeolab/anchor) - \[Python\] - ⚓ Find
    bimodal, unimodal, and multimodal features in your data

-   [bonvoyage](https://github.com/yeolab/bonvoyage) - \[Python\] - 📐
    Transform percentage-based units into a 2d space to evaluate changes
    in distribution with both magnitude and direction.

-   [Cell_BLAST](https://github.com/gao-lab/Cell_BLAST) - \[Python\] - A
    BLAST-like toolkit for scRNA-seq data querying and automated
    annotation.

-   [CellCNN](https://github.com/eiriniar/CellCnn) - \[Python\] -
    Representation Learning for detection of phenotype-associated cell
    subsets

-   [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) -
    \[Linux Binary\] - Cell Ranger is a set of analysis pipelines that
    process Chromium single-cell RNA-seq output to align reads, generate
    gene-cell matrices and perform clustering and gene expression
    analysis. *Software requires registration with 10xgenomics.*

-   [Clustergrammer](https://github.com/maayanlab/clustergrammer) -
    \[Python, JavaScript\] - Interative web-based heatmap for
    visualizing and analyzing high dimensional biological data,
    including single-cell RNA-seq. Clustergrammer can be used within a
    Jupyter notebook as an interative widget that can be shared using
    GitHub and NBviewer, see [example
    notebook](http://nbviewer.jupyter.org/github/MaayanLab/CCLE_Clustergrammer/blob/master/notebooks/Clustergrammer_CCLE_Notebook.ipynb).

-   [Clustergrammer2](https://github.com/ismms-himc/clustergrammer2) -
    \[Python, JavaScript\] - Interative WebGL web-based heatmap for
    visualizing and analyzing single-cell high-dimensional and
    location-based biological data. Clustergrammer can be used within a
    Jupyter notebook as an interative widget that can be shared using
    GitHub and NBviewer, see [case
    studies](https://clustergrammer.readthedocs.io/case_studies.html).

-   [cyclum](https://github.com/KChen-lab/cyclum) - \[python\] - Cyclum
    is a novel AutoEncoder approach that characterizes circular
    trajectories in the high-dimensional gene expression space. Applying
    Cyclum to removing cell-cycle effects leads to substantially
    improved delineations of cell subpopulations, which is useful for
    establishing various cell atlases and studying tumor
    heterogeneity. [bioRxiv](https://www.biorxiv.org/content/10.1101/625566v1)

-   [dropkick](https://github.com/KenLauLab/dropkick) - \[Python\] -
    Automated cell filtering for single-cell RNA sequencing data.

-   [dynamo](https://github.com/aristoteleo/dynamo-release) -
    \[Python\] - Inclusive model of expression dynamics with scSLAM-seq
    and multiomics, vector field reconstruction and potential landscape
    mapping.

-   [Falco](https://github.com/VCCRI/Falco/) - \[AWS cloud\] - [Falco: A
    quick and flexible single-cell RNA-seq processing framework on the
    cloud](http://www.biorxiv.org/content/early/2016/07/15/064006.abstract).

-   [FastProject](https://github.com/yoseflab/fastproject) -
    \[Python\] - Signature analysis on low-dimensional projections of
    single-cell expression data.

-   [flotilla](https://github.com/yeolab/flotilla) - \[Python\] -
    Reproducible machine learning analysis of gene expression and
    alternative splicing data

-   [GPfates](https://github.com/Teichlab/GPfates) - \[Python\] - Model
    transcriptional cell fates as mixtures of Gaussian Processes

-   [GSEApy](https://github.com/zqfang/GSEApy) - \[Python\] - GSEApy:
    Gene Set Enrichment Analysis in Python. GSEApy is a Python/Rust
    implementation for GSEA and wrapper for Enrichr. GSEApy can be used
    for RNA-seq, ChIP-seq, Microarray data. It can be used for
    convenient GO enrichment and to produce publication quality figures
    in python.

-   [HTSeq](https://github.com/htseq/htseq) - \[Python\] - A Python
    library to facilitate programmatic analysis of data from
    high-throughput sequencing (HTS) experiments. A popular component
    of HTSeq is htseq-count, a script to quantify gene expression in
    bulk and single-cell RNA-Seq and similar experiments.

-   [ICGS](https://github.com/nsalomonis/altanalyze) - \[Python\] -
    Iterative Clustering and Guide-gene Selection (Olsson et al. Nature
    2016). Identify discrete, transitional and mixed-lineage states from
    diverse single-cell transcriptomics platforms. Integrated FASTQ
    pseudoalignment /quantification (Kallisto), differential expression,
    cell-type prediction and optional cell cycle exclusion analyses.
    Specialized methods for processing BAM and 10X Genomics spares
    matrix files. Associated single-cell splicing PSI methods
    (MultIPath-PSI). Apart of the AltAnalyze toolkit along with
    accompanying visualization methods (e.g., heatmap, t-SNE,
    SashimiPlots, network graphs). Easy-to-use graphical user and
    commandline interfaces.

-   [InMoose](https://github.com/epigenelabs/inmoose) - \[Python\] -
    InMoose is the Integrated Multi Omic Open Source Environment. It is
    a collection of tools for the analysis of omic data. Allows for
    batch effect correction, cohort QC, Differential Expression Analysis
    and Consensus Clustering.

-   [ivis](https://github.com/beringresearch/ivis) - \[Python or R\] -
    Structure-preserving dimensionality reduction in single-cell
    datasets.

-   [kb-python](https://github.com/pachterlab/kb_python) - \[Python\]
    - kb-python is a python package for processing single-cell
    RNA-sequencing. It wraps
    the [kallisto \| bustools](https://www.kallistobus.tools) single-cell
    RNA-seq command line tools in order to unify multiple processing
    workflows.

-   [knn-smoothing](https://github.com/yanailab/knn-smoothing) -
    \[python or R or matlab\] - The algorithm is based on the
    observation that across protocols, the technical noise exhibited by
    UMI-filtered scRNA-Seq data closely follows Poisson statistics.
    Smoothing is performed by first identifying the nearest neighbors of
    each cell in a step-wise fashion, based on variance-stabilized and
    partially smoothed expression profiles, and then aggregating their
    transcript counts.

-   [MIMOSCA](https://github.com/asncd/MIMOSCA) - \[python\] - A
    repository for the design and analysis of pooled single cell RNA-seq
    perturbation experiments (Perturb-seq).

-   [nimfa](https://github.com/ccshao/nimfa) - \[Python\] - Nimfa is a
    Python scripting library which includes a number of published matrix
    factorization algorithms, initialization methods, quality and
    performance measures and facilitates the combination of these to
    produce new strategies. The library represents a unified and
    efficient interface to matrix factorization algorithms and methods.

-   [novoSpaRc](https://github.com/rajewsky-lab/novosparc) -
    \[Python\] - Predict locations of single cells in space by solely
    using single-cell RNA sequencing data. An existing reference
    database of marker genes is not required, but significantly enhances
    performance if
    available. [bioRxiv](https://www.biorxiv.org/content/early/2018/10/30/456350).

-   [outrigger](https://github.com/YeoLab/outrigger) - \[Python\] -
    Outrigger is a program to calculate alternative splicing scores of
    RNA-Seq data based on junction reads and a *de novo*, custom
    annotation created with a graph database, especially made for
    single-cell analyses.

-   [PyGMNormalize](https://github.com/ficusss/PyGMNormalize) -
    \[Python\] - Python implementation
    of [edgeR](https://www.ncbi.nlm.nih.gov/pubmed/19910308) normalization
    method for count matrices.

-   [RAPIDS-singlecell](https://github.com/scverse/rapids_singlecell) -
    \[Python\] - A GPU-accelerated tool leveraging RAPIDS for scRNA
    analysis. Seamless scverse compatibility for efficient single-cell
    data processing and analysis. Replcates features from Scanpy, while
    also incorporating select functionalities from Squidpy and
    Decoupler.

-   [rMATS](http://rnaseq-mats.sourceforge.net/) - \[Python\] - RNA-Seq
    Multavariate Analysis of Transcript Splicing.

-   [Scanpy](https://github.com/theislab/scanpy) - \[Python\] - Scanpy
    provides computationally efficient tools that scale up to very large
    data sets and enables simple integration of advanced machine
    learning algorithms.

-   [scbean](https://github.com/jhu99/scbean) - \[Python\] - Scbean
    integrates a range of models for single-cell data analysis,
    including dimensionality reduction, removing batch effects, and
    transferring well-annotated cell type labels from scRNA-seq to
    scATAC-seq and spatial resolved transcriptomics, and joint-analysis
    of paired multimodal single-cell data.

-   [SCCAF](https://github.com/SCCAF/sccaf) - \[Python\] Single Cell
    Clustering Assessment Framework (SCCAF) is a method for automated
    identification of putative cell types from single cell data by
    iteratively applying clustering and a machine learning
    approach. [Putative cell type discovery from single-cell gene
    expression data](https://www.nature.com/articles/s41592-020-0825-9)

-   [schist](https://github.com/dawe/schist) - \[Python\] - schist is a
    scanpy-compatible python library which implements Nested Stochastic
    Block Models to identify cell groups in single cell experiments.

-   [scVI](https://github.com/YosefLab/scVI) - \[python\] - scVI is a
    ready-to-use scalable framework for the probabilistic representation
    and analysis of gene expression in single cells (batch correction,
    visualization, clustering, and differential expression). [Deep
    generative modeling for single-cell
    transcriptomics](https://www.nature.com/articles/s41592-018-0229-2)

-   [scTDA](https://github.com/RabadanLab/scTDA) - \[Python\] - scTDA is
    an object oriented python library for topological data analysis of
    high-throughput single-cell RNA-seq data. It includes tools for the
    preprocessing, analysis, and exploration of single-cell RNA-seq data
    based on topological representations.

-   [scTCRseq](https://github.com/ElementoLab/scTCRseq) - \[python\] -
    Map T-cell receptor (TCR) repertoires from single cell RNAseq.

-   [singlet](https://github.com/iosonofabio/singlet) - \[Python\] -
    Single cell RNA-Seq analysis with phenotypes.

-   [SPRING](https://github.com/AllonKleinLab/SPRING) - \[matlab,
    javascript, python\] - SPRING is a collection of pre-processing
    scripts and a web browser-based tool for visualizing and interacting
    with high dimensional data. SPRING was developed for single cell
    RNA-Seq data but can be applied more generally.

-   [scTOP](https://github.com/Emergent-Behaviors-in-Biology/scTOP) -
    \[Python\] - Single-cell type order parameters. Physics-inspired
    method of processing single-cell RNA-seq and identifying cell fate,
    motivated by the epigenetic landscape.

**Quality control**

-   [gene_network_evaluation](https://github.com/EngreitzLab/gene_network_evaluation/) -
    \[Python\] - A flexible framework to evaluate the plausibility of
    gene programs inferred from single-cell genomic data. The assessment
    is broken down into themes such as goodness of fit (ability to
    explain the data), co-regulation, mechanistic interactions etc.
    Under each theme, multiple evaluation tasks are conceptualised and
    implemented using appropriate statistical tests.

**Marker and differential gene expression identification**

-   [GPseudoClust](https://github.com/magStra/GPseudoClust) -
    \[Python\] - Software that clusters genes for pseudotemporally
    ordered data and quantifies the uncertainty in cluster allocations
    arising from the uncertainty in the pseudotime ordering.

-   [GiniClust](https://github.com/lanjiangboston/GiniClust) -
    \[Python/R\] - GiniClust is a clustering method implemented in
    Python and R for detecting rare cell-types from large-scale
    single-cell gene expression data. GiniClust can be applied to
    datasets originating from different platforms, such as multiplex
    qPCR data, traditional single-cell RNAseq or newly emerging
    UMI-based single-cell RNAseq, e.g. inDrops and Drop-seq.

-   [Phenotype Cover](https://github.com/euxhenh/phenotype-cover) -
    \[Python\] - Provides two algorithms for marker selection (G-PC,
    CEM-PC) introduced in [Multiset multicover methods for
    discriminative marker
    selection](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(22)00229-6).
    Most marker selection methods focus on differential expression (DE)
    analysis. Although such methods work well for data with a few
    non-overlapping marker sets, they are not appropriate for large
    atlas-size datasets where several cell types and tissues are
    considered. To address this, we define the phenotype cover (PC)
    problem for marker selection and present algorithms that can improve
    the discriminative power of marker sets.

**Cell clustering**

-   [BackSPIN](https://github.com/linnarsson-lab/BackSPIN) -
    \[Python\] - Biclustering algorithm developed taking into account
    intrinsic features of single-cell RNA-seq experiments.

-   [dropClust](https://github.com/debsin/dropClust) - \[R/Python\] -
    Efficient clustering of ultra-large scRNA-seq data.

**Dimension reduction**

-   [torchdr](https://github.com/TorchDR/TorchDR) - \[python\] -
    Dimensionality reduction toolbox using PyTorch, featuring various
    algorithms such as TSNE, UMAP, and more. Supports GPU acceleration
    to maximize computational efficiency.

-   [scvis](https://bitbucket.org/jerry00/scvis-dev) - \[python\]
    - [Interpretable dimensionality reduction of single cell
    transcriptome data with deep generative
    models](https://doi.org/10.1101/178624)

-   [ZIFA](https://github.com/epierson9/ZIFA) - \[Python\] -
    Zero-inflated dimensionality reduction algorithm for single-cell
    data.

-   [scPRINT](https://github.com/cantinilab/scPRINT) - \[python\] -
    scPRINT is pretrained on 50M cells and generates multiple cell
    embeddings from single cell RNAseq profiles. [scPRINT: pre-training
    on 50 million cells allows robust gene network
    predictions](https://www.biorxiv.org/content/10.1101/2024.07.29.605556v1)

**Archetypal analysis**

-   [scAAnet](https://github.com/AprilYuge/scAAnet_latest) -
    \[Python\] - scAAnet performs non-linear archetypal analysis through
    autoencoder networks to identify shared gene expression programs
    (GEPs) among heterogenous cell populations and infer relative
    activity of each GEP across cells.

**Batch-effect removal**

-   [BatchEffectRemoval](https://github.com/ushaham/BatchEffectRemoval) -
    \[Python\] - [Removal of Batch Effects using Distribution-Matching
    Residual Networks](https://doi.org/10.1093/bioinformatics/btx196)

-   [ResPAN](https://github.com/AprilYuge/ResPAN) - \[Python\] - ResPAN
    is a light structured **Res**idual autoencoder and mutual nearest
    neighbor **P**aring guided **A**dversarial **N**etwork for scRNA-seq
    batch correction.

-   [TASC](https://github.com/scrna-seq/TASC) - \[C++, python\] - To
    account for cell-to-cell technical differences, we propose a
    statistical framework, TASC (Toolkit for Analysis of Single Cell
    RNA-seq), an empirical Bayes approach to reliably model the
    cell-specific dropout rates and amplification bias by use of
    external RNA spike-ins. TASC incorporates the technical parameters,
    which reflect cell-to-cell batch effects, into a hierarchical
    mixture model to estimate the biological variance of a gene and
    detect differentially expressed genes. More importantly, TASC is
    able to adjust for covariates to further eliminate confounding that
    may originate from cell size and cell cycle differences.

-   [UNCURL](https://github.com/yjzhang/uncurl_python) - \[Python\] -
    Unsupervised and semi-supervised sampling effect removal for
    single-cell RNA-seq data.

**Cell projection and unimodal integration**

-   [Monet](https://github.com/flo-compbio/monet) - \[python\] - A
    package for analyzing and integrating scRNA-Seq data using PCA-based
    latent spaces.

**Simulation\*\*\***

-   [dropsim](https://github.com/marchinilab/dropsim) - \[R\] -
    Simulating droplet based scRNA-seq data.

-   [powsimR](https://bvieth.github.io/powsimR/) - \[R\] - Power
    analysis is essential to optimize the design of RNA-seq experiments
    and to assess and compare the power to detect differentially
    expressed
    genes. [PowsimR](https://academic.oup.com/bioinformatics/article/33/21/3486/3952669) is
    a flexible tool to simulate and evaluate differential expression
    from bulk and especially single-cell RNA-seq data making it suitable
    for a priori and posterior power analyses.

-   [splatter](http://bioconductor.org/packages/splatter/) - \[R\] -
    Splatter is a package for the simulation of single-cell RNA
    sequencing count data. It provides a simple interface for creating
    complex simulations that are reproducible and well-documented.

-   [symsim](https://github.com/YosefLab/SymSim) - \[R\] - SymSim
    (Synthetic model of multiple variability factors for Simulation) is
    an R package for simulation of single cell RNA-Seq data.

**Pseudotime and trajectory inference**

-   [CoSpar](https://cospar.readthedocs.io/) - \[python\] - CoSpar is a
    toolkit for dynamic inference by integrating state and lineage
    information. It gains superior robustness and accuracy by exploiting
    both the local coherence and sparsity of differentiation
    transitions, i.e., neighboring initial states share similar yet
    sparse fate outcomes. When only state information is available,
    CoSpar also improves upon existing dynamic inference methods by
    imposing sparsity and coherence.

-   [ECLAIR](https://github.com/GGiecold/ECLAIR) - \[python\] - ECLAIR
    stands for Ensemble Clustering for Lineage Analysis, Inference and
    Robustness. Robust and scalable inference of cell lineages from gene
    expression data.

-   [MERLoT](https://github.com/soedinglab/merlot) - \[R/python\] -
    Reconstructing complex lineage trees from scRNA-seq data using
    MERLoT.

-   [ouijaflow](http://www.github.com/kieranrcampbell/ouijaflow) -
    \[python\] - [A descriptive marker gene approach to single-cell
    pseudotime inference](https://doi.org/10.1101/060442)

-   [Palantir](https://github.com/dpeerlab/Palantir/) - \[Python\]
    - [Characterization of cell fate probabilities in single-cell data
    with Palantir](https://www.nature.com/articles/s41587-019-0068-4)

-   [SCDIFF](https://github.com/phoenixding/scdiff) - \[Python,
    JavaScript\] - SCDIFF is a single-cell trajectory inference method
    with interactive visualizations powered by D3.js. SCDIFF utilized
    the TF regulatory information to mitigate the impact of enormous
    single-cell RNA-seq noise (such as drop-out). With the TF regulatory
    information, SCDIFF is also able to predict the TFs (and their
    activation time), which drive the cells to different cell fates.
    Such predictive power has been [experimentally
    validated](https://genome.cshlp.org/content/28/3/383).

-   [SCIMITAR](https://github.com/dimenwarper/scimitar) - \[Python\] -
    Single Cell Inference of Morphing Trajectories and their Associated
    Regulation module (SCIMITAR) is a method for inferring biological
    properties from a pseudotemporal ordering. It can also be used to
    obtain progression-associated genes that vary along the trajectory,
    and genes that change their correlation structure over the
    trajectory; progression co-associated genes.

-   [scVelo](https://github.com/theislab/scvelo) - \[Python\] - scVelo
    is a scalable toolkit for RNA velocity analysis in single cells. It
    generalizes the concept of RNA velocity by relaxing previously made
    assumptions with a dynamical model. It allows to identify putative
    driver genes, infer a latent time, estimate reaction rates of
    transcription, splicing and degradation, and detect competing
    kinetics.

-   [TopSLAM](https://github.com/mzwiessele/topslam) - \[python\] -
    Extracting and using probabilistic Waddington\'s landscape
    recreation from single cell gene expression measurements

-   [VELOCYTO](http://velocyto.org/) - \[Python, R\] - Estimating RNA
    velocity in single cell RNA sequencing datasets.

**Cell type identification and classification**

-   [scANVI](https://github.com/YosefLab/scVI) - \[python\] -
    single-cell ANnotation using Variational Inference (scANVI) is a
    semi-supervised variant of scVI designed to leverage any available
    cell state annotations --- for instance when only one data set in a
    cohort is annotated, or when only a few cells in a single data set
    can be labeled using marker genes. [Harmonization and Annotation of
    Single-cell Transcriptomics data with Deep Generative
    Models](https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full)

-   [DeepSort](https://github.com/ZJUFanLab/DeepSort) - \[python\] - A
    reference-free cell-type annotation tool for single-cell RNA-seq
    data using deep learning with a weighted graph neural network, which
    is learned based on the most comprehensive single-cell
    transcriptomics atlases involving 764,741 cells across 88 tissues of
    human and
    mouse. [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.05.13.094953v1)

-   [ImmClassifier](https://github.com/xliu-uth/ImmClassifier) -
    \[R,python,Docker\] - A cell type annotation algorithm that employs
    a knowledge-based approach to annotating cells based on their
    underlying ontology and multitudes of previously-published data. By
    encoding immune cell hierarchy in a neural network, ImmClassifier is
    able to identify fine-grained cell types with high accuracy. By
    running in Docker the tool is
    platform-agnostic. [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.23.002758v1)

-   [Celltypist](https://pypi.org/project/celltypist-dev/) -
    \[Python\] - Celltypist is an automated cell type annotation tool
    for scRNA-seq datasets on the basis of logistic regression
    classifiers optimized by the stochastic gradient descent algorithm.
    Celltypist provides several different models for predictions, with a
    current focus on immune sub-populations, in order to assist in the
    accurate classification of different cell types and subtypes.

-   [scPRINT](https://github.com/cantinilab/scPRINT) - \[python\] -
    scPRINT is pretrained on 50M cells to predict multiple cell labels
    de novo, from any single cell RNAseq profile. [scPRINT: pre-training
    on 50 million cells allows robust gene network
    predictions](https://www.biorxiv.org/content/10.1101/2024.07.29.605556v1)

**Doublet Identification**

-   [AMULET](https://github.com/UcarLab/AMULET) - \[shell, Python, R\] -
    A count based method for detecting multiplets from single nucleus
    ATAC-seq (snATAC-seq) data. [Genome
    Biology](https://doi.org/10.1186/s13059-021-02469-x)

-   [demuxlet](https://github.com/statgen/demuxlet) - \[shell\]
    - [Multiplexed droplet single-cell RNA-sequencing using natural
    genetic variation](https://www.nature.com/articles/nbt.4042)

-   [DoubletDetection](https://github.com/JonathanShor/DoubletDetection) -
    \[R, Python\] - A Python3 package to detect doublets (technical
    errors) in single-cell RNA-seq count matrices. An [R
    implementation](https://github.com/TomKellyGenetics/DoubletDetection) is
    in development.

-   [Scrublet](https://github.com/swolock/scrublet) - \[Python\] -
    Computational identification of cell doublets in single-cell
    transcriptomic
    data. [BioRxiv](https://www.biorxiv.org/content/early/2018/07/09/357368)

-   [solo](https://github.com/calico/solo) - \[Python\] - Doublet
    detection via semi-supervised deep learning.

**Cell subsampling**

-   [geosketch](https://github.com/brianhie/geosketch) - \[Python\] -
    Method to subsample massive scRNA-seq datasets while preserving rare
    cell states. Resulting "sketch" accelerates clustering,
    visualization, and integration
    analyses. [Paper](https://doi.org/10.1016/j.cels.2019.05.003)

**Feature (Gene) imputation**

-   [MAGIC](https://github.com/KrishnaswamyLab/MAGIC) - \[R, Python,
    MATLAB\] - Markov Affinity-based Graph Imputation of Cells (MAGIC).
    A diffusion-based imputation method reveals gene-gene interactions
    in single-cell RNA-sequencing data.
    On [BioRviv](http://www.biorxiv.org/content/early/2017/02/25/111591) and [published
    in Cell](https://doi.org/10.1016/j.cell.2018.05.061).

-   [scPRINT](https://github.com/cantinilab/scPRINT) - \[python\] -
    scPRINT is pretrained on 50M cells to denoise and perform zero
    imputation of any single cell RNAseq profile. [scPRINT: pre-training
    on 50 million cells allows robust gene network
    predictions](https://www.biorxiv.org/content/10.1101/2024.07.29.605556v1)

**Variant calling**

-   [cb_sniffer](https://github.com/sridnona/cb_sniffer) - \[python\] -
    Mutation barcode caller, calls mutant and ref barcodes from 10x
    single cell data.

-   [cerebra](https://github.com/czbiohub/cerebra) - \[python\] -
    Cerebra is a tool for high-throughput summarizing of vcf entries
    following traditional variant calling for a sequencing experiment.
    Helps to extract relevant mutation information from among tens of
    thousands of vcf lines.

-   [monovar](https://bitbucket.org/hamimzafar/monovar) - \[python\] -
    Monovar is a single nucleotide variant (SNV) detection and
    genotyping algorithm for single-cell DNA sequencing data. It takes a
    list of bam files as input and outputs a vcf file containing the
    detected SNVs.

-   [SSrGE](https://github.com/lanagarmire/SSrGE) - \[python\] - SSrGE
    is an approach to identify SNVs correlated with Gene Expression
    using multiple regularized linear regressions. It contains its own
    pipeline to infer SNVs from scRNA-seq reads and is able to identify
    and sort genes and SNVs for a given cell subgroup. Nature
    Communication (2017) [Using single nucleotide variations in
    single-cell RNA-seq to identify subpopulations and
    genotype-phenotype
    linkage](https://www.nature.com/articles/s41467-018-07170-5).

**Epigenomics**

-   [AtacWorks](https://github.com/clara-parabricks/AtacWorks) -
    \[python\] - AtacWorks is a deep learning tool to denoise and
    identify accessible chromatin regions from low-coverage, low cell
    count, or low-quality ATAC-seq data. AtacWorks can denoise signal
    and identify peaks from rare cellular subtypes in a mixed
    population. [Biorxiv](https://www.biorxiv.org/content/10.1101/829481v2)

-   [BIRD](https://github.com/WeiqiangZhou/BIRD) - \[C++/R\] - BIRD is a
    tool for predicting chromatin accessibility and inferring regulatory
    element activities in single cells using scRNA-seq. [Global
    prediction of chromatin accessibility using small-cell-number and
    single-cell
    RNA-seq](https://academic.oup.com/nar/article/47/19/e121/5552071)

-   [DeepCpg](https://github.com/cangermueller/deepcpg) - \[python\] -
    DeepCpG is a deep neural network for predicting the methylation
    state of CpG dinucleotides in multiple cells. It allows to
    accurately impute incomplete DNA methylation profiles, to discover
    predictive sequence motifs, and to quantify the effect of sequence
    mutations.

-   [EpiScanpy](https://github.com/colomemaria/epiScanpy) - \[python\] -
    EpiScanpy is the epigenomic extension of scRNA-seq analysis tool
    Scanpy. It analyses single-cell open chromatin (scATAC-seq) and
    single-cell DNA methylation (for example scBS-seq) data. [EpiScanpy:
    integrated single-cell epigenomic
    analysis](https://www.nature.com/articles/s41467-021-25131-3)

-   [SCALE](https://github.com/jsxlei/SCALE) - \[python\] - SCALE is a
    deeplearning tool combining GMM with VAE for single-cell ATAC-seq
    analysis (visualization, clustering, imputation, batch effect
    removal, downstream analysis for celltype-specific TFs). [SCALE
    method for single-cell ATAC-seq analysis via latent feature
    extraction](https://www.nature.com/articles/s41467-019-12630-7)

-   [scbs](https://github.com/LKremer/scbs) - \[python\] - A command
    line tool for the analysis of Single-Cell Bisulfite-Sequencing data.
    scbs makes it easy to obtain a cell×region methylation matrix
    (≈count matrix) from raw single-cell methylation files and enables
    efficient storage, quality control and visualization. Furthermore,
    scbs allows you to scan the entire genome for variably or
    differentially methylated regions (VMRs or DMRs), and implements a
    novel approach for more accurate quantification of methylation in
    genomic
    intervals. [bioRχiv](https://doi.org/10.1101/2022.06.15.496318)

-   [scE2G](https://github.com/EngreitzLab/scE2G) - \[python\] - Family
    of models to predict enhancer-gene regulation based on single-cell
    data. These models use features derived from either single-cell
    ATAC-seq or multiomicRNA and ATAC-seq data in supervised classifiers
    trained on a new harmonized CRISPR perturbation dataset including
    over 13,000 examples.

**Multi-assay data integration**

-   [CITE-seq-Count](https://github.com/Hoohm/CITE-seq-Count) -
    \[python\] Cite-seq-Count is a python package that deals with
    Cellular Indexing of Transcriptomes and Epitopes by Sequencing
    (CITE-seq) and cell hashing
    data. [CITE-seq](https://cite-seq.com/) is a multimodal single cell
    phenotyping method that allows for immunophenotyping of cells with a
    potentially limitless number of markers and unbiased transcriptome
    analysis using existing single-cell sequencing approaches.

-   [Cobolt](https://github.com/epurdom/cobolt) - \[python\] - Cobolt is
    a novel method that not only allows for analyzing the data from
    joint-modality platforms, but provides a coherent framework for the
    integration of multiple datasets measured on different
    modalities. [Cobolt: integrative analysis of multimodal single-cell
    sequencing
    data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02706-x)

-   [gimVI](https://github.com/YosefLab/scVI) - \[python\] - gimVI is a
    joint generative model for imputation of missing genes in spatial
    transcriptomics assay from unpaired scRNA-seq data. [A joint model
    of unpaired data from scRNA-seq and spatial transcriptomics for
    imputing missing gene expression
    measurements](https://arxiv.org/pdf/1905.02269.pdf)

-   [GLUE](https://github.com/gao-lab/GLUE) - \[python\] - GLUE
    (Graph-Linked Unified Embedding) is a deep learning method for
    unpaired single-cell multi-omics data integration and regulatory
    inference
    ([Paper](https://www.nature.com/articles/s41587-022-01284-4)).

-   [MATCHER](https://github.com/jw156605/MATCHER) - \[python\]
    - [MATCHER: An algorithm for integrating single cell transcriptomic
    and epigenomic data using manifold alignment. MATCHER takes multiple
    types of single cell measurements performed on distinct single cells
    and infers single cell multi-omic
    profiles.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1269-0)

-   [MultiVI](https://github.com/scverse/scvi-tools) - \[python\] -
    MultiVI is a probabilistic framework that leverages deep neural
    networks to jointly analyze scRNA, scATAC and multiomic (scRNA +
    scATAC) data. [MultiVI: deep generative model for the integration of
    multi-modal
    data](https://www.biorxiv.org/content/10.1101/2021.08.20.457057v1)

-   [MOFA](https://github.com/bioFAM/MOFA) - \[python, R\] - Multi‐Omics
    Factor Analysis, a framework for unsupervised integration of
    multi‐omics data sets. MOFA is a method for disentangling the
    different sources of heterogeneity in bulk and single-cell
    multi-omics data sets. It identifies the latent factors that drive
    unique and shared variability in the different assays. The factors
    can be used for visualisation, pseudotime reconstruction,
    imputation, among other
    functionalities. [Paper](http://msb.embopress.org/content/14/6/e8124)

-   [SCALEX](https://github.com/jsxlei/SCALEX) - \[python\] - SCALEX
    provides a VAE framework for integration of heterogeneous
    single-cell data by disentangling batch-invariant components from
    batch-related variations and projecting the batch-invariant
    components into a generalized, low-dimensional cell-embedding
    space. [Construction of continuously expandable single-cell atlases
    through integration of heterogeneous datasets in a generalized
    cell-embedding
    space](https://www.biorxiv.org/content/10.1101/2021.04.06.438536v1)

-   [scarf](https://scarf.readthedocs.io) - \[python\] - 🧣 Toolkit for
    highly memory efficient analysis of single-cell RNA-Seq, scATAC-Seq
    and CITE-Seq data. Analyze atlas scale datasets with millions of
    cells on
    laptop. [Preprint](https://www.biorxiv.org/content/10.1101/2021.05.02.441899v1)

-   [scDART](https://github.com/PeterZZQ/scDART) - \[python\] - scDART
    is a deep learning framework that integrates scRNA-seq and
    scATAC-seq data and learns cross-modalities relationships
    simultaneously. [scDART: integrating unmatched scRNA-seq and
    scATAC-seq data and learning cross-modality relationship
    simultaneously](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02706-x)

-   [SISUA](https://github.com/trungnt13/sisua) - \[python\] - In this
    study, we propose models based on the Bayesian generative approach,
    where protein quantification available as CITE-seq counts from the
    same cells are used to constrain the learning process, thus forming
    a semi-supervised model. The generative model is based on the deep
    variational autoencoder (VAE) neural network
    architecture. [bioRxiv](https://www.biorxiv.org/content/10.1101/631382v1)

-   [TotalVI](https://github.com/YosefLab/scVI) - \[python\] - Total
    Variational Inference (totalVI) is a coupled generative model and
    inference procedure for CITE-seq data. TotalVI deals with
    modelisation of the background noise of protein measurements,
    harmonization of multiple CITE-seq experiments and imputation of
    missing proteins. [A Joint Model of RNA Expression and Surface
    Protein Abundance in Single
    Cells](https://www.biorxiv.org/content/biorxiv/early/2019/10/07/791947.full)

**Rare cell detection**

-   [FiRE](https://github.com/princethewinner/FiRE) - \[python, R,
    C++\] - Finder of rare entities (FiRE) helps identify rare cell
    types in voluminous single-cell datasets. Design of FiRE is inspired
    by the observation that rareness estimation of a particular data
    point is the flip side of measuring the density around it. In
    principle, FiRE uses the Sketching technique, a variant of locality
    sensitive hashing, to assign rareness score to every
    cell. [Paper](https://www.nature.com/articles/s41467-018-07234-6)

**Single cell large model**

-   [geneformer](https://github.com/jkobject/geneformer) - \[Python\] a
    single-cell large model training on 30 million human single-cell
    transcriptomics, supporting batch integration, gene dosage
    sensitivity predictions, chromatin dynamics prediction, network
    dynamics prediction, etc. manuscript open access: [Transfer learning
    enables predictions in network
    biology](https://www.nature.com/articles/s41586-023-06139-9)

-   [scGPT](https://github.com/bowang-lab/scGPT) - \[Python\] a
    single-cell large model training on 33 million human single-cell
    transcriptomics, supporting single-cell annotation, batch
    integration, perturbation prediction manuscript open access: [scGPT:
    toward building a foundation model for single-cell multi-omics using
    generative AI](https://www.nature.com/articles/s41592-024-02201-0)

-   [scFoundation](https://github.com/biomap-research/scFoundation) -
    \[Python\] a single-cell large model training on 50 million human
    single-cell transcriptomics with 100 million parameters, supporting
    single-cell clustering, drug response prediction, perturbation
    prediction, single-cell annotation, gene module inference, etc.
    manuscript open access: [Large-scale foundation model on single-cell
    transcriptomics](https://www.nature.com/articles/s41592-024-02305-7)

-   [CellPLM](https://github.com/OmicsML/CellPLM?tab=readme-ov-file) -
    \[Python\] the first single-Cell Pre-trained Language Model that
    encodes cell-cell relations and it consistently outperforms existing
    pre-trained and non-pre-trained models in diverse downstream tasks,
    with 100x higher inference speed compared to existing pre-trained
    models, training on 9 million scRNA-seq cells and 2 million SRT
    cells. manuscript open access: [CellPLM: Pre-training of Cell
    Language Model Beyond Single
    Cells](https://openreview.net/pdf?id=BKXvPDekud)

**Other applications**

-   [BASIC](http://ttic.uchicago.edu/~aakhan/BASIC/) - \[python\] -
    BASIC is a semi-de novo assembly method to determine the full-length
    sequence of the BCR in single B cells from scRNA-seq data.

-   [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) - \[python, R,
    snakemake\] - An automatic data handling pipeline for
    drop-seq/scrb-seq data. It runs from raw fastq.gz data until the
    final count matrix with QC plots along the way.

-   [ffq](https://github.com/pachterlab/ffq) - \[python\] - Fetch run
    and metadata information for single-cell genomics datasets.

-   [gget](https://github.com/pachterlab/gget) - \[Python\] - gget is a
    free, open-source command-line tool and Python package that enables
    efficient querying of genomic databases. gget consists of a
    collection of separate but interoperable modules, each designed to
    facilitate one type of database querying in a single line of code.

-   [immunarch](https://github.com/immunomind/immunarch) - \[R\] - R
    Package for Fast and Painless Exploration of Single-cell and Bulk
    T-cell/Antibody Immune Repertoires

-   [SCope](https://github.com/aertslab/SCope) - \[python\] - SCope is a
    fast visualization tool for large-scale and high dimensional
    scRNA-seq datasets.
    Publication [here](https://doi.org/10.1016/j.cell.2018.05.057).

-   [scTE](https://github.com/JiekaiLab/scTE) - \[python\] - Quantifying
    transposable element expression from single-cell sequencing data.

-   [SCIFER](https://github.com/hallvaaw/SCIFER) - \[shell, python\] -
    Approach for analysis of LINE-1 mRNA expression in single cells at a
    single locus resolution.

-   [Sinto](https://github.com/timoast/sinto) - \[python\] - A toolkit
    for working with aligned single-cell reads. Includes functions to
    split BAM files by cell barcode, add cell barcodes as read tags,
    move cell barcode information to read groups, and create a
    scATAC-seq fragment file from a BAM file.

-   [sircel](https://github.com/pachterlab/Sircel) - \[python\] - sircel
    (pronounced \"circle\") separates reads in a fastq file based on
    barcode sequences that occur at known positions of reads. This is an
    essential first step in analyzing single-cell genomics data from
    experiments such as Drop-Seq. Barcode sequences often contain
    deletion and/or mismatch errors that arise during barcode synthesis
    and sequencing, and we have designed our barcode recovery approach
    with these issues in mind. In addition to identifying barcodes in an
    unbiased manner, sircel also quantifies their
    abundances. [doi](https://doi.org/10.1101/136242)

-   [Snakemake single-cell-rna-seq
    workflow](https://github.com/snakemake-workflows/single-cell-rna-seq) -
    \[python, R, snakemake\] - An automated pipeline for single cell
    RNA-seq analysis.

-   [Wishbone](https://github.com/ManuSetty/wishbone) - \[python\]
    - [Wishbone is an algorithm to identify bifurcating developmental
    trajectories from single cell data. Wishbone can applied to both
    single cell RNA-seq and mass cytometry
    datasets.](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3569.html)

**Spatial transcriptomics**

-   [cell2location](https://github.com/BayraktarLab/cell2location) -
    \[Python\] A Bayesian model that perform spatial deconvolution in
    SRT data and create cellular maps of diverse tissues based on
    negative binomial distribution. manuscript open
    access: [Cell2location maps fine-grained cell types in spatial
    transcriptomics](https://www.nature.com/articles/s41587-021-01139-4)

-   [DestVI](https://github.com/scverse/scvi-tools) - \[Python\] A
    spatial deconvolution method designed with single-cell Latent
    Variable Model, scLVM(variational auto-encoder architecture).
    manuscript open access: [DestVI identifies continuums of cell types
    in spatial transcriptomics
    data](https://www.nature.com/articles/s41587-022-01272-8)

-   [DSTG](https://github.com/Su-informatics-lab/DSTG) - \[Python\] A
    spatial deconvolution method designed with graph-based convolutional
    networks. manuscript open access: [DSTG: deconvoluting spatial
    transcriptomics data through graph-based artificial
    intelligence](https://academic.oup.com/bib/article/22/5/bbaa414/6105942)

-   [Merfishtools](https://merfishtools.github.io/) - \[Python\] -
    MERFISHtools implement a Bayesian framework for accurately
    predicting gene or transcript expression from MERFISH data.

-   [NMFreg](https://github.com/tudaga/NMFreg_tutorial) - \[Python\] -
    The method is proposed
    in [Slide-seq](https://science.sciencemag.org/content/363/6434/1463) paper
    and reconstructs expression of each Slide-seq bead as a weighted
    combination of metagene factors, each corresponding to the
    expression signature of an individual cell type, defined from
    scRNA-seq.

-   [PASTE](https://github.com/raphael-group/paste) - \[Python\] A
    spatial alignment tool for aligning homogeneous spatial
    transcriptomic slices based on optimal transport and euclidean
    distance. manuscript open access: [Alignment and integration of
    spatial transcriptomics
    data](https://www.nature.com/articles/s41592-022-01459-6)

-   [PASTE2](https://github.com/raphael-group/paste2) - \[Python\] A
    spatial alignment tool for aligning homogeneous spatial
    transcriptomic slices based on the partial extension of the Fused
    Gromov-Wasserstein optimal transport. manuscript open
    access: [PASTE2: Partial Alignment of Multi-slice Spatially Resolved
    Transcriptomics
    Data](https://pmc.ncbi.nlm.nih.gov/articles/PMC9881963/pdf/nihpp-2023.01.08.523162v1.pdf)

-   [SLAT](https://github.com/gao-lab/SLAT) - \[Python\] SLAT is to
    align both homogeneous and heterogeneous (the first work) single
    cell spatial omics data by employing a graph alignment framework
    consists of LGCN and adversarial discriminator. manuscript open
    access: [Spatial-linked alignment tool (SLAT) for aligning
    heterogenous
    slices](https://www.nature.com/articles/s41467-023-43105-5)

-   [SpaGCN](https://github.com/jianhuupenn/SpaGCN) - \[Python\] SpaGCN
    is a graph convolutional network to integrate gene expression and
    histology to identify spatial domains and spatially variable genes.
    manuscript open access: [SpaGCN: Integrating gene expression,
    spatial location and histology to identify spatial domains and
    spatially variable genes by graph convolutional
    network](https://www.nature.com/articles/s41592-021-01255-8)

-   [SpatialDe](https://github.com/Teichlab/SpatialDE) - \[Python\]
    - [SpatialDE](https://www.nature.com/articles/nmeth.4636) is a
    statistical test to identify genes with spatial patterns of
    expression variation from multiplexed imaging or spatial
    RNA-sequencing data.

-   [SpatialPrompt](https://github.com/swainasish/SpatialPrompt) -
    \[Python\] [SpatialPrompt](https://www.nature.com/articles/s42003-024-06349-5) is
    a spot deconvolution and domain identification method for spatially
    resolved transcriptomics datasets. Main advantage of this tool is,
    it is highly scalable for large datasets.

-   [Splotch](https://github.com/tare/Splotch) - \[Python\] Splotch is a
    hierarchical generative probabilistic model for analyzing Spatial
    Transcriptomics data.

-   [squidpy](https://github.com/scverse/squidpy) - \[Python\] - Squidpy
    is a Python package for the analysis and visualization of spatial
    molecular data. It provides tools to process, analyze, and visualize
    spatial transcriptomics data, including spatially resolved
    transcriptomics and spatial proteomics. [Squidpy: a scalable
    framework for spatial omics data
    analysis](https://www.nature.com/articles/s41592-021-01264-7).

-   [Starspace](https://starspace.readthedocs.io/en/latest/) -
    \[Python\] - Defines a schema for gene or protein expression data
    containing spatially localized information. Converts data from a
    variety of assay types, including Spatial Transcriptomics, CODEX,
    In-situ Sequencing, MERFISH, osmFISH, and starMAP. Demonstrates how
    to visualize and interact with these data using common analysis
    packages, and convert the formats into loom and anndata objects, for
    downstream analysis in R and Python.

-   [STAGATE](https://github.com/QIFEIDKN/STAGATE) - \[Python\] STAGATE
    is designed for spatial clustering and denoising expressions of
    spatial resolved transcriptomics (ST) data by learning
    low-dimensional latent embeddings with both spatial information and
    gene expressions via a graph attention auto-encoder(GATE).
    manuscript open access: [Deciphering spatial domains from spatially
    resolved transcriptomics with an adaptive graph attention
    auto-encoder](https://www.nature.com/articles/s41467-022-29439-6)

-   [STAligner](https://github.com/zhoux85/STAligner) - \[Python\]
    STAligner is designed for alignment and integration of spatially
    resolved transcriptomics data, it employs a graph attention
    auto-encoder neural network(GATE) to extract spatially aware
    embedding and introduces the triplet loss to update the spot
    embedding to reduce the distance from the anchor to positive spot.
    manuscript open access: [STAligner enables the integration and
    alignment of multiple spatial transcriptomics
    datasets](https://www.nature.com/articles/s43588-023-00543-x)

-   [Tangram](https://github.com/broadinstitute/Tangram) - \[Python\]
    Tangram is used to map single-cell (or single-nucleus) gene
    expression data onto spatial gene expression data designed with
    optimizing a specially designed mapping objective loss. manuscript
    open access: [Deep learning and alignment of spatially resolved
    single-cell transcriptomes with
    Tangram](https://www.nature.com/articles/s41592-021-01264-7)

**Tutorials and workflows**

-   [Analysis of single cell RNA-seq
    data](https://github.com/cellgeni/scRNA.seq.course) - \[R and
    Python\] -
    The [course](https://www.singlecellcourse.org/index.html) is taught
    through the University of Cambridge Bioinformatics training unit,
    but the material found on these pages is meant to be used for anyone
    interested in learning about computational analysis of scRNA-seq
    data.

-   [Aaron Lun\'s Single Cell workflow on
    Bioconductor](http://bioconductor.org/help/workflows/simpleSingleCell/) -
    \[R\] - This article describes a computational workflow for basic
    analysis of scRNA-seq data using software packages from the
    open-source Bioconductor project.

-   [ATAC-Seq Pipeline](https://github.com/tobiasrausch/ATACseq) -
    \[Shell and R\] - [Chromatin accessibility landscape of pediatric
    T-lymphoblastic leukemia and human T-cell
    precursors](https://www.embopress.org/doi/full/10.15252/emmm.202012104).

-   [Bioconductor2016 Single-cell-RNA-sequencing workshop by Sandrine
    Dudoit lab](https://github.com/drisso/bioc2016singlecell) - \[R\] -
    SCONE, clusterExperiment, and slingshot tutorial.

-   [BiomedCentral Single Cell Omics
    collectin](http://www.biomedcentral.com/collections/singlecellomics) -
    collection of papers describing techniques for single-cell analysis
    and protocols.

-   [Clustering 3K PBMCs with Scanpy in
    Galaxy](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-scanpy-pbmc3k/tutorial.html) -
    Galaxy Training Material.

-   [**CRUK CI** Introduction to single-cell RNA-seq data
    analysis](https://github.com/bioinformatics-core-shared-training/UnivCambridge_ScRnaSeq_Nov2021/) \| [website](https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/)

-   [CSHL Single Cell Analysis -
    Bioinformatics](https://github.com/YeoLab/single-cell-bioinformatics/) course
    materials - Uses Shalek 2013 and Macaulay 2016 datasets to teach
    machine learning to biologists

-   [**Dana-Farber Cancer Institute** Trajectory inference across
    conditions: differential expression and differential
    progression](https://github.com/kstreet13/bioc2020trajectories) \| [website](https://kstreet13.github.io/bioc2020trajectories/)

-   [**Dan Beiting** DIY
    Transcriptomics](https://github.com/DIYtranscriptomics/DIYtranscriptomics.github.io) \| [website](https://diytranscriptomics.com/) -
    A hybrid course covering best practices for bulk and single cell
    RNA-seq data analysis, with a primary focus on empowering students
    to be independent in the use of lightweight and open-source software
    and the R/bioconductor environment.

-   [**EBI** Single cell RNA-seq
    tutorial](https://gitlab.com/mperalc/scRNA-seq_workshop_2021) \| [website](https://mperalc.gitlab.io/scRNA-seq_workshop_2021/)

-   [**ELIXIR EXCELERATE** Single RNA-seq data analysis with
    R](https://github.com/NBISweden/excelerate-scRNAseq) \| [website](https://nbisweden.github.io/excelerate-scRNAseq/)

-   [Festival of Genomics California Single Cell
    Workshop](https://kdkorthauer.github.io/FestivalWorkshopVignettes/) -
    \[R\] - Explores basic workflow from exploratory data analysis to
    normalization and downstream analyses using a dataset of 1679 cells
    from the Allen Brain Atlas.

-   [Gilad Lab Single Cell Data
    Exploration](http://jdblischak.github.io/singleCellSeq/analysis/) -
    R-based exploration of single cell sequence data. Lots of
    experimentation.

-   [GPU accelerated single-cell analysis using
    RAPIDS](https://github.com/clara-parabricks/rapids-single-cell-examples) -
    NVIDIA tutorials on using RAPIDS (<https://www.rapids.ai/>) to
    accelerate single-cell analysis on GPUs.

-   [**Harvard Chan Bioinformatics Core** Single-cell RNA-seq data
    analysis
    workshop](https://github.com/hbctraining/scRNA-seq_online) \| [website](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html)

-   [Harvard STEM Cell Institute Single Cell Workshop
    2015](http://hms-dbmi.github.io/scw/) - workshop on common
    computational analysis techniques for scRNA-seq data from
    differential expression to subpopulation identification and network
    analysis. [See course description for more
    information](http://scholar.harvard.edu/jeanfan/classes/single-cell-workshop-2015)

-   [Hemberg Lab scRNA-seq course
    materials](http://hemberg-lab.github.io/scRNA.seq.course/index.html)

-   [kallistobustools](https://github.com/pachterlab/kallistobustools) -
    kallisto \| bustools workflow for pre-processing single-cell RNA-seq
    data.

-   [**NBIS** Single cell RNA-seq analysis
    workshop](https://github.com/nbisweden/workshop-scRNAseq) \| [website](https://nbisweden.github.io/workshop-scRNAseq/)

-   [**MGC/BioSB** Course - Single Cell
    Analysis](https://github.com/LeidenCBC/MGC-BioSB-SingleCellAnalysis2021)

-   [nf-core/scrnaseq](https://nf-co.re/scrnaseq) - nf-core/scrnaseq is
    a bioinformatics best-practice analysis pipeline for processing 10x
    Genomics single-cell RNA-seq data. The pipeline is built using
    Nextflow, a workflow tool to run tasks across multiple compute
    infrastructures in a very portable manner. It uses
    Docker/Singularity containers making installation trivial and
    results highly reproducible.

-   [nf-core/scflow](https://nf-co.re/scflow) - nf-core/scflow is a
    bioinformatics pipeline for scalable, reproducible, best-practice
    analyses of single-cell/nuclei RNA-sequencing data. The pipeline is
    built using Nextflow, a workflow tool to run tasks across multiple
    compute infrastructures in a very portable manner.

-   [nf-core/scnanoseq](https://nf-co.re/scnanoseq) - nf-core/scnanoseq
    is a Nextflow analysis pipeline for processing 10X Genomics
    single-cell/nuclei RNA-seq data derived from the Oxford Nanopore
    long-read sequencer. The pipeline has been designed to be scalable
    to large datasets (including PromethION data), portable and
    reproducible.

-   [Orchestrating Single-Cell Analysis with
    Bioconductor](https://bioconductor.org/books/release/OSCA/) -
    \[R\] - This blogdown book describes a comprehensive and
    reproducible workflow for the analysis of single-cell RNA-sequencing
    data.

-   [Pre-processing of 10X Single-Cell RNA Datasets in
    Galaxy](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-preprocessing-tenx/tutorial.html) -
    Galaxy Training Material.

-   [Theis Lab Single Cell
    Tutorial](https://github.com/theislab/single-cell-tutorial) - The
    main part of this repository is a case study where the
    best-practices established in the manuscript are applied to a mouse
    intestinal epithelium regions dataset from Haber et al., Nature
    551 (2018) available from the GEO under GSE92332.

-   [Using Seurat (v1.2) for unsupervised clustering and biomarker
    discovery](http://www.satijalab.org/seurat/get_started_v1_2.html) -
    301 single cells across diverse tissues from (Pollen et al., Nature
    Biotechnology, 2014). Original tutorial using Seurat 1.2

-   [Using Seurat (v1.2) for spatial inference in single-cell
    data](http://www.satijalab.org/seurat/get_started_v1_2.html) - 851
    single cells from zebrafish embryogenesis (Satija\*, Farrell\* et
    al., Nature Biotechnology, 2015). Original tutorial using Seurat 1.2

-   [Seurat (v3.0) - Guided Clustering
    Tutorial](https://satijalab.org/seurat/v3.0/immune_alignment.html) -
    new tutorial using Seurat 3.0

-   [**SIB** Single-cell
    Transcriptomics](https://github.com/sib-swiss/single-cell-training/) \| [website](https://sib-swiss.github.io/single-cell-training/latest/)

-   [**SIB NBIS/SciLifeLab** Advanced topics in Single Cell
    Omics](https://github.com/NBISweden/single-cell_sib_scilifelab_2021) \| [website](https://nbisweden.github.io/single-cell_sib_scilifelab_2021/)

-   [**SIB** Advanced topics in single-cell
    transcriptomics](https://github.com/fmicompbio/adv_scrnaseq_2020)

-   [**WEHI** Single cell RNA-seq analysis
    workshop](https://github.com/yunshun/SingleCellWorkshop/) \| [website](https://yunshun.github.io/SingleCellWorkshop/)

-   [**Wellcome Sanger Institute** Analysis of single cell RNA-seq
    data](https://github.com/hemberg-lab/scRNA.seq.course) \| [website](https://www.singlecellcourse.org/)

**Web portals, apps, and databases**

**Web portals and databases**

-   [10X Genomics
    datasets](https://support.10xgenomics.com/single-cell/datasets) -
    10x genomics public datasets, including 1.3M cell mouse brain
    dataset.

-   [ASAP](http://asap.epfl.ch/) - Automated Single-cell Analysis
    Pipeline (deposited
    in [BioRXiv](http://biorxiv.org/content/early/2016/12/22/096222) on
    December 22, 2016).

-   [cellBrowser](https://github.com/maximilianh/cellBrowser) -
    \[Python, Javascript\] Python pipeline and Javascript scatter plot
    library for single-cell datasets. [Demo](https://cells.ucsc.edu/)

-   [CellView](https://github.com/mohanbolisetty/CellView) - CellView is
    an R Shiny web application that allows knowledge-based and
    hypothesis-driven exploration of processed single cell
    transcriptomic
    data. [ref](https://www.biorxiv.org/content/early/2017/04/04/123810).

-   [Cell_BLAST](http://cblast.gao-lab.org/) - A Web portal powered by
    Cell_BLAST (scRNA-seq querying tool) and ACA (scRNA-seq database).

-   [CELLxGENE](https://cellxgene.cziscience.com/) - CELLxGENE is a
    suite of tools that help scientists to find, download, explore,
    analyze, annotate, and publish single cell datasets. It includes
    several powerful tools with various features to help you to engage
    with single cell data.

-   [conquer](http://imlspenticton.uzh.ch:3838/conquer/) - A repository
    of consistently processed, analysis-ready single-cell RNA-seq data
    sets.

-   [Curated Database of single-cell
    studies](http://www.nxn.se/single-cell-studies/) - Available as
    a tsv download. Over 500 single cell transcriptomics studies have
    been published to date. Many of these have data available, but the
    links between data, study, and systems studied can be hard to
    identify through literature search. This manuscript describes a
    nearly exhaustive and manually curated database of single cell
    transcriptomics studies with descriptions of what kind of data and
    what biological systems have been
    studied. [bioRxiv](https://doi.org/10.1101/742304).

-   [D^3^E](http://www.sanger.ac.uk/sanger/GeneRegulation_D3E/) -
    Discrete Distributional Differential Expression (D^3^E) is a tool
    for identifying differentially-expressed genes, based on single-cell
    RNA-seq data.

-   [dseqr](https://docs.dseqr.com/) - Dseqr runs end-to-end
    multi-sample single-cell and bulk RNA-seq analyses using a user
    friendly web app built around best practices from the OSCA handbook.
    Features include pseudobulk differential expression analysis,
    automated cluster annotation, reference mapping with Azimuth, Gene
    Ontology analysis, and drug connectivity mapping. Projects can
    either be analysed online or locally using the [dseqr R
    package](https://github.com/hms-dbmi/dseqr).

-   [EBI Single Cell Expression
    Atlas](https://www.ebi.ac.uk/gxa/sc/home) - The Single Cell
    Expression Atlas contains uniformly re-analysed single cell
    expression data across different species and provides interactive
    visualizations to explore that data.

-   [Galaxy Single Cell Omics
    Workbench](https://singlecell.usegalaxy.eu/) - dedicated Galaxy
    server for analyzing single cell data.

-   [IRIS3](https://bmbl.bmi.osumc.edu/iris3/) - IRIS3 (integrated
    cell-type-specific regulon inference server from single-cell
    RNA-Seq) is an easy-to-use server empowered by over 20
    functionalities to support comprehensive interpretations and
    graphical visualizations of identified cell-type-specific regulons.

-   [JingleBells](http://jinglebells.bgu.ac.il/) - A repository of
    standardized single cell RNA-Seq datasets for analysis and
    visualization in IGV at the single cell level. Currently focused on
    immune cells (<http://www.jimmunol.org/content/198/9/3375.long>).

-   [SCPortalen](http://single-cell.clst.riken.jp/) - SCPortalen: human
    and mouse single-cell centric
    database. [ref](https://doi.org/10.1093/nar/gkx949)

-   [scRNA.seq.datasets](https://hemberg-lab.github.io/scRNA.seq.datasets) -
    Collection of public scRNA-Seq datasets used by [Hemberg
    Lab](http://www.sanger.ac.uk/science/groups/hemberg-group)

-   [scRNASeqDB](https://bioinfo.uth.edu/scrnaseqdb/) - A database
    aggregating human single-cell RNA-seq
    datasets. [ref](http://biorxiv.org/content/early/2017/01/31/104810)

-   [Single Cell
    Portal](https://portals.broadinstitute.org/single_cell) - The
    Single-Cell Portal was developed to facilitate open data and open
    science in Single-cell Genomics. The portal currently focuses on
    sharing scientific results interactively, and sharing associated
    datasets.

-   [Single-Cell Tumor Immune Atlas
    project](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/Tumor-Immune-Cell-Atlas) -
    \[R\] - We generated a [single-cell tumor immune
    atlas](https://singlecellgenomics-cnag-crg.shinyapps.io/TICA/),
    jointly analyzing \>500,000 cells from 217 patients and 13 cancer
    types, providing the basis for a patient stratification based on
    immune cell compositions.

-   [V-SVA](https://github.com/nlawlor/V-SVA) - An R Shiny application
    for detecting and annotating hidden sources of variation in single
    cell RNA-seq data.

-   [WOT](https://broadinstitute.github.io/wot/) - Waddington Optimal
    Transport (wot) uses time-course data to infer how the probability
    distribution of cells in gene-expression space evolves over time, by
    using the mathematical approach of optimal transport.

**Interactive visualization and analysis**

-   [Cellar](https://github.com/euxhenh/cellar) - \[Python\] - Cellar is
    an easy to use, interactive, and comprehensive software tool for the
    assignment of cell types in single-cell studies. It supports
    preprocessing, dimensionality reduction, clustering, differential
    expression & enrichment analysis, spatial transcriptomics, label
    transfer, semi-supervised clustering and more. A live web server
    running Cellar is
    available [here](https://cellar.cmu.hubmapconsortium.org/app/cellar). [NatureComm](https://www.nature.com/articles/s41467-022-29744-0).

-   [cellBrowser](https://github.com/maximilianh/cellBrowser) -
    \[Python, Javascript\] Python pipeline and Javascript scatter plot
    library for single-cell datasets. [Demo](https://cells.ucsc.edu/)

-   [Cirrocumulus](https://github.com/klarman-cell-observatory/cirrocumulus) -
    Cirrocumulus is an interactive visualization tool for large-scale
    single-cell genomics (e.g. sc/snRNA-seq, spatial) data.

-   [CReSCENT](http://crescent.cloud/) - \[R, Javascript, Python\] -
    CReSCENT: CanceR Single Cell ExpressioN Toolkit ([Mohanraj et al.
    2020](https://doi.org/10.1093/nar/gkaa437)), is an intuitive and
    scalable web portal incorporating a containerized pipeline execution
    engine for standardized analysis of cancer scRNA-seq data and
    associated metadata. CReSCENT uses public data sets and
    preconfigured pipelines that are accessible to computational biology
    non-experts and are user-editable to allow optimization, comparison,
    and reanalysis for specific experiments. Users can also upload their
    own scRNA-seq data for analysis and results can be kept private or
    shared with other users.

-   [FASTGenomics](https://fastgenomics.org) - \[Python, R\] -
    FASTGenomics is an online platform to share single-cell RNA
    sequencing data and analyses using reproducible workflows. Gene
    expression data can be shared meeting European data protection
    standards (GDPR). FASTGenomics enables the user to upload their own
    data and generate customized and reproducible workflows for the
    exploration and analysis of gene expression data ([Scholz et al.
    2018](https://doi.org/10.1101/272476)). Follow us
    on [Twitter](https://twitter.com/FastGenomics).

-   [iS-CellR](https://github.com/immcore/iS-CellR) - iS-CellR
    (Interactive platform for Single-cell RNAseq) is a web-based Shiny
    app that integrates the Seurat package with Shiny\'s reactive
    programming framework to provide comprhensive analysis and
    interactive visualization of single-cell RNAseq
    data. [Paper](https://doi.org/10.1093/bioinformatics/bty517)

-   [PIVOT](https://github.com/qinzhu/PIVOT) - Platform for Interactive
    analysis and Visualization Of Transcriptomics
    data. [ref](https://doi.org/10.1186/s12859-017-1994-0)

-   [scRNAseqApp](https://bioconductor.org/packages/scRNAseqApp/) - The
    scRNAseqApp is a Shiny app package designed for interactive
    visualization of single-cell data. It is an enhanced version derived
    from the [ShinyCell](https://github.com/SGDDNB/ShinyCell),
    repackaged to accommodate multiple datasets. The app enables users
    to visualize data containing various types of information
    simultaneously, facilitating comprehensive analysis. Additionally,
    it includes a user management system to regulate database
    accessibility for different users.

-   [singleCellTK](https://bioconductor.org/packages/singleCellTK) - The
    singleCellTK is an R/Shiny package and GUI for analyzing and
    visualizing scRNA-Seq through a web interface. Analysis modules
    include data summary and filtering, dimensionality reduction and
    clustering, batch correction, differential expression analysis,
    pathway activity analysis, and power analysis. [Comprehensive
    generation, visualization, and reporting of quality control metrics
    for single-cell RNA sequencing
    data](https://doi.org/10.1038/s41467-022-29212-9).

-   [spatialGE-web](https://spatialge.moffitt.org) - A web application
    providing point-and-click access to the methods in the [spatialGE R
    package](https://github.com/FridleyLab/spatialGE) and other tools
    (including STdeconvolve, InSituType, SpaGCN). The web application
    requires no coding experience. User accounts can be created to
    safely keep samples and results organized within projects.

-   [STREAM](http://stream.pinellolab.org/) - STREAM is an interactive
    computational pipeline for reconstructing complex celluar
    developmental trajectories from sc-qPCR, scRNA-seq or scATAC-seq
    data. [preprint](https://doi.org/10.1101/302554).

-   [Vitessce](http://vitessce.io) - \[JavaScript, Python, R\] - A
    framework for integrative visualization of multi-modal single-cell
    data, supporting microscopy images, cell segmentations,
    dimensionality reduction scatterplots, gene expression heatmaps, and
    genome browser tracks. Vitessce is packaged as a [React
    component](http://github.com/vitessce/vitessce), [Jupyter
    widget](http://github.com/vitessce/vitessce-python), and [R
    htmlwidget](http://github.com/vitessce/vitessce-r).

**Big data approach overview**

-   [Single-cell Transcriptome Study as Big
    Data](http://www.sciencedirect.com/science/article/pii/S1672022916000437)

**Experimental design**

-   [Design and computational analysis of single-cell RNA-sequencing
    experiments](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y)

-   [How to design a single-cell RNA-sequencing experiment: pitfalls,
    challenges and perspectives](https://doi.org/10.1093/bib/bby007)

-   [Sensei](https://kchen-lab.github.io/sensei/table_beta.html). [How
    to determine the number of patients required to ascertain a
    cell-type abundance change estimated from scRNA-seq
    experiment](https://www.biorxiv.org/content/10.1101/2020.05.31.126565v1).

**People**

Gender bias at conferences is a well known problem
(<http://www.sciencemag.org/careers/2015/07/countering-gender-bias-conferences>).
Creating a list of potential speakers can help mitigate this bias and a
community of people developing and maintaining helps to further
diversify this list beyond smaller networks.

**Female**

-   [Rhonda Bacher (University of Wisconsin-Madison,
    USA)](https://twitter.com/rbacher)

-   [Barbara Di Camillo (Information Engineering Department, University
    of Padova, Italy](https://www.dei.unipd.it/~dicamill/)

-   [Jinmiao Chen (Singapore Immunology Network, A\*STAR,
    Singapore)](https://www.a-star.edu.sg/sign/PEOPLE/Principal-Investigators/Investigator-Details?givenName=Jinmiao&lastName=CHEN)

-   [Jean Fan (Johns Hopkins University, USA)](https://jef.works/team)

-   [Brooke Fridley (Children\'s Mercy Hospital,
    USA)](https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/labs-and-research-teams/fridley-lab/)

-   [Sandrine Dudoit (UC Berkeley,
    USA)](http://www.stat.berkeley.edu/~sandrine/)

-   [Lana X. Garmire, (University of Hawaii Cancer Center,
    USA)](http://garmiregroup.org/)

-   [Laleh Haghverdi (EMBL,
    Germany)](https://www.embl.de/research/units/genome_biology/huber/members/index.php?s_personId=CP-60025160)

-   [Stephanie Hicks (Johns Hopkins Bloomberg School of Public Health,
    USA)](http://www.stephaniehicks.com/)

-   [Christina Kendziorski (University of Wisconsin--Madison,
    USA)](https://www.biostat.wisc.edu/~kendzior/)

-   [Keegan Korthauer (Dana Farber Cancer Institute,
    USA)](http://kkorthauer.org/)

-   [Smita Krishnaswamy (Yale University)](https://krishnaswamylab.org)

-   [Ning Leng (Morgridge Institute for Research,
    USA)](https://www.biostat.wisc.edu/~ningleng/)

-   [Elisabetta Mereu (Centre for Genomic Regulation,
    Barcelona)](http://www.cnag.cat/elisabetta-mereu)

-   [Samantha Morris (Depts of Dev. Bio. and Genetics, Washington
    University, St. Louis)](http://morrislab.wustl.edu/)

-   [Alicia Oshlack (Murdoch Children\'s Research Institute,
    Australia)](https://www.mcri.edu.au/users/dr-alicia-oshlack)

-   [Dana Pe\'er (Memorial Sloan Kettering Cancer Center,
    USA)](https://www.mskcc.org/research/ski/labs/dana-pe-er)

-   [Emma Pierson (Stanford University,
    USA)](http://cs.stanford.edu/people/emmap1/)

-   [Aviv Regev (Broad Institute,
    USA)](https://www.broadinstitute.org/scientific-community/science/core-faculty-labs/regev-lab/regev-lab-home)

-   [Sarah Snelling (University of Oxford,
    UK)](https://www.ndorms.ox.ac.uk/research/research-groups/snelling)

-   [Charlotte Soneson (Institute of Molecular Life Sciences, University
    of Zurich)](https://csoneson.github.io/)

-   [Sarah Teichmann (Wellcome Trust Sanger Institute,
    UK)](http://www.teichlab.org/)

-   [Barbara Treutlein (ETH Zurich, CH)](https://bsse.ethz.ch/qdb)

-   [Catalina Vallejos (The Alan Turing Institute & UCL,
    UK)](https://sites.google.com/view/catalinavallejos)

-   [Sanja Vickovic (New York Genome Center & Columbia University,
    USA)](https://www.vickovictechinnovation.com/)

**Male**

-   [Stein Aerts (KU Leuven Center for Human Genetics,
    Belgium)](http://aertslab.org/)

-   [Ahmet Coskun (Georgia Tech,
    USA)](https://coskunlab.bme.gatech.edu/)

-   [Bart DePlancke (EPFL, School of Life sciences, Institute of
    Bioengineering, Switzerland)](https://deplanckelab.epfl.ch/)

-   [Raphael Gottardo (Fred Hutchinson Cancer Research Center,
    USA)](https://www.fredhutch.org/en/labs/profiles/gottardo-raphael.html)

-   [Chung Chau Hon (RIKEN Centre for Integrative Medical Sciences,
    Yokohama)](http://www.riken.jp/en/research/labs/ims/genom_inf/)

-   [Yanxiang Deng (University of Pennsylvania,
    USA)](https://yanxiangdenglab.org/)

-   [Martin Hemberg (Sanger Institute,
    UK)](http://www.sanger.ac.uk/science/groups/hemberg-group)

-   [John Hickey (Duke University,
    USA)](https://sites.google.com/view/john-w-hickey/)

-   [Holger Heyn (Centre for Genomic Regulation,
    Barcelona)](http://www.cnag.crg.eu/holger-heyn)

-   [Peter Kharchenko (Department of Biomedical Informatics, Harvard
    Medical School, USA)](http://pklab.med.harvard.edu/)

-   [Sten Linnarson (Karolinska Institutet,
    Sweden)](http://linnarssonlab.org/)

-   [Aaron Lun (Cancer Research UK,
    UK)](http://www.cruk.cam.ac.uk/users/aaron-lun)

-   [John Marioni (EBI, UK)](http://www.ebi.ac.uk/research/marioni)

-   [Davis McCarthy (EBI, UK)](https://sites.google.com/site/davismcc/)

-   [John Reid (MRC Biostatistics Unit, Cambridge University,
    UK)](http://johnreid.github.io/)

-   [Mark Robinson (Institute of Molecular Life Sciences, University of
    Zurich)](http://www.imls.uzh.ch/en/research/robinson.html)

-   [Yvan Saeys (Vlaams Instituut voor Biotechnologie, Ghent,
    Belgium)](http://www.vib.be/en/research/scientists/Pages/Yvan-Saeys-Lab.aspx)

-   [Rickard Sandberg (Karolinska Institutet,
    SE)](https://sandberglab.se/)

-   [Neville Sanjana (New York Genome Center &
    NYU)](http://sanjanalab.org/)

-   [Rahul Satija (New York Genome Center)](http://satijalab.org/)

-   [Peter Sims (Columbia University, Department of Systems
    Biology)](http://www.columbia.edu/~pas2182/index.php/home-top.html)

-   [Oliver Stegle (EBI, UK)](http://www.ebi.ac.uk/research/stegle)

-   [Fabian Theis (Institute of Computational Biology, Helmholtz Zentrum
    München)](https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/2494/index.html)

-   [Cole Trapnell (University of Washington, Department of Genome
    Sciences)](http://cole-trapnell-lab.github.io/)

-   [Itai Yanai (New York University, School of Medicine, Institute for
    Computational Medicine, USA)](https://yanailab.org/)
