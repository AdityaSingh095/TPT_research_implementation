# TPT-Based Molecular Reaction Graph Embeddings

Research implementation and extension of **“Clustering Molecular Energy Landscapes by Adaptive Network Embedding”** by Paula Mercurio and Di Liu (2024), adapted to a large molecular reaction network.

The project combines **Transition Path Theory (TPT)**, thermodynamically weighted random walks, and Skip-Gram embeddings to learn molecular representations that capture reaction-network dynamics. The original source-target TPT formulation is further extended to a **multi-pair TPT framework** for broader coverage of a highly disconnected reaction graph.

## Results

| Metric                           |              Result |
| -------------------------------- | ------------------: |
| Molecules                        |         **123,174** |
| Reaction edges                   |         **770,844** |
| TPT source-target pairs          |              **20** |
| Multi-pair TPT coverage          |           **54.2%** |
| Active TPT edges                 | **329,884 (42.8%)** |
| Random walks                     |       **1,231,740** |
| Embedding dimension              |             **128** |
| Embedded molecules               |  **123,174 (100%)** |
| Intra-category cosine similarity |          **0.4717** |
| Inter-category cosine similarity |          **0.3775** |
| Embedding separation             |         **+0.0942** |

The multi-pair formulation produced non-zero TPT committor values across roughly **66.8K molecules**, while random-walk sampling generated embeddings for the complete 123K-node graph.

## Method

### 1. Thermodynamic Reaction Graph

The molecular reaction network contains **123,174 molecular nodes and 770,844 directed reaction edges**.

Reaction edges are weighted using reaction free energy (`ΔG`), kinetic competition between outgoing reactions, and a robustness term favouring energetically viable transitions.

This converts the reaction network into a weighted Markov graph where transition probabilities reflect both energetic favourability and competing reaction routes.

### 2. Entropy-Based Weight Optimisation

The raw thermodynamic weights span many orders of magnitude, causing most transitions to become effectively invisible during random walks.

Temperature scaling is therefore applied to the edge weights.

The scaling parameter `τ` is optimised with **Optuna over 50 trials** by maximising the Shannon entropy of the resulting transition probabilities.

Best observed temperature:

```text
τ ≈ 19.99
```

After scaling:

```text
Dominant edges (> 1e-2)       : 276,456  (35.86%)
Visible edges (1e-8 – 1e-2)  : 494,388  (64.14%)
Edges below 1e-8             : 0
```

This prevents the random walk from collapsing onto only a small number of dominant reactions.

## Transition Path Theory

For a source molecule `A` and target molecule `B`, TPT computes the probability that a molecule lies along a reactive trajectory between them.

The implementation calculates:

* stationary distribution `π`
* forward committor `q+`
* reactive probability current
* effective current `f+`

The committor is obtained by solving a sparse linear system over the reaction graph.

GPU acceleration using **CuPy / cupyx sparse solvers** is used when available, with SciPy sparse solvers as the CPU fallback.

## Multi-Pair TPT Extension

Standard TPT considers one source-target pair. This becomes restrictive for the reaction graph because a large fraction of the network is split across disconnected or weakly connected regions.

The implementation therefore introduces a **20-pair TPT extension**.

### Pair Selection

The graph is first projected to an undirected weighted network and partitioned using **Louvain community detection**.

Candidate source-target pairs are then selected using forward/backward graph reachability so that each pair represents a connected reaction region.

For every selected pair:

1. Productive and sink regions are identified.
2. The forward committor is solved.
3. Reactive currents are calculated.
4. Effective currents are extracted.
5. Results are aggregated across pairs.

Final coverage:

```text
q+ > 0 nodes     : ~66.6K / 123.2K
TPT coverage     : ~54%
Active TPT edges : 329,884 / 770,844
```

The remaining nodes largely belong to reaction components that cannot participate in a valid source-target transition for the selected TPT pairs.

## TPT-Guided Random Walks

The aggregated TPT quantities are incorporated into the transition probabilities used for graph sampling.

Each node receives:

```text
10 walks/node
9 transitions/walk
```

producing:

```text
1,231,740 total random walks
```

The sampling combines the underlying thermodynamic transition weight with the aggregated effective TPT current.

This biases walks toward dynamically relevant reaction pathways while retaining graph-wide exploration.

## Molecular Embeddings

The sampled reaction trajectories are treated as sequences and used to train a **Skip-Gram Word2Vec model**.

Configuration:

```text
Embedding dimension : 128
Context window      : 5
Epochs              : 20
Minimum count       : 1
```

Result:

```text
Embedded molecules : 123,174 / 123,174
Embedding coverage : 100%
```

Each molecule is therefore represented by a **128-dimensional vector encoding its reaction-network context**.

Embeddings are exported both as a complete dictionary and as individual NumPy vectors for downstream use.

## Similarity Propagation

Cosine similarity is computed between molecular embeddings and propagated over edges carrying positive TPT effective current.

This enables retrieval of molecules related to a reference molecule through both:

* embedding-space similarity
* reaction-pathway connectivity

The implementation uses batched matrix operations for cosine computation rather than per-node similarity loops.

## Embedding Evaluation

Embedding structure is evaluated by separating molecules into:

```text
Sink       : q+ = 0
Passive    : q+ > 0, f+ = 0
Active     : q+ > 0, f+ > 0
```

Observed embedding statistics:

```text
Mean intra-category cosine similarity : 0.4717
Mean inter-category cosine similarity : 0.3775
Separation                            : +0.0942
```

The positive separation indicates that molecules with similar TPT behaviour are closer in the learned embedding space than molecules belonging to different pathway categories.

## Visualisation

The 128-D molecular representations are projected using:

* PCA
* t-SNE
* UMAP

The notebook generates visualisations for:

* TPT coverage categories
* committor probability
* stationary probability
* effective current
* embedding quality diagnostics

## Molecular Clustering

A downstream clustering pipeline applies **UMAP + HDBSCAN** to the learned molecular representations.

UMAP first reduces the 128-D embeddings to a lower-dimensional clustering space. HDBSCAN then discovers dense molecular groups without requiring predefined labels.

Cluster analysis includes:

* fuzzy cluster membership
* centroid-based cosine distance
* membership entropy
* intra-cluster cosine compactness
* inter-cluster centroid similarity

The experiment retains **30 major molecular clusters** for downstream analysis.

## Pipeline

```text
Molecular Reaction Graph
        │
        ▼
Thermodynamic Edge Weighting
   ΔG + Competition + Robustness
        │
        ▼
Entropy-Based Temperature Scaling
        │
        ▼
Weighted Markov Transition Graph
        │
        ▼
Louvain Graph Partitioning
        │
        ▼
20 Source-Target Pairs
        │
        ▼
Transition Path Theory
   ├── Stationary Distribution π
   ├── Forward Committor q+
   └── Effective Current f+
        │
        ▼
Multi-Pair TPT Aggregation
        │
        ▼
TPT-Biased Random Walks
        │
        ▼
1.23M Reaction Trajectories
        │
        ▼
Skip-Gram / Word2Vec
        │
        ▼
128-D Molecular Embeddings
        │
        ├── Similarity Propagation
        ├── PCA / t-SNE / UMAP
        └── UMAP + HDBSCAN Clustering
```

## Generated Artifacts

The pipeline generates several reusable outputs:

```text
reaction_graph_t_scaled.pkl
    Temperature-scaled reaction graph

reaction_graph_tpt_paper3.pkl
    Reaction graph enriched with TPT quantities and embeddings

molecule_embeddings_paper3.pkl
    {SMILES -> 128-D embedding} dictionary

embeddings_index_paper3.pkl
    SMILES-to-embedding-file lookup

embeddings_individual_paper3/
    Individual NumPy embedding files

mol_cluster_artifacts.pkl
    Molecular clustering results and fuzzy memberships
```

The enriched graph stores:

### Node Attributes

```text
q_plus
    Multi-pair averaged forward committor

pi
    Stationary probability

embedding
    128-D molecular embedding

similarity_to_source
    Propagated cosine similarity
```

### Edge Attributes

```text
weight
    Thermodynamic transition weight

restriction
    Competing-reaction penalty

robustness
    Energetic viability indicator

effective_current
    Aggregated multi-pair TPT reactive current
```

## Requirements

Main dependencies:

```text
numpy
scipy
networkx
gensim
scikit-learn
matplotlib
optuna
umap-learn
hdbscan
python-louvain
```

Optional GPU acceleration:

```text
cupy
cupyx
```

## Running the Notebook

The implementation is contained in:

```text
Reaction_graph_embedding.ipynb
```

The notebook expects the initial molecule-centric reaction graph:

```text
reaction_graph_molecule_centric.pkl
```

Run the notebook sequentially to:

```text
1. Construct thermodynamic edge weights
2. Optimise transition-weight scaling
3. Analyse graph connectivity
4. Compute multi-pair TPT
5. Generate TPT-guided random walks
6. Train molecular embeddings
7. Evaluate embedding structure
8. Visualise the latent space
9. Cluster molecular representations
```

GPU acceleration is automatically used for sparse TPT calculations when CuPy is available.

## Research Basis

This implementation builds on:

**Paula Mercurio and Di Liu.
“Clustering Molecular Energy Landscapes by Adaptive Network Embedding.”
Journal of Materials Informatics, 2024.
DOI: 10.20517/jmi.2023.40**

The original work combines network embedding with energy-landscape sampling and Transition Path Theory.

This implementation adapts the framework to a **123K-node molecular reaction network** and adds several extensions for large, disconnected reaction graphs:

* thermodynamic reaction weighting using `ΔG` and kinetic competition
* entropy-optimised temperature scaling
* automated graph-aware source-target selection
* **multi-pair TPT aggregation**
* GPU-accelerated sparse committor computation
* TPT-guided Skip-Gram molecular embeddings
* embedding-space evaluation and molecular clustering

## Tech Stack

**Python · NetworkX · SciPy · CuPy · NumPy · Gensim · Word2Vec · Optuna · Louvain · UMAP · HDBSCAN · Scikit-learn**
