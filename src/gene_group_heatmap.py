from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def load_geneset(geneset_file) -> list:
    with open(geneset_file) as f:
        return [line.strip().replace('"', "") for line in f]


def main():
    cluster_number = 1

    geneset = load_geneset(
        f"data/raw/geneset/totalcell_cluster{cluster_number}_genes.txt"
    )

    metadata = pd.read_csv("data/processed/metadata.csv")
    metadata["timepoint"] = metadata["timepoint"].astype("category")

    total_anndata = sc.read_h5ad("data/processed/anndata/sepsis.h5ad")
    total_anndata.obs = pd.merge(
        total_anndata.obs, metadata, left_on="batch", right_on="GEO"
    )
    base_genes = total_anndata.var.index
    target_genes = [gene for gene in geneset if gene in base_genes]
    # target_genes = [
    #     "HLA-DRB1",
    #     "HLA-DQB1",
    #     "PMAIP1",
    #     "CD7",
    #     "GNAS",
    #     "FKBP8",
    #     "TGFB1",
    #     "CTSW",
    #     "CD83",
    # ]
    # target_genes = [
    #     "RUFY1",
    #     "SELP",
    #     "ARPC1B",
    #     "LAMTOR1",
    #     "MPP1",
    #     "CMTM5",
    #     "SNAP23",
    #     "ITGA2B",
    #     "MYH9",
    # ]
    Path("figures/geneset").mkdir(parents=True, exist_ok=True)

    with plt.rc_context():
        fig = sc.pl.matrixplot(
            total_anndata,
            target_genes,
            groupby="timepoint",
            dendrogram=True,
            show=False,
            return_fig=True,
            swap_axes=True,
        )
        fig.savefig(
            f"figures/geneset/totalcell_cluster{cluster_number}_genes_selected_timepoint.png"
        )


if __name__ == "__main__":
    main()
