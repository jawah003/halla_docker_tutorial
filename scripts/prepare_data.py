#!/usr/bin/env python3
"""Build the tutorial's HAllA inputs from the curated gut microbiome-metabolome collection.

The committed files under `data/` were produced by this script, so you do not need
to run it to follow the tutorial. Run it if you want to change the filters, use a
different dataset from the collection, or see exactly how the inputs were derived.

Source (clone it yourself if you want to re-run this):
    https://github.com/borenstein-lab/microbiome-metabolome-curated-data

Usage:
    python scripts/prepare_data.py --source /path/to/microbiome-metabolome-curated-data
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

DATASET = "KIM_ADENOMAS_2020"

# --- Filters -----------------------------------------------------------------
MIN_PREVALENCE = 0.50   # keep genera detected in at least this fraction of samples
N_TAXA = 30             # then keep this many, by mean relative abundance
N_METABOLITES = 40      # keep this many metabolites, by variance
MAX_PER_SUPERPATHWAY = 8  # ...but no more than this many from any one superpathway
PSEUDOCOUNT = 1e-5      # added before log10 of relative abundance
SHANNON_FDR_ALPHA = 0.05


def sanitise(name):
    """Make a feature name safe to use as a HAllA feature ID.

    HAllA joins the feature names of a cluster with ';' in sig_clusters.txt, so a
    feature whose own name contains ';' gets silently split into two names that
    match nothing. Both axes of this dataset are affected: GTDB lineages are
    ';'-delimited, and metabolite names like "docosapentaenoate (n3 DPA; 22:5n3)"
    carry a semicolon of their own. Commas are safe, so ';' becomes ','.
    """
    name = str(name).replace(";", ",")
    name = re.sub(r"\s+", " ", name).strip()
    return name


def parse_lineage(lineage):
    """Split a GTDB-style 'd__..;p__..;..;g__..' string into a rank -> name dict."""
    ranks = {"d": "Domain", "p": "Phylum", "c": "Class", "o": "Order",
             "f": "Family", "g": "Genus"}
    out = {v: "" for v in ranks.values()}
    for part in str(lineage).split(";"):
        part = part.strip()
        if "__" not in part:
            continue
        prefix, _, value = part.partition("__")
        if prefix in ranks:
            out[ranks[prefix]] = value.strip()
    return out


def taxon_id(lineage):
    """Short, readable, unique-enough ID for a genus.

    56 of the 499 genera in this dataset have an empty genus field (the 16S read
    could be placed no deeper than family), so those fall back to the deepest
    named rank. This is the kind of mess you get in real taxonomy tables.
    """
    ranks = parse_lineage(lineage)
    if ranks["Genus"]:
        return sanitise(ranks["Genus"])
    for rank in ("Family", "Order", "Class", "Phylum"):
        if ranks[rank]:
            return sanitise(f"{ranks[rank]}_unclassified")
    return sanitise(lineage)


def shannon(row):
    p = np.asarray(row, dtype=float)
    p = p[p > 0]
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def pick_metabolites(mtb, mtb_map):
    """High-confidence, annotated metabolites: most variable, capped per superpathway."""
    conf = mtb_map["High.Confidence.Annotation"].astype(str).str.upper().isin({"TRUE", "T"})
    annotated = mtb_map["Superpathway"].notna() & (mtb_map["Superpathway"].astype(str) != "")
    eligible = mtb_map.loc[conf & annotated, ["Compound", "Superpathway", "Subpathway", "HMDB", "KEGG"]]
    eligible = eligible[eligible["Compound"].isin(mtb.columns)]

    variance = mtb[eligible["Compound"]].var().sort_values(ascending=False)
    eligible = eligible.set_index("Compound").loc[variance.index]

    kept, per_path = [], {}
    for compound, row in eligible.iterrows():
        path = row["Superpathway"]
        if per_path.get(path, 0) >= MAX_PER_SUPERPATHWAY:
            continue
        kept.append(compound)
        per_path[path] = per_path.get(path, 0) + 1
        if len(kept) == N_METABOLITES:
            break
    return eligible.loc[kept]


def write_halla_matrix(df, path):
    """Write a HAllA input file: features as rows, samples as columns, first header cell '#'."""
    out = df.copy()
    out.index.name = "#"
    out.to_csv(path, sep="\t", float_format="%.6g")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True,
                    help="path to a clone of borenstein-lab/microbiome-metabolome-curated-data")
    ap.add_argument("--out", default=None, help="output directory (default: ../data)")
    args = ap.parse_args()

    src = Path(args.source).expanduser() / "data" / "processed_data" / DATASET
    if not src.is_dir():
        sys.exit(f"Dataset directory not found: {src}")
    out_dir = Path(args.out) if args.out else Path(__file__).resolve().parent.parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    genera = pd.read_csv(src / "genera.tsv", sep="\t", index_col=0)
    mtb = pd.read_csv(src / "mtb.tsv", sep="\t", index_col=0)
    mtb_map = pd.read_csv(src / "mtb.map.tsv", sep="\t")
    sample_meta = pd.read_csv(src / "metadata.tsv", sep="\t")

    samples = [s for s in genera.index if s in set(mtb.index)]
    genera, mtb = genera.loc[samples], mtb.loc[samples]
    print(f"{DATASET}: {len(samples)} samples with both microbiome and metabolome data")

    # --- Y: taxa -------------------------------------------------------------
    prevalence = (genera > 0).mean()
    prevalent = prevalence[prevalence >= MIN_PREVALENCE].index
    named = [c for c in prevalent if parse_lineage(c)["Genus"]]
    ranked = genera[named].mean().sort_values(ascending=False)
    taxa_cols = list(ranked.index[:N_TAXA])
    print(f"taxa: {len(prevalence)} genera -> {len(prevalent)} at >={MIN_PREVALENCE:.0%} prevalence "
          f"-> {len(named)} with a named genus -> top {len(taxa_cols)} by mean abundance")

    ids = [taxon_id(c) for c in taxa_cols]
    if len(set(ids)) != len(ids):
        sys.exit("Duplicate taxon IDs after sanitising; disambiguate before continuing.")

    Y = np.log10(genera[taxa_cols] + PSEUDOCOUNT).T
    Y.index = ids
    write_halla_matrix(Y, out_dir / "Y_taxa.txt")

    taxa_meta = pd.DataFrame([parse_lineage(c) for c in taxa_cols], index=ids)
    taxa_meta["lineage"] = [sanitise(c) for c in taxa_cols]
    taxa_meta["label"] = [f"{r['Genus']} ({r['Family']})" if r["Family"] else r["Genus"]
                          for r in (parse_lineage(c) for c in taxa_cols)]
    taxa_meta["mean_rel_abundance"] = genera[taxa_cols].mean().values
    taxa_meta["mean_log10_abundance"] = Y.mean(axis=1).values
    taxa_meta.index.name = "feature"
    taxa_meta.to_csv(out_dir / "taxa_metadata.txt", sep="\t", float_format="%.6g")

    n_present = (genera[taxa_cols] > 0).sum()
    pd.DataFrame({
        "feature": ids,
        "prevalence_pct": (n_present / len(samples) * 100).round(2).values,
        "n_present": n_present.values,
        "n_samples": len(samples),
    }).to_csv(out_dir / "taxa_prevalence.txt", sep="\t", index=False)

    # --- X: metabolites ------------------------------------------------------
    chosen = pick_metabolites(mtb, mtb_map)
    mtb_ids = [sanitise(c) for c in chosen.index]
    if len(set(mtb_ids)) != len(mtb_ids):
        sys.exit("Duplicate metabolite IDs after sanitising; disambiguate before continuing.")
    print(f"metabolites: {mtb.shape[1]} measured -> top {len(mtb_ids)} by variance, "
          f"<={MAX_PER_SUPERPATHWAY} per superpathway")

    # Already log-scaled and median-centred in the source collection, so no transform.
    X = mtb[list(chosen.index)].T
    X.index = mtb_ids
    write_halla_matrix(X, out_dir / "X_metabolites.txt")

    # Numeric annotation: how each metabolite tracks whole-community diversity.
    # Computed from the FULL genus table, not the 30 selected taxa.
    diversity = genera.apply(shannon, axis=1)
    rho, pval = [], []
    for compound in chosen.index:
        r, p = stats.spearmanr(mtb[compound], diversity)
        rho.append(r)
        pval.append(p)
    fdr = multipletests(pval, method="fdr_bh")[1]

    mtb_meta = chosen.reset_index(drop=True)
    mtb_meta.index = mtb_ids
    mtb_meta["label"] = [sanitise(c) for c in chosen.index]
    mtb_meta["shannon_rho"] = np.round(rho, 4)
    mtb_meta["shannon_fdr"] = fdr
    # Both the coefficient and its FDR are written out unmasked. The notebook masks
    # rho by fdr at render time, so you can change the threshold without re-running
    # this script.
    mtb_meta.index.name = "feature"
    mtb_meta.to_csv(out_dir / "metabolite_metadata.txt", sep="\t", float_format="%.6g")
    print(f"        {int((fdr < SHANNON_FDR_ALPHA).sum())}/{len(fdr)} correlate with Shannon "
          f"diversity at FDR < {SHANNON_FDR_ALPHA}")

    # --- Sample metadata (not used by HAllA; handy for subsetting) -----------
    keep = [c for c in ("Sample", "Study.Group", "Age", "Gender", "Race") if c in sample_meta.columns]
    sample_meta = sample_meta[keep].set_index("Sample").reindex(samples)
    sample_meta.index.name = "sample"
    sample_meta.to_csv(out_dir / "sample_metadata.txt", sep="\t")

    print(f"\nwrote X {X.shape[0]}x{X.shape[1]} and Y {Y.shape[0]}x{Y.shape[1]} to {out_dir}")
    for name, count in sample_meta["Study.Group"].value_counts().items():
        print(f"  {name}: {count} samples")


if __name__ == "__main__":
    main()
