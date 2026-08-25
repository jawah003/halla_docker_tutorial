# Where this data came from

The four data files in this directory are a small subset of one public dataset, reshaped into the
format HAllA expects. They are committed so you can run the tutorial without downloading anything.

## Source

**Dataset:** `KIM_ADENOMAS_2020` from the Curated Gut Microbiome-Metabolome Data Resource.

- Resource: <https://github.com/borenstein-lab/microbiome-metabolome-curated-data> (MIT licence)
- Resource paper: Muller E, Algavi YM, Borenstein E. *The gut microbiome-metabolome dataset
  collection: a curated resource for integrative meta-analysis.* npj Biofilms and Microbiomes 8, 79
  (2022). <https://doi.org/10.1038/s41522-022-00345-5>
- Original study: Kim M, Vogtmann E, Ahlquist DA, et al. *Fecal Metabolomic Signatures in Colorectal
  Adenoma Patients Are Associated with Gut Microbiota and Early Events of Colorectal Cancer
  Pathogenesis.* mBio 11(1) (2020). <https://doi.org/10.1128/mBio.03186-19>

If you publish anything derived from this data, cite the original study and the resource paper.

## What the study measured

240 participants (102 controls, 102 with colorectal adenomas, 36 with carcinoma), each contributing
one stool sample profiled two ways: 16S rRNA sequencing for the microbiome and Metabolon's
untargeted LC-MS platform for the metabolome.

## What `prepare_data.py` did to it

Reshaping (both matrices):

- Transposed from the source layout (samples as rows) to HAllA's layout (features as rows).
- Replaced `;` with `,` in every feature name. HAllA joins the members of a cluster with `;` in
  `sig_clusters.txt`, so a feature whose own name contains a semicolon gets silently split into
  fragments that match nothing downstream. Both axes were affected: GTDB lineages are
  semicolon-delimited, and metabolite names like `docosapentaenoate (n3 DPA; 22:5n3)` carry one too.
  The original strings are preserved in the `label` column of the metadata files.

Taxa (`Y_taxa.txt`, 30 features):

- Started from `genera.tsv`, the genus-level relative abundance table (499 genera).
- Kept genera detected in at least 50% of samples (54 remained).
- Dropped lineages with no genus name assigned — 16S reads that could only be placed to family
  (49 remained).
- Kept the 30 most abundant of those, by mean relative abundance.
- Transformed to `log10(relative abundance + 1e-5)`. Absent genera therefore sit at `-5`.

Metabolites (`X_metabolites.txt`, 40 features):

- Started from `mtb.tsv` (462 metabolites), keeping only those flagged
  `High.Confidence.Annotation` with a named Metabolon superpathway.
- Kept the 40 most variable, capped at 8 per superpathway so no single class dominates the figure.
- No transformation: the resource publishes these values already log-scaled.

## The metadata files

Neither is used by HAllA. Both are read by the notebook when it draws the annotated hallagram, and
both are indexed by feature ID so they join to the rows of the matrices above.

- `taxa_metadata.txt` — taxonomic ranks parsed from the GTDB lineage, the full lineage, a readable
  display label, and mean abundance.
- `metabolite_metadata.txt` — Metabolon superpathway and subpathway, HMDB and KEGG identifiers, a
  readable display label, and `shannon_rho` / `shannon_fdr`: the Spearman correlation of each
  metabolite with per-sample Shannon diversity, with a Benjamini-Hochberg FDR. Diversity was
  computed from the *full* 499-genus table, not the 30 selected taxa. This column exists to
  demonstrate a numeric annotation strip; it is not a finding.
- `taxa_prevalence.txt` — percentage of the 240 samples in which each genus was detected. Drives the
  marginal barplot.
- `sample_metadata.txt` — study group, age, sex and race per sample. Not used by the tutorial;
  included so you can subset the samples yourself.

## Regenerating

```bash
git clone https://github.com/borenstein-lab/microbiome-metabolome-curated-data
python scripts/prepare_data.py --source microbiome-metabolome-curated-data
```
