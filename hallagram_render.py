"""HAllAgram rendering engine: metadata colour strips, custom labels, marginal barplot.

A standalone matplotlib/seaborn renderer for HAllA output. It reads the files HAllA
writes into its output directory and draws the hallagram itself, so you can annotate
it in ways the built-in `--hallagram` cannot: categorical and numeric metadata strips
on either axis, human-readable tick labels, and a left-margin barplot.

It does NOT import or require the `halla` package, so it runs in your own Python
environment rather than inside the HAllA container.

Public API:
    build_hallagram_inputs(out_dir) -> dict
    render_native_hallagram(inputs, ...) -> (fig, legend_fig)

Requires: pandas, numpy, scipy, matplotlib, seaborn.
"""
import os
import re

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import matplotlib.patheffects as path_effects
import seaborn as sns

# ---------------------------------------------------------------------------
# Palette / colour helpers
# ---------------------------------------------------------------------------
# 23-colour expanded palette (Paired-like + Set3-like + black)
_CB_PALETTE = [
    '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
    '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6',
    '#6A3D9A', '#FFFF99', '#B15928', '#8DD3C7',
    '#FB8072', '#80B1D3', '#D9D9D9', '#BC80BD',
    '#CCEBC5', '#FFED6F', '#FDB462', '#FCCDE5',
    '#BEBADA', '#FB9A99', 'black',
]
_MISSING_COLOR = '#888888'
_COLORMAP_OPTIONS = ['RdBu_r', 'RdBu', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdYlBu', 'Spectral', 'Viridis', 'Plasma']


def _rainbow_palette(n):
    """Generate an R-style rainbow palette of *n* hex colours."""
    cmap = matplotlib.colormaps.get_cmap('rainbow')
    return [matplotlib.colors.to_hex(cmap(i / max(n - 1, 1))) for i in range(n)]


def _is_missing(value):
    if pd.isna(value):
        return True
    return str(value).strip().lower() in {'', 'na', 'nan', 'none', 'missing', 'null'}


def _is_numeric_metadata(metadata_df, column):
    values = metadata_df[column]
    non_missing = [v for v in values if not _is_missing(v)]
    if not non_missing:
        return False
    return pd.to_numeric(pd.Series(non_missing), errors='coerce').notna().all()


def _build_color_map(metadata_df, ordered_features, column, palette_offset=0, global_cat_color=None,
                     category_order=None, category_colors=None):
    """Return (rgba_array [n_features, 4], category_to_hex dict).

    category_colors (optional dict {category: color}) overrides the palette for those
    categories; any category not listed falls back to the global/palette assignment.
    """
    values = metadata_df.loc[[str(f) for f in ordered_features if str(f) in metadata_df.index], column]
    categories = sorted({str(v) for v in values if not _is_missing(v)})
    if category_order:
        ordered = [str(v) for v in category_order if str(v) in categories]
        categories = ordered + [cat for cat in categories if cat not in ordered]

    cat_color = {}
    if global_cat_color is not None:
        for cat in categories:
            cat_color[cat] = global_cat_color.get(cat, _MISSING_COLOR)
    else:
        n_cats = len(categories)
        n_pal  = len(_CB_PALETTE)
        if n_cats <= n_pal:
            palette = _CB_PALETTE
        else:
            palette = _rainbow_palette(n_cats)
        for i, cat in enumerate(categories):
            cat_color[cat] = palette[(palette_offset + i) % len(palette)]

    if category_colors:
        for cat in categories:
            if cat in category_colors:
                cat_color[cat] = category_colors[cat]

    # NOTE: a missing ('NA') legend entry is added later by the renderer *only* when the
    # column actually has missing values (see _col_has_missing), so it is not injected here.

    rgba = np.zeros((len(ordered_features), 4), dtype=float)
    for i, f in enumerate(ordered_features):
        key = str(f)
        if key in metadata_df.index:
            v = metadata_df.loc[key, column]
            hex_c = _MISSING_COLOR if _is_missing(v) else cat_color[str(v)]
        else:
            hex_c = _MISSING_COLOR
        rgba[i] = matplotlib.colors.to_rgba(hex_c)

    return rgba, cat_color


def _col_has_missing(metadata_df, ordered_features, column):
    """True if any plotted feature is absent from the metadata or has a missing value."""
    for f in ordered_features:
        key = str(f)
        if key not in metadata_df.index or _is_missing(metadata_df.loc[key, column]):
            return True
    return False


def _build_numeric_color_map(metadata_df, ordered_features, column, cmap_name, center_zero=False):
    """Return (rgba_array [n_features, 4], vmin, vmax).

    If center_zero, the colour scale is symmetrized about 0 so that 0 maps to the
    colormap midpoint (correct warm=positive / cold=negative semantics for signed values).
    """
    numeric_values = np.full(len(ordered_features), np.nan, dtype=float)
    for i, feature in enumerate(ordered_features):
        key = str(feature)
        if key not in metadata_df.index:
            continue
        value = metadata_df.loc[key, column]
        if _is_missing(value):
            continue
        try:
            numeric_values[i] = float(value)
        except (TypeError, ValueError):
            continue

    valid = numeric_values[~np.isnan(numeric_values)]
    if len(valid) == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin, vmax = float(np.min(valid)), float(np.max(valid))

    if center_zero and len(valid) > 0:
        m = max(abs(vmin), abs(vmax))
        if m > 0:
            vmin, vmax = -m, m

    if vmin == vmax:
        norm = matplotlib.colors.Normalize(vmin=vmin - 0.5, vmax=vmax + 0.5, clip=True)
    else:
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    if cmap_name not in matplotlib.colormaps:
        raise ValueError(f"Invalid numeric metadata colormap: {cmap_name}. Valid options: {_COLORMAP_OPTIONS}")
    cmap = matplotlib.colormaps.get_cmap(cmap_name)

    rgba = np.zeros((len(ordered_features), 4), dtype=float)
    missing_rgba = matplotlib.colors.to_rgba(_MISSING_COLOR)
    for i, value in enumerate(numeric_values):
        rgba[i] = missing_rgba if np.isnan(value) else cmap(norm(value))

    return rgba, vmin, vmax


# ---------------------------------------------------------------------------
# HAllA output parsing
# ---------------------------------------------------------------------------
def parse_performance_alpha(perf_path):
    """Extract the FDR alpha threshold used in the HAllA run."""
    if os.path.exists(perf_path):
        with open(perf_path, 'r') as f:
            for line in f:
                if 'fdr alpha' in line:
                    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
                    if match:
                        return float(match.group(0))
                    print(f"Invalid FDR alpha value in {perf_path}: {line.strip()}")
                    return 0.05
    return 0.05


def build_hallagram_inputs(out_dir):
    """Build native hallagram inputs from standard HAllA output files."""
    x_table        = os.path.join(out_dir, "X.tsv")
    y_table        = os.path.join(out_dir, "Y.tsv")
    assoc_table    = os.path.join(out_dir, "all_associations.txt")
    perf_table     = os.path.join(out_dir, "performance.txt")
    sig_table      = os.path.join(out_dir, "sig_clusters.txt")
    x_linkage_file = os.path.join(out_dir, "X_linkage.npy")
    y_linkage_file = os.path.join(out_dir, "Y_linkage.npy")

    required_files = [x_table, y_table, assoc_table, perf_table, sig_table]
    missing_required = [p for p in required_files if not os.path.exists(p)]
    if missing_required:
        raise FileNotFoundError(
            f"Missing required HAllA output files: {', '.join(missing_required)}"
        )

    x_features = pd.read_table(x_table, index_col=0).index.to_numpy()
    y_features = pd.read_table(y_table, index_col=0).index.to_numpy()
    x_map = {name: i for i, name in enumerate(x_features.tolist())}
    y_map = {name: i for i, name in enumerate(y_features.tolist())}

    associations = pd.read_table(assoc_table)
    alpha        = parse_performance_alpha(perf_table)

    sim_table        = np.zeros((len(x_features), len(y_features)))
    fdr_reject_table = np.zeros((len(x_features), len(y_features)), dtype=bool)

    q_col = 'q-values' if 'q-values' in associations.columns else 'p-values'
    for _, row in associations.iterrows():
        xi = x_map[row['X_features']]
        yi = y_map[row['Y_features']]
        sim_table[xi, yi]        = row['association']
        fdr_reject_table[xi, yi] = row[q_col] < alpha

    sig_clusters      = pd.read_table(sig_table)
    significant_blocks = []
    for _, row in sig_clusters.iterrows():
        xb = [x_map[n] for n in str(row['cluster_X']).split(';') if n in x_map]
        yb = [y_map[n] for n in str(row['cluster_Y']).split(';') if n in y_map]
        if xb and yb:
            significant_blocks.append([xb, yb])

    if os.path.exists(x_linkage_file) and os.path.exists(y_linkage_file):
        clust_x_idx = sch.to_tree(np.load(x_linkage_file)).pre_order()
        clust_y_idx = sch.to_tree(np.load(y_linkage_file)).pre_order()
    else:
        clust_x_idx = list(range(len(x_features)))
        clust_y_idx = list(range(len(y_features)))

    return {
        'x_features':        x_features,
        'y_features':        y_features,
        'sim_table':         sim_table,
        'fdr_reject_table':  fdr_reject_table,
        'significant_blocks': significant_blocks,
        'clust_x_idx':       clust_x_idx,
        'clust_y_idx':       clust_y_idx,
    }


# ---------------------------------------------------------------------------
# Core rendering function
# ---------------------------------------------------------------------------
def render_native_hallagram(
    inputs,
    dpi,
    cmap,
    trim,
    block_num,
    x_metadata_df=None,
    y_metadata_df=None,
    x_metadata_columns=None,
    y_metadata_columns=None,
    x_tick_label_col=None,
    y_tick_label_col=None,
    x_dataset_label=None,
    y_dataset_label=None,
    x_label_size=0,
    y_label_size=0,
    metadata_font_size=0,
    numeric_metadata_cmap='RdBu_r',
    category_orders=None,
    category_colors=None,
    numeric_center_zero=None,
    na_legend_labels=None,
    numeric_cmaps=None,
    hide_na_legend=None,
    cbar_label='Association',
    y_barplot=None,
    legend_panel_width_in=3.0,
    legend_gap_in=0.35,
    barplot_width_in=1.3,
    legend_scale=1.0,
    legend_titles=None,
):
    if not inputs['significant_blocks']:
        raise ValueError("No significant clusters found; hallagram cannot be generated.")

    # 1. Data preparation
    eff_bn     = len(inputs['significant_blocks']) if block_num < 0 else block_num
    top_blocks = inputs['significant_blocks'][:eff_bn]

    if trim:
        inc_x = sorted({idx for blk in top_blocks for idx in blk[0]})
        inc_y = sorted({idx for blk in top_blocks for idx in blk[1]})
    else:
        inc_x = list(range(len(inputs['x_features'])))
        inc_y = list(range(len(inputs['y_features'])))

    cx_idx = np.asarray([i for i in inputs['clust_x_idx'] if i in inc_x])
    cy_idx = np.asarray([i for i in inputs['clust_y_idx'] if i in inc_y])

    # Transpose so rows=Y features (vertical) and cols=X features (horizontal)
    clust_sim  = inputs['sim_table'][cx_idx, :][:, cy_idx].T
    clust_fdr  = inputs['fdr_reject_table'][cx_idx, :][:, cy_idx].T
    clust_x_ft = np.asarray(inputs['x_features'])[cx_idx]
    clust_y_ft = np.asarray(inputs['y_features'])[cy_idx]

    # After transpose: n_rows = Y features, n_cols = X features
    n_rows, n_cols = clust_sim.shape

    x_ori2clust = {orig: i for i, orig in enumerate(cx_idx)}
    y_ori2clust = {orig: i for i, orig in enumerate(cy_idx)}

    # Proportional Font Sizing
    auto_y = max(4, min(12, 120 // max(1, n_rows)))
    auto_x = max(4, min(12, 120 // max(1, n_cols)))

    y_font_size = y_label_size if y_label_size > 0 else auto_y
    x_font_size = x_label_size if x_label_size > 0 else auto_x
    metadata_label_size = metadata_font_size if metadata_font_size > 0 else 7

    # 2. Validate metadata
    x_metadata_columns = x_metadata_columns or []
    y_metadata_columns = y_metadata_columns or []
    category_orders = category_orders or {}
    category_colors = category_colors or {}
    numeric_center_zero = set(numeric_center_zero or [])
    na_legend_labels = na_legend_labels or {}
    numeric_cmaps = numeric_cmaps or {}
    hide_na_legend = set(hide_na_legend or [])
    legend_titles = legend_titles or {}

    if x_metadata_columns and x_metadata_df is None:
        raise ValueError("X metadata columns selected but no X metadata file provided.")
    if y_metadata_columns and y_metadata_df is None:
        raise ValueError("Y metadata columns selected but no Y metadata file provided.")

    if x_metadata_df is not None and x_metadata_columns:
        missing = [f for f in clust_x_ft if str(f) not in x_metadata_df.index]
        if missing:
            raise ValueError(f"X metadata missing features: {', '.join(str(m) for m in missing[:10])}")
    if y_metadata_df is not None and y_metadata_columns:
        missing = [f for f in clust_y_ft if str(f) not in y_metadata_df.index]
        if missing:
            raise ValueError(f"Y metadata missing features: {', '.join(str(m) for m in missing[:10])}")

    x_numeric_cols = {col for col in x_metadata_columns if _is_numeric_metadata(x_metadata_df, col)}
    y_numeric_cols = {col for col in y_metadata_columns if _is_numeric_metadata(y_metadata_df, col)}

    # 3. Build Global Colour Map
    all_unique_classes = set()
    if x_metadata_df is not None and x_metadata_columns:
        for col in x_metadata_columns:
            if col in x_numeric_cols:
                continue
            vals = x_metadata_df.loc[[str(f) for f in clust_x_ft if str(f) in x_metadata_df.index], col]
            all_unique_classes.update({str(v) for v in vals if not _is_missing(v)})
    if y_metadata_df is not None and y_metadata_columns:
        for col in y_metadata_columns:
            if col in y_numeric_cols:
                continue
            vals = y_metadata_df.loc[[str(f) for f in clust_y_ft if str(f) in y_metadata_df.index], col]
            all_unique_classes.update({str(v) for v in vals if not _is_missing(v)})

    all_unique_classes = sorted(list(all_unique_classes))
    global_cat_color = None

    # If total unique classes fits in our palette, assign universally
    if len(all_unique_classes) <= len(_CB_PALETTE):
        global_cat_color = {}
        for i, cat in enumerate(all_unique_classes):
            global_cat_color[cat] = _CB_PALETTE[i]

    # 4. Build colour arrays
    x_strips, x_legend = [], {}
    for ci, col in enumerate(x_metadata_columns):
        has_missing = _col_has_missing(x_metadata_df, clust_x_ft, col)
        na_label = na_legend_labels.get(col, 'NA')
        if col in x_numeric_cols:
            cmn = numeric_cmaps.get(col, numeric_metadata_cmap)
            rgba, vmin, vmax = _build_numeric_color_map(
                x_metadata_df, clust_x_ft, col, cmn,
                center_zero=(col in numeric_center_zero),
            )
            x_legend[col] = {'type': 'numeric', 'vmin': vmin, 'vmax': vmax, 'cmap': cmn,
                             'col': col, 'has_missing': has_missing, 'na_label': na_label}
        else:
            rgba, cc = _build_color_map(
                x_metadata_df, clust_x_ft, col, palette_offset=ci * 2, global_cat_color=global_cat_color,
                category_order=category_orders.get(col), category_colors=category_colors.get(col),
            )
            x_legend[col] = {'type': 'categorical', 'colors': cc,
                             'col': col, 'has_missing': has_missing, 'na_label': na_label}
        x_strips.append((col, rgba))

    y_strips, y_legend = [], {}
    for ci, col in enumerate(y_metadata_columns):
        has_missing = _col_has_missing(y_metadata_df, clust_y_ft, col)
        na_label = na_legend_labels.get(col, 'NA')
        if col in y_numeric_cols:
            cmn = numeric_cmaps.get(col, numeric_metadata_cmap)
            rgba, vmin, vmax = _build_numeric_color_map(
                y_metadata_df, clust_y_ft, col, cmn,
                center_zero=(col in numeric_center_zero),
            )
            y_legend[col] = {'type': 'numeric', 'vmin': vmin, 'vmax': vmax, 'cmap': cmn,
                             'col': col, 'has_missing': has_missing, 'na_label': na_label}
        else:
            rgba, cc = _build_color_map(
                y_metadata_df, clust_y_ft, col, palette_offset=len(x_metadata_columns) * 2 + ci * 2, global_cat_color=global_cat_color,
                category_order=category_orders.get(col), category_colors=category_colors.get(col),
            )
            y_legend[col] = {'type': 'categorical', 'colors': cc,
                             'col': col, 'has_missing': has_missing, 'na_label': na_label}
        y_strips.append((col, rgba))

    n_xs = len(x_strips)
    n_ys = len(y_strips)

    # 5. Custom tick labels
    def _get_tick_labels(meta_df, features, col):
        if meta_df is None or col is None or col not in meta_df.columns:
            return [str(f) for f in features]
        labels = []
        for f in features:
            key = str(f)
            if key in meta_df.index:
                labels.append(str(meta_df.loc[key, col]))
            else:
                labels.append(key)
        return labels

    # After transpose: rows = Y features, cols = X features
    row_labels = _get_tick_labels(y_metadata_df, clust_y_ft, y_tick_label_col)
    col_labels = _get_tick_labels(x_metadata_df, clust_x_ft, x_tick_label_col)

    # 6. Figure (heatmap only; the BrBG scale + metadata legends are drawn in the
    #    unified legend column on the right; there is no attached left colorbar)
    cell_size = 0.3
    fig_w = max(5.0, cell_size * n_cols + 4.0)
    fig_h = max(5.0, cell_size * n_rows + 4.0)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    ax_heat = fig.add_subplot(1, 1, 1)

    # 7. Draw heatmap (no seaborn colorbar; scale shown in the unified legend)
    vmax = np.max(clust_sim)
    vmin = np.min(clust_sim)
    if vmin < 0 and vmax > 0:
        vmax = max(abs(vmin), vmax)
        vmin = -vmax

    sns.set_style('whitegrid')
    sns.heatmap(
        clust_sim,
        xticklabels=col_labels,
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        square=False,
        zorder=3,
        cbar=False,
        ax=ax_heat,
    )
    ax_heat.yaxis.tick_right()
    ax_heat.yaxis.set_label_position('right')
    ax_heat.tick_params(axis='y', length=0)
    ax_heat.set_yticklabels(row_labels, rotation=0, ha='left', fontsize=y_font_size)
    ax_heat.tick_params(axis='x', rotation=90, labelsize=x_font_size)
    for tick in ax_heat.get_xticklabels():
        tick.set_ha('center')
    ax_heat.set_xlim(0, n_cols)
    ax_heat.set_ylim(n_rows, 0)
    ax_heat.margins(0)

    if n_ys:
        ax_heat.set_yticklabels([])
        ax_heat.tick_params(axis='y', length=0, pad=0)

    # Significance dots
    for i in range(n_cols):
        for j in range(n_rows):
            if clust_fdr[j, i]:
                ax_heat.scatter(i + .5, j + .5, c='black', marker='o', zorder=3, s=25)
                ax_heat.scatter(i + .5, j + .5, c='white', marker='o', zorder=3, s=10)

    # Block borders + rank labels
    # After transpose: cx = X feature positions (columns), cy = Y feature positions (rows)
    bw = 1.65
    ts = 10
    for rank, block in enumerate(top_blocks):
        cx = [x_ori2clust[i] for i in block[0] if i in x_ori2clust]
        cy = [y_ori2clust[i] for i in block[1] if i in y_ori2clust]
        if not cx or not cy:
            continue
        ax_heat.vlines(
            [min(cx), max(cx) + 1], min(cy), max(cy) + 1,
            color='black', linewidths=bw, zorder=4, capstyle='projecting',
        )
        ax_heat.hlines(
            [min(cy), max(cy) + 1], min(cx), max(cx) + 1,
            color='black', linewidths=bw, zorder=4, capstyle='projecting',
        )
        text_size = (min(max(cx) - min(cx), max(cy) - min(cy)) + 1) * ts
        txt = ax_heat.text(
            np.mean(cx) + 0.5, np.mean(cy) + 0.5, str(rank + 1),
            size=text_size, color='white', ha='center', va='center', weight='bold',
        )
        txt.set_path_effects([
            path_effects.Stroke(linewidth=3, foreground='black'),
            path_effects.Normal(),
        ])

    fig.tight_layout(pad=0.1)
    pos = ax_heat.get_position()

    strip_w_f = 0.20 / fig_w
    strip_h_f = 0.20 / fig_h
    gap_w_f   = 0.06 / fig_w
    gap_h_f   = 0.06 / fig_h

    # 8. Y metadata strips (right side, aligned with rows = Y features)
    last_right_ax = None
    for si, (col_name, rgba) in enumerate(y_strips):
        x0   = pos.x1 + gap_w_f + si * (strip_w_f + gap_w_f)
        ax_s = fig.add_axes([x0, pos.y0, strip_w_f, pos.height])
        img  = rgba.reshape(n_rows, 1, 4)
        ax_s.imshow(img, aspect='auto', interpolation='nearest', extent=[0, 1, n_rows, 0], origin='upper')
        ax_s.set_xlim(0, 1)
        ax_s.set_ylim(n_rows, 0)
        ax_s.margins(0)
        ax_s.text(0.5, 1.01, col_name, transform=ax_s.transAxes,
                  fontsize=metadata_label_size, va='bottom', ha='center', rotation=90)
        ax_s.set_xticks([])
        ax_s.xaxis.set_ticks_position('none')

        ax_s.yaxis.set_visible(True)
        ax_s.yaxis.tick_left()
        ax_s.set_yticks(np.arange(0.5, n_rows, 1.0))
        ax_s.set_yticklabels([])
        ax_s.tick_params(axis='y', which='both', length=3, width=0.5, direction='out', left=True, right=False, pad=0)
        for sp in ax_s.spines.values():
            sp.set_visible(False)
        last_right_ax = ax_s

    if last_right_ax is not None:
        last_right_ax.yaxis.tick_right()
        last_right_ax.yaxis.set_label_position('right')
        last_right_ax.set_yticks(np.arange(0.5, n_rows, 1.0))
        last_right_ax.set_yticklabels(row_labels, fontsize=y_font_size, rotation=0, ha='left')
        last_right_ax.tick_params(axis='y', length=0, pad=3)

    # 9. X metadata strips (below, aligned with columns = X features)
    last_below_ax = None
    for si, (col_name, rgba) in enumerate(x_strips):
        y0   = pos.y0 - gap_h_f - (si + 1) * strip_h_f - si * gap_h_f
        ax_s = fig.add_axes([pos.x0, y0, pos.width, strip_h_f])
        img  = rgba.reshape(1, n_cols, 4)
        ax_s.imshow(img, aspect='auto', interpolation='nearest', extent=[0, n_cols, 1, 0], origin='upper')
        ax_s.set_xlim(0, n_cols)
        ax_s.set_ylim(0, 1)
        ax_s.margins(0)
        ax_s.text(1.005, 0.5, col_name, transform=ax_s.transAxes, fontsize=metadata_label_size, va='center', ha='left', rotation=0)

        ax_s.xaxis.set_visible(True)
        ax_s.xaxis.tick_top()
        ax_s.set_xticks(np.arange(0.5, n_cols, 1.0))
        ax_s.set_xticklabels([])
        ax_s.tick_params(axis='x', which='both', length=3, width=0.5, direction='out', top=True, bottom=False, pad=0)

        ax_s.set_yticks([])
        ax_s.yaxis.set_visible(False)
        for sp in ax_s.spines.values():
            sp.set_visible(False)
        last_below_ax = ax_s

    if last_below_ax is not None:
        ax_heat.tick_params(axis='x', labelbottom=False)
        last_below_ax.xaxis.tick_bottom()
        last_below_ax.set_xticks(np.arange(0.5, n_cols, 1))
        last_below_ax.set_xticklabels(col_labels, rotation=90, fontsize=x_font_size, ha='center')
        last_below_ax.tick_params(axis='x', length=0, pad=2)

    # 10. X and Y Custom Dataset Labels
    if x_dataset_label:
        target_ax = last_below_ax if last_below_ax is not None else ax_heat
        target_ax.set_xlabel(x_dataset_label, fontsize=12, weight='bold', labelpad=x_font_size + 10)

    if y_dataset_label:
        target_ax = last_right_ax if last_right_ax is not None else ax_heat
        target_ax.set_ylabel(y_dataset_label, fontsize=12, weight='bold', labelpad=y_font_size + 10, rotation=-90)
        target_ax.yaxis.set_label_position('right')

    # 11. Left-margin marginal barplot (optional; e.g. taxa prevalence %)
    if y_barplot:
        bvals   = y_barplot.get('values', {}) or {}
        b_label = y_barplot.get('label', '')
        b_color = y_barplot.get('color', '#4C72B0')
        row_v = np.array([float(bvals.get(str(f), np.nan)) for f in clust_y_ft], dtype=float)
        finite = row_v[np.isfinite(row_v)]
        b_max = y_barplot.get('vmax') or (float(np.nanmax(finite)) if finite.size else 1.0)
        b_w_f = y_barplot.get('width_in', barplot_width_in) / fig_w
        b_gap_f = 0.12 / fig_w
        ax_bar = fig.add_axes([pos.x0 - b_gap_f - b_w_f, pos.y0, b_w_f, pos.height])
        ycen = np.arange(0.5, n_rows, 1.0)
        ax_bar.barh(ycen, np.nan_to_num(row_v, nan=0.0), height=0.9,
                    color=b_color, edgecolor='none', zorder=3)
        ax_bar.set_ylim(n_rows, 0)
        ax_bar.set_xlim(b_max, 0)               # 0 adjacent to heatmap; bars grow left
        ax_bar.set_yticks([])
        ax_bar.set_xticks([0, b_max])
        ax_bar.tick_params(axis='x', labelsize=max(6, metadata_label_size * 0.55), length=3, pad=1)
        for sp in ('top', 'right', 'left'):
            ax_bar.spines[sp].set_visible(False)
        if b_label:
            ax_bar.set_xlabel(b_label, fontsize=metadata_label_size)

    # 12. Unified legend column on the right: BrBG colorbar (top), then X, then Y.
    #     legend_scale enlarges legend fonts + colorbars together; legend_titles overrides titles.
    fig.canvas.draw()
    _rend = fig.canvas.get_renderer()
    _lab_ax = last_right_ax if last_right_ax is not None else ax_heat
    _labs = [t for t in _lab_ax.get_yticklabels() if t.get_text()]
    right_px = (max(t.get_window_extent(_rend).x1 for t in _labs)
                if _labs else _lab_ax.get_window_extent(_rend).x1)
    x_leg = right_px / fig.bbox.width + legend_gap_in / fig_w
    leg_w = legend_panel_width_in * legend_scale / fig_w

    # (display title, legend-dict) per group; titles honor the legend_titles overrides
    legend_groups = [(cbar_label or 'Association',
                      {'type': 'cbar', 'vmin': vmin, 'vmax': vmax, 'cmap': cmap})]
    for col in x_metadata_columns:
        legend_groups.append((legend_titles.get(col, f"X: {col}"), x_legend[col]))
    for col in y_metadata_columns:
        legend_groups.append((legend_titles.get(col, f"Y: {col}"), y_legend[col]))

    fs_entry = max(6.0, metadata_label_size * 0.55) * legend_scale
    fs_title = max(7.0, metadata_label_size * 0.62) * legend_scale
    line_in  = fs_title * 1.6 / 72.0             # inches per title line (for multi-line titles)
    y_in = pos.y1 * fig_h                        # start at top of heatmap (inches)
    gap_in = 0.22 * legend_scale
    for title, cc in legend_groups:
        n_lines = str(title).count('\n') + 1
        show_na = cc.get('has_missing') and cc.get('col') not in hide_na_legend
        if cc['type'] in ('cbar', 'numeric'):
            bar_h = 0.16 * legend_scale
            pad_below = 0.24 * legend_scale       # room for the colorbar's tick labels
            h_in = n_lines * line_in + bar_h + pad_below + 0.10 * legend_scale
            y_in -= h_in
            ax_l = fig.add_axes([x_leg, (y_in + pad_below) / fig_h, leg_w, bar_h / fig_h])
            cmn = cc.get('cmap', numeric_metadata_cmap)
            sm = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(vmin=cc['vmin'], vmax=cc['vmax']),
                cmap=matplotlib.colormaps.get_cmap(cmn))
            sm.set_array([])
            cb = fig.colorbar(sm, cax=ax_l, orientation='horizontal')
            cb.set_ticks([cc['vmin'], cc['vmax']])
            cb.set_ticklabels([f"{cc['vmin']:.2g}", f"{cc['vmax']:.2g}"])
            cb.ax.tick_params(labelsize=fs_entry, length=2, pad=1)
            ax_l.set_title(title, fontsize=fs_title, fontweight='bold', pad=3, loc='left')
            if cc['type'] == 'numeric' and show_na:
                ax_l.text(1.03, 0.5, cc.get('na_label', 'NA'), transform=ax_l.transAxes,
                          va='center', ha='left', fontsize=fs_entry,
                          bbox=dict(facecolor=_MISSING_COLOR, edgecolor='#cccccc',
                                    linewidth=0.5, boxstyle='square,pad=0.2'))
        else:
            handles = [Patch(facecolor=hex_c, edgecolor='#cccccc', linewidth=0.5, label=str(catv))
                       for catv, hex_c in cc['colors'].items()]
            if show_na:
                handles.append(Patch(facecolor=_MISSING_COLOR, edgecolor='#cccccc',
                                     linewidth=0.5, label=cc.get('na_label', 'NA')))
            n = len(handles)
            ncol = 2 if n > 16 else 1
            nrow = int(np.ceil(n / ncol))
            row_in = fs_entry * 1.7 / 72.0
            h_in = n_lines * line_in + 0.14 * legend_scale + row_in * max(1, nrow)
            y_in -= h_in
            ax_l = fig.add_axes([x_leg, y_in / fig_h, leg_w, h_in / fig_h])
            ax_l.axis('off')
            leg = ax_l.legend(handles=handles, loc='upper left', frameon=False,
                              fontsize=fs_entry, ncol=ncol, title=title,
                              title_fontproperties={'weight': 'bold', 'size': fs_title},
                              handlelength=1.0, handleheight=0.9, handletextpad=0.4,
                              borderpad=0.2, labelspacing=0.25, columnspacing=1.0)
            leg._legend_box.align = 'left'
        y_in -= gap_in

    return fig, None
