"""
Generate HTML rescue report for 520 consistently removed genes
across 5 CRISPR-KO libraries aligned to T2T-CHM13.

Run from the scripts/ directory:
    python generate_rescue_report.py

Output: results/rescue_report_520_genes.html
"""

import os
import base64
import warnings
from io import BytesIO
from datetime import date

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
BASE         = os.path.join(SCRIPT_DIR, "..")
BASE_REPORT  = os.path.join(BASE, "data", "guiderefine_output", "T2T-CHM13")
BASE_CSV     = os.path.join(BASE, "data", "removed_genes_survey")
OUT_PATH     = os.path.join(BASE, "results", "rescue_report_520_genes.html")

REPORTS = {
    "Avana":    "avana_library_full_report.xlsx",
    "Brunello": "broadgpp-brunello-library-contents_full_report.xlsx",
    "TKOv3":    "tkov3_guide_sequence_full_report.xlsx",
    "Yusa":     "yusa_hcrispr_ko_grnas_full_report.xlsx",
    "Jacquere": "Jacquere_PerGuideAnnotations_Quota4_full_report.xlsx",
}

LIB_COLORS = {
    "Avana":    "#4C72B0",
    "Brunello": "#DD8452",
    "TKOv3":    "#55A868",
    "Yusa":     "#C44E52",
    "Jacquere": "#8172B2",
}

TIER_COLORS = {
    "Tier 1": "#27ae60",
    "Tier 2": "#f39c12",
    "Tier 3": "#e74c3c",
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading removed gene list...")
rg = pd.read_csv(os.path.join(BASE_CSV, "removed_genes_all_library.csv"))
removed_genes = sorted(set.intersection(
    set(rg["avana_removed_genes_all"].dropna()),
    set(rg["brunello_removed_genes_all"].dropna()),
    set(rg["toronto_removed_genes_all"].dropna()),
    set(rg["yusa_removed_genes_all"].dropna()),
    set(rg["jacquere_removed_genes_all"].dropna()),
))
print(f"  {len(removed_genes)} genes in intersection")


def load_col(lib_name, col):
    fpath = os.path.join(BASE_REPORT, REPORTS[lib_name])
    df = pd.read_excel(fpath, sheet_name="Genes Report")
    return df.groupby("gene")[col].sum().reindex(removed_genes, fill_value=0)


print("Loading full reports (this may take ~1 min)...")
counted  = pd.DataFrame({k: load_col(k, "counted total sgRNA")    for k in REPORTS})
not_exon = pd.DataFrame({k: load_col(k, "sgRNA not target exon")  for k in REPORTS})
actual   = pd.DataFrame({k: load_col(k, "actual total sgRNA")     for k in REPORTS})
print("  Done.")

# ---------------------------------------------------------------------------
# Derived metrics
# ---------------------------------------------------------------------------
zero_all  = (counted == 0).all(axis=1)
has_2plus = (counted >= 2).any(axis=1)
has_1only = (~zero_all) & (~has_2plus)

tier = pd.Series("Tier 3", index=counted.index)
tier[has_1only]  = "Tier 2"
tier[has_2plus]  = "Tier 1"

best_lib   = counted.idxmax(axis=1)
best_lib[zero_all] = "None"
max_guides = counted.max(axis=1)

intron_total_per_gene = not_exon.sum(axis=1)
n_intron_affected     = (not_exon > 0).any(axis=1).sum()
zero_due_to_intron    = zero_all & ((not_exon > 0).any(axis=1))

# Missing gene symbols per library
missing_genes = {}
for lib_name in REPORTS:
    fpath = os.path.join(BASE_REPORT, REPORTS[lib_name])
    df_lib = pd.read_excel(fpath, sheet_name="Genes Report")
    missing_genes[lib_name] = [g for g in removed_genes if g not in set(df_lib["gene"])]

# Tier 1 rescue table (sorted by total guides)
tier1_tbl = counted[has_2plus].copy()
tier1_tbl["Max (single lib)"]  = tier1_tbl.max(axis=1)
tier1_tbl["Best Library"]      = tier1_tbl.idxmax(axis=1)
tier1_tbl["Total (all libs)"]  = tier1_tbl[list(REPORTS)].sum(axis=1)
tier1_tbl = tier1_tbl.sort_values("Total (all libs)", ascending=False)

tier2_genes = sorted(tier[tier == "Tier 2"].index.tolist())
tier3_genes = sorted(tier[tier == "Tier 3"].index.tolist())

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def gene_tags(gene_list, badge_class="gene-tag"):
    return " ".join(f'<span class="{badge_class}">{g}</span>' for g in gene_list)


def df_to_html(df, table_id=""):
    return df.to_html(
        table_id=table_id,
        classes="table table-striped table-bordered table-sm table-hover",
        border=0,
        index=True,
    )


# ---------------------------------------------------------------------------
# Figure 1 – Guide count distribution per library (grouped bar)
# ---------------------------------------------------------------------------
print("Generating Figure 1...")
fig1, ax = plt.subplots(figsize=(10, 5))
libs = list(REPORTS.keys())
x = np.arange(len(libs))
width = 0.15
buckets = {
    "0 guides":  (counted == 0).sum(),
    "1 guide":   (counted == 1).sum(),
    "2 guides":  (counted == 2).sum(),
    "3 guides":  (counted == 3).sum(),
    "≥4 guides": (counted >= 4).sum(),
}
bcolors = ["#e74c3c", "#f39c12", "#2196F3", "#27ae60", "#1a5276"]
for i, (label, vals) in enumerate(buckets.items()):
    offset = (i - 2) * width
    bars = ax.bar(x + offset, [vals[l] for l in libs], width,
                  label=label, color=bcolors[i], edgecolor="white", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(libs, fontsize=12)
ax.set_ylabel("Number of genes (out of 520)", fontsize=11)
ax.set_title("Surviving exon-targeting guide distribution for 520 removed genes", fontsize=12, pad=12)
ax.legend(title="Guide count", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=True)
ax.set_ylim(0, 310)
sns.despine()
fig1.tight_layout()
fig1_b64 = fig_to_b64(fig1)

# ---------------------------------------------------------------------------
# Figure 2 – Heatmap of guide counts (520 genes × 5 libraries)
# ---------------------------------------------------------------------------
print("Generating Figure 2 (heatmap)...")
sort_order = pd.DataFrame({
    "tier_ord": tier.map({"Tier 1": 0, "Tier 2": 1, "Tier 3": 2}),
    "neg_total": -counted.sum(axis=1),
}).sort_values(["tier_ord", "neg_total"]).index

hmap = counted.loc[sort_order]
n_t1 = (tier == "Tier 1").sum()
n_t2 = (tier == "Tier 2").sum()

cmap_hm   = ListedColormap(["#f5f5f5", "#aed6f1", "#2196F3", "#27ae60", "#1a5276"])
bounds_hm = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
norm_hm   = BoundaryNorm(bounds_hm, cmap_hm.N)

fig2, ax2 = plt.subplots(figsize=(5, 20))
im = ax2.imshow(hmap.values, aspect="auto", cmap=cmap_hm, norm=norm_hm)
ax2.set_xticks(range(len(libs)))
ax2.set_xticklabels(libs, fontsize=10, rotation=30, ha="right")
ax2.set_yticks(range(len(sort_order)))
ax2.set_yticklabels(sort_order, fontsize=3.5)
ax2.set_title("Surviving guides (520 genes × 5 libraries)\nsorted by rescue tier", fontsize=10, pad=10)
ax2.axhline(n_t1 - 0.5,        color=TIER_COLORS["Tier 2"], lw=1.5, ls="--")
ax2.axhline(n_t1 + n_t2 - 0.5, color=TIER_COLORS["Tier 3"], lw=1.5, ls="--")
cbar = fig2.colorbar(im, ax=ax2, ticks=[0, 1, 2, 3, 4], shrink=0.3)
cbar.set_label("Surviving guides", fontsize=8)
# Tier labels on right margin
for label, ypos, col in [
    ("Tier 1", n_t1 / 2,              TIER_COLORS["Tier 1"]),
    ("Tier 2", n_t1 + n_t2 / 2,       TIER_COLORS["Tier 2"]),
    ("Tier 3", n_t1 + n_t2 + 37,      TIER_COLORS["Tier 3"]),
]:
    ax2.text(5.2, ypos, label, fontsize=8, color=col, fontweight="bold",
             va="center", transform=ax2.get_yaxis_transform())
fig2.tight_layout()
fig2_b64 = fig_to_b64(fig2)

# ---------------------------------------------------------------------------
# Figure 3 – Tier breakdown donut
# ---------------------------------------------------------------------------
print("Generating Figure 3 (donut)...")
tier_counts = tier.value_counts().reindex(["Tier 1", "Tier 2", "Tier 3"])
fig3, ax3 = plt.subplots(figsize=(6, 5))
wedge_props = dict(width=0.5, edgecolor="white", linewidth=2)
wedges, _, autotexts = ax3.pie(
    tier_counts.values,
    labels=[
        f"Tier 1 – Rescuable\n(≥2 guides, n={tier_counts['Tier 1']})",
        f"Tier 2 – Marginal\n(1 guide, n={tier_counts['Tier 2']})",
        f"Tier 3 – Unrescuable\n(0 guides, n={tier_counts['Tier 3']})",
    ],
    colors=[TIER_COLORS["Tier 1"], TIER_COLORS["Tier 2"], TIER_COLORS["Tier 3"]],
    autopct="%1.1f%%", startangle=90, wedgeprops=wedge_props,
    textprops={"fontsize": 9},
)
for at in autotexts:
    at.set_fontsize(9)
ax3.set_title("Rescue tier breakdown\n520 consistently removed genes", fontsize=11, pad=12)
fig3.tight_layout()
fig3_b64 = fig_to_b64(fig3)

# ---------------------------------------------------------------------------
# Figure 4 – Top intronic-affected genes (stacked bar)
# ---------------------------------------------------------------------------
print("Generating Figure 4 (intronic)...")
top_intron_idx = intron_total_per_gene.sort_values(ascending=False).head(35).index
top_intron_df  = not_exon.loc[top_intron_idx]

fig4, ax4 = plt.subplots(figsize=(10, 7))
bottom = np.zeros(len(top_intron_idx))
for lib in libs:
    vals = top_intron_df[lib].values
    ax4.barh(range(len(top_intron_idx)), vals, left=bottom,
             color=LIB_COLORS[lib], label=lib, edgecolor="white", linewidth=0.4)
    bottom += vals
ax4.set_yticks(range(len(top_intron_idx)))
ax4.set_yticklabels(top_intron_idx, fontsize=8)
ax4.invert_yaxis()
ax4.set_xlabel("Intronic guides removed (total across all libraries)", fontsize=10)
ax4.set_title("Top 35 genes by intronic guide removal\n(already excluded from counted total sgRNA)", fontsize=11)
ax4.legend(title="Library", bbox_to_anchor=(1.01, 1), loc="upper left")
sns.despine()
fig4.tight_layout()
fig4_b64 = fig_to_b64(fig4)

# ---------------------------------------------------------------------------
# Figure 5 – Best rescue library per gene (Tier 1 stacked bar)
# ---------------------------------------------------------------------------
print("Generating Figure 5 (best rescue library)...")
best_lib_t1 = best_lib[has_2plus].value_counts()
fig5, ax5 = plt.subplots(figsize=(6, 4))
bars = ax5.bar(best_lib_t1.index, best_lib_t1.values,
               color=[LIB_COLORS[l] for l in best_lib_t1.index],
               edgecolor="white", linewidth=0.5)
ax5.bar_label(bars, fontsize=10, fontweight="bold")
ax5.set_ylabel("Number of Tier 1 genes", fontsize=10)
ax5.set_title("Best single rescue library per Tier 1 gene\n(library with most surviving guides)", fontsize=11)
ax5.set_ylim(0, best_lib_t1.max() + 20)
sns.despine()
fig5.tight_layout()
fig5_b64 = fig_to_b64(fig5)

# ---------------------------------------------------------------------------
# Build HTML
# ---------------------------------------------------------------------------
print("Building HTML report...")
today_str = date.today().strftime("%B %d, %Y")

# Per-library stats rows
lib_stat_rows = ""
for lib in ["TKOv3", "Jacquere", "Brunello", "Yusa", "Avana"]:
    s = counted[lib]
    lib_stat_rows += (
        f"<tr><td><strong>{lib}</strong></td>"
        f"<td>{s.mean():.2f}</td>"
        f"<td>{(s==0).sum()}</td>"
        f"<td>{(s==1).sum()}</td>"
        f"<td>{(s==2).sum()}</td>"
        f"<td>{(s==3).sum()}</td>"
        f"<td>{(s>=4).sum()}</td>"
        f"<td>{int(s.max())}</td></tr>\n"
    )

# Missing gene rows
missing_rows = ""
for lib in ["Avana", "Brunello", "TKOv3", "Yusa", "Jacquere"]:
    m = missing_genes[lib]
    genes_str = ", ".join(f"<code>{g}</code>" for g in m) if m else "<em>none</em>"
    missing_rows += f"<tr><td>{lib}</td><td>{len(m)}</td><td>{genes_str}</td></tr>\n"

# Tier 2 gene tags
tier2_tags  = gene_tags(tier2_genes)
tier3_tags  = gene_tags(tier3_genes)
intron_tags = gene_tags(sorted(zero_due_to_intron[zero_due_to_intron].index))

# Tier 1 rescue table HTML (all cols except Best Library already in index)
tier1_show = tier1_tbl[list(REPORTS) + ["Max (single lib)", "Best Library", "Total (all libs)"]].copy()
tier1_html = tier1_show.to_html(
    classes="table table-striped table-bordered table-sm",
    border=0, index=True
)

# ---------------------------------------------------------------------------
HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Rescue Report – 520 Unrepresented CRISPR-KO Genes</title>
  <link rel="stylesheet"
    href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
  <style>
    body          {{ font-family: "Segoe UI", Arial, sans-serif; background: #f4f6f9; }}
    .hero         {{ background: linear-gradient(135deg,#1a237e,#3949ab);
                    color:white; padding:3rem 2rem; margin-bottom:2rem; }}
    .section      {{ background:white; border-radius:10px; padding:2rem;
                    margin-bottom:2rem; box-shadow:0 2px 8px rgba(0,0,0,.07); }}
    .stat-box     {{ background:white; border-radius:10px; padding:1.5rem;
                    text-align:center; box-shadow:0 2px 8px rgba(0,0,0,.1); }}
    .stat-number  {{ font-size:2.5rem; font-weight:700; }}
    .t1           {{ color:#27ae60; }}
    .t2           {{ color:#f39c12; }}
    .t3           {{ color:#e74c3c; }}
    .gene-tag     {{ display:inline-block; background:#e8eaf6; color:#1a237e;
                    padding:2px 7px; border-radius:4px; margin:2px 1px;
                    font-size:.78rem; font-family:monospace; }}
    .info-box     {{ background:#e3f2fd; border-left:5px solid #2196F3;
                    padding:1rem 1.2rem; border-radius:0 8px 8px 0; margin-bottom:1rem; }}
    .warn-box     {{ background:#fff3cd; border-left:5px solid #ffc107;
                    padding:1rem 1.2rem; border-radius:0 8px 8px 0; margin-bottom:1rem; }}
    .ok-box       {{ background:#e8f5e9; border-left:5px solid #4CAF50;
                    padding:1rem 1.2rem; border-radius:0 8px 8px 0; margin-bottom:1rem; }}
    .err-box      {{ background:#fce4ec; border-left:5px solid #e74c3c;
                    padding:1rem 1.2rem; border-radius:0 8px 8px 0; margin-bottom:1rem; }}
    h2            {{ border-left:5px solid #3949ab; padding-left:12px;
                    margin-bottom:1.5rem; color:#1a237e; }}
    h3            {{ color:#3949ab; margin-top:1.5rem; }}
    .fig-caption  {{ font-style:italic; color:#666; font-size:.85rem; margin-top:.4rem; }}
    .table        {{ font-size:.82rem; }}
    footer        {{ border-top:1px solid #dee2e6; }}
  </style>
</head>
<body>

<div class="hero">
  <div class="container">
    <h1 class="display-5 fw-bold">
      Rescue Report: 520 Consistently Unrepresented Genes
    </h1>
    <p class="lead mb-1">
      Cross-library survival analysis of exon-targeting guides across
      Avana, Brunello, TKOv3, Yusa, and Jacquere
      &mdash; T2T-CHM13 reference genome
    </p>
    <small>Generated: {today_str}</small>
  </div>
</div>

<div class="container">

<!-- ── Summary stats ─────────────────────────────────────────────────── -->
<div class="row g-3 mb-4">
  <div class="col-6 col-md-3">
    <div class="stat-box">
      <div class="stat-number" style="color:#3949ab">520</div>
      <div class="text-muted small">Genes removed from<br>all 5 libraries (intersection)</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-box">
      <div class="stat-number t1">{has_2plus.sum()}</div>
      <div class="text-muted small">Tier 1 – Rescuable<br>(&ge;2 guides in &ge;1 library)</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-box">
      <div class="stat-number t2">{has_1only.sum()}</div>
      <div class="text-muted small">Tier 2 – Marginal<br>(exactly 1 guide, any library)</div>
    </div>
  </div>
  <div class="col-6 col-md-3">
    <div class="stat-box">
      <div class="stat-number t3">74</div>
      <div class="text-muted small">Tier 3 – Unrescuable<br>(0 guides in all 5 libraries)</div>
    </div>
  </div>
</div>

<!-- ── Section 1: Background ─────────────────────────────────────────── -->
<div class="section">
  <h2>1. Background &amp; Methodology</h2>
  <p>
    These 520 protein-coding genes are absent from all five widely-used CRISPR-KO
    libraries when their sgRNAs are evaluated against the T2T-CHM13 reference genome
    using the <strong>guiderefine</strong> pipeline. This report quantifies how many
    <em>exon-targeting, specific</em> guides remain available per gene in each library,
    and proposes a tiered rescue strategy.
  </p>
  <div class="info-box">
    <strong>Guide counting method:</strong> All counts are taken from the
    <code>counted total sgRNA</code> column of each library's
    <code>full_report.xlsx</code> (Genes Report sheet). This column already
    subtracts <em>all</em> disposal categories:
    <ul class="mb-0 mt-1">
      <li>Multi-target guides</li>
      <li>Single-mismatch off-target guides</li>
      <li>PAM-distal double-mismatch guides</li>
      <li>Guides not targeting anywhere in the T2T genome</li>
      <li><strong>Guides not targeting an exon (intronic guides)</strong></li>
    </ul>
    Values are summed across duplicate gene symbols (isoforms with the same HGNC symbol).
    Genes absent from a library report are set to 0.
  </div>
</div>

<!-- ── Section 2: Per-library performance ────────────────────────────── -->
<div class="section">
  <h2>2. Per-Library Surviving Guide Summary</h2>
  <div class="row align-items-center">
    <div class="col-lg-7">
      <table class="table table-bordered table-sm">
        <thead class="table-dark">
          <tr>
            <th>Library</th><th>Mean</th>
            <th>0 guides</th><th>1 guide</th><th>2 guides</th>
            <th>3 guides</th><th>&ge;4 guides</th><th>Max</th>
          </tr>
        </thead>
        <tbody>
          {lib_stat_rows}
        </tbody>
      </table>
      <p class="text-muted small">
        <strong>TKOv3</strong> is the strongest source: highest mean (1.46),
        fewest zero-guide genes (137), and the only library reaching 4 surviving guides.
        It is the primary rescue library for 174 of the 342 Tier 1 genes.
      </p>
    </div>
    <div class="col-lg-5">
      <img src="data:image/png;base64,{fig3_b64}" class="img-fluid" alt="Tier donut">
      <p class="fig-caption text-center">Figure 1. Rescue tier breakdown.</p>
    </div>
  </div>

  <img src="data:image/png;base64,{fig1_b64}" class="img-fluid mt-3" alt="Guide distribution">
  <p class="fig-caption">
    Figure 2. Number of the 520 genes in each guide-count bucket per library.
    The red bars (0 guides) dominate in Avana and Yusa.
    TKOv3 and Jacquere have the most genes with &ge;2 guides.
  </p>
</div>

<!-- ── Section 3: Best rescue library ───────────────────────────────── -->
<div class="section">
  <h2>3. Best Rescue Library Per Gene (Tier 1)</h2>
  <div class="row align-items-center">
    <div class="col-md-5">
      <img src="data:image/png;base64,{fig5_b64}" class="img-fluid" alt="Best library">
      <p class="fig-caption">
        Figure 3. Library providing the highest guide count
        for each of the {has_2plus.sum()} Tier 1 genes.
      </p>
    </div>
    <div class="col-md-7">
      <div class="ok-box">
        <strong>TKOv3</strong> is the best single source for
        <strong>174 / {has_2plus.sum()} Tier 1 genes (51%).</strong><br>
        Brunello covers a further 66, Avana 49, Jacquere 33, Yusa 20.
      </div>
      <p class="mt-2">
        For genes where multiple libraries offer surviving guides, combining guides
        from the two best sources is strongly recommended to reach the
        target of <strong>&ge;4 independent guides per gene</strong>.
      </p>
    </div>
  </div>
</div>

<!-- ── Section 4: Heatmap ────────────────────────────────────────────── -->
<div class="section">
  <h2>4. Per-Gene Surviving Guide Heatmap</h2>
  <p>
    All 520 genes sorted by tier (Tier 1 top), then by total guides descending.
    Dashed lines mark tier boundaries.
    This panel is suitable as a <strong>supplementary figure</strong> in the manuscript.
  </p>
  <div class="text-center">
    <img src="data:image/png;base64,{fig2_b64}" class="img-fluid"
         style="max-height:950px" alt="Heatmap">
  </div>
  <p class="fig-caption text-center">
    Figure 4. Heatmap of surviving exon-targeting guides.
    White = 0, light blue = 1, medium blue = 2, green = 3, dark green = 4.
    Dashed lines separate tiers.
  </p>
</div>

<!-- ── Section 5: Intronic guides ───────────────────────────────────── -->
<div class="section">
  <h2>5. Intronic Guide Impact</h2>
  <div class="ok-box">
    <strong>Intronic guides are already excluded</strong> from
    <code>counted total sgRNA</code>. The heatmap above reflects
    <em>exon-targeting guides only</em>. No further filtering is needed
    for the rescue analysis.
  </div>
  <p>
    Nevertheless, intronic removal is the <em>direct cause</em> of zero surviving guides
    for <strong>{zero_due_to_intron.sum()} genes</strong> &mdash; genes that had original
    library design counts but whose guides all fell in introns when mapped to T2T-CHM13.
    These genes appear in Tier 3 not because no guide was ever designed, but because
    none of the existing designs target a coding exon in the T2T reference.
    <strong>New guide design targeting annotated CDS regions is required for these genes.</strong>
  </p>
  <p><strong>Intronic-driven Tier 3 genes ({zero_due_to_intron.sum()}):</strong></p>
  <div class="mb-3">{intron_tags}</div>

  <p>
    <strong>{n_intron_affected} of 520 genes (23%)</strong> had at least one intronic
    guide removed across any library. Total removed: Avana&nbsp;{int(not_exon['Avana'].sum())}
    &middot; Brunello&nbsp;{int(not_exon['Brunello'].sum())}
    &middot; TKOv3&nbsp;{int(not_exon['TKOv3'].sum())}
    &middot; Yusa&nbsp;{int(not_exon['Yusa'].sum())}
    &middot; Jacquere&nbsp;{int(not_exon['Jacquere'].sum())}.
  </p>
  <img src="data:image/png;base64,{fig4_b64}" class="img-fluid" alt="Intronic bar chart">
  <p class="fig-caption">
    Figure 5. Top 35 genes ranked by total intronic guides removed across all libraries.
    DMPK, CLDN5, CLN3, FMO2, SLAIN1, ATXN3 had the most original guides lost to intron targeting.
  </p>
</div>

<!-- ── Section 6: Tier 3 genes ──────────────────────────────────────── -->
<div class="section">
  <h2>6. Tier 3 &mdash; 74 Completely Unrescuable Genes</h2>
  <div class="err-box">
    These 74 genes have <strong>0 surviving exon-targeting guides in all 5 libraries.</strong>
    They cannot be rescued by pooling from existing libraries.
    New guide design using the pre-computed spacer CDS data is required.
  </div>
  <div class="mb-3">{tier3_tags}</div>
  <div class="warn-box">
    <strong>Suggested action:</strong> Use
    <code>all_spacers_cds_data_499_removed_genes.csv</code> to identify candidate
    spacers within the CDS of these genes. Run all candidates through the
    guiderefine pipeline (T2T-CHM13 alignment + off-target scoring)
    before synthesis and validation.
  </div>
</div>

<!-- ── Section 7: Tier 2 genes ──────────────────────────────────────── -->
<div class="section">
  <h2>7. Tier 2 &mdash; {has_1only.sum()} Marginal Genes (1 surviving guide)</h2>
  <div class="warn-box">
    Each of these genes has exactly <strong>one</strong> surviving exon-targeting guide
    in at least one library. A single guide is insufficient for robust KO validation
    (off-target effects cannot be separated from on-target phenotype).
    <strong>Supplement with at least one newly designed guide per gene.</strong>
  </div>
  <div class="mb-3">{tier2_tags}</div>
</div>

<!-- ── Section 8: Tier 1 rescue table ──────────────────────────────── -->
<div class="section">
  <h2>8. Tier 1 &mdash; {has_2plus.sum()} Rescuable Genes (&ge;2 surviving guides)</h2>
  <p>
    Sorted by total guides available across all libraries (descending).
    The <em>Best Library</em> column indicates the primary source.
    Where multiple libraries contribute guides, combine them to reach &ge;4 per gene.
  </p>
  {tier1_html}
  <p class="fig-caption small">
    Counts shown are exon-targeting surviving guides from each library's full report.
  </p>
</div>

<!-- ── Section 9: Gene symbol issues ───────────────────────────────── -->
<div class="section">
  <h2>9. Gene Symbol Annotation Gaps</h2>
  <p>
    The following genes were absent from the guiderefine Genes Report for specific
    libraries, most likely because the library was designed using deprecated HGNC symbols.
    They default to 0 in this analysis.
    <strong>Jacquere is the only library with full coverage of all 520 genes.</strong>
  </p>
  <table class="table table-bordered table-sm">
    <thead class="table-dark">
      <tr><th>Library</th><th>Missing</th><th>Gene symbols</th></tr>
    </thead>
    <tbody>{missing_rows}</tbody>
  </table>
  <div class="warn-box">
    Before concluding that these genes have 0 guides, cross-reference the original
    library files using sequence-level matching against the current CDS annotation.
    The guides may exist under an older gene symbol.
  </div>
</div>

<!-- ── Section 10: Research figure suggestions ──────────────────────── -->
<div class="section">
  <h2>10. Suggested Research Figure Layout</h2>
  <p>
    For the manuscript figure illustrating the cross-library rescue strategy,
    the following panel arrangement is recommended:
  </p>

  <h3>Panel A &mdash; Guide scarcity motivation</h3>
  <div class="ok-box">
    Adapt <strong>Figure 2</strong> (grouped bar chart) to show the proportion of 520 genes
    with 0/1/&ge;2 surviving guides per library. This motivates the cross-library
    approach by showing that no single library is sufficient.
  </div>

  <h3>Panel B &mdash; Rescue tier overview</h3>
  <div class="ok-box">
    Use <strong>Figure 1</strong> (donut chart) with the three-tier colour code
    (green / amber / red) to summarise the rescue landscape at a glance.
  </div>

  <h3>Panel C &mdash; Best rescue source</h3>
  <div class="ok-box">
    Use <strong>Figure 3</strong> (stacked bar) to show that TKOv3 is the dominant
    rescue source (51% of Tier 1 genes), with Brunello and Avana as secondary sources.
  </div>

  <h3>Panel D &mdash; Per-gene heatmap (supplementary)</h3>
  <div class="ok-box">
    <strong>Figure 4</strong> (heatmap) as a supplementary figure.
    For the main figure, consider showing only Tier 1 genes (&le;342 rows)
    to improve legibility.
    Label a handful of biologically notable genes (e.g., ATXN3, DMPK, CLN3, FOXO3).
  </div>

  <h3>Proposed rescue workflow (for methods / figure legend)</h3>
  <ol>
    <li>
      <strong>Tier 1 (&ge;2 guides in &ge;1 library):</strong>
      Use the library with the highest surviving guide count as primary source.
      Combine guides from multiple libraries to reach a minimum of 4 independent
      guides per gene. TKOv3 is the recommended first choice.
    </li>
    <li>
      <strong>Tier 2 (exactly 1 guide):</strong>
      Retain the single available guide. Design at least one additional guide
      <em>de novo</em> from the CDS spacer dataset and validate both independently.
    </li>
    <li>
      <strong>Tier 3 (0 guides in all libraries):</strong>
      Design all guides <em>de novo</em> using
      <code>all_spacers_cds_data_499_removed_genes.csv</code>.
      Apply the full guiderefine pipeline (T2T-CHM13 alignment,
      off-target scoring, exon verification) before synthesis.
      Note that 51 of these 74 genes are in Tier 3 specifically because their
      existing library guides targeted introns; they may be easier to recover
      with new T2T-aware designs.
    </li>
  </ol>
</div>

</div><!-- /container -->

<footer class="text-center text-muted py-4 mt-2">
  <small>
    Generated from guiderefine full reports &middot; T2T-CHM13 reference
    &middot; {today_str}
  </small>
</footer>
</body>
</html>"""

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
with open(OUT_PATH, "w", encoding="utf-8") as fh:
    fh.write(HTML)

print(f"\nReport saved to: {os.path.abspath(OUT_PATH)}")
