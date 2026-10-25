#!/usr/bin/env python3
"""
Figure 2 (extension) — frozen cross-sectional evidence.

Panels (analysis_extension_plan.md WP2 / manuscript_v4 storyboard):
  A  TAL programme schematic: nephron position + six frozen genes grouped by role
  B  Six-gene expression context across epithelial cell types (descriptive)
  C  GSE183276 donor-level TAL scores (n=10 AKI, n=17 Reference)
  D  GSE210622 donor-level TAL scores (n=8 AKI, n=4 non-AKI; mixed acquisition)
  E  Cross-sectional effect forest, evidence blocks separated (no pooled P)
  F  Six-gene direction + LODO stability summary

All numbers are read verbatim from the frozen results/ CSVs and are never
recomputed here. Statistical units: donor (C-D), sample/donor (E, tissue
blocks). No cell-level P values, no pooled cohort P value.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, patches

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(BASE, "results")
EXT = os.path.join(BASE, "results", "extension")
OUT = os.path.join(BASE, "figures", "extension")

CTRL = "#2166AC"   # blue = control
AKI  = "#B2182B"   # red  = AKI
CTRL_L = "#D1E5F0"
AKI_L  = "#FDDBC7"
GREY = "#444444"
TAL_C = "#D95F02"
GOLD = "#E0A800"

TAL_GENES = ["SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10"]

def panel_letter(ax, letter):
    ax.text(-0.05, 1.03, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="right", va="bottom")

# ------------------------------------------------------------------ data ---
ev = pd.read_csv(os.path.join(RES, "formal_evidence.csv"))
fx = pd.read_csv(os.path.join(RES, "formal_effects.csv"))
fg = pd.read_csv(os.path.join(RES, "formal_genes.csv"))
bs = pd.read_csv(os.path.join(RES, "bulk_scores.csv"))
be = pd.read_csv(os.path.join(RES, "bulk_effects.csv"))
bg = pd.read_csv(os.path.join(RES, "bulk_genes.csv"))

prim = ev[ev["tier"] == "primary_ge50"].copy()
DATA_CTRL = {"GSE183276": "Reference", "GSE210622": "Reference"}
CTRL_NAME = {"GSE183276": "Reference", "GSE210622": "non-AKI"}
FACET_TITLE = {
    "GSE183276": "GSE183276\nAKI (n = 10) vs Reference (n = 17)",
    "GSE210622": "GSE210622 \u2014 postmortem AKI vs mixed-source non-AKI\nAKI (n = 8) vs non-AKI (n = 4)",
}
cohorts = {}
for ds in ["GSE183276", "GSE210622"]:
    sub = prim[(prim["dataset"] == ds) & (prim["group"].isin(["AKI", DATA_CTRL[ds]]))].copy()
    sub["group"] = sub["group"].replace({DATA_CTRL[ds]: "control"})
    cohorts[ds] = sub

eff_cell = fx[(fx["tier"] == "primary_ge50") & (fx["model"] == "unadjusted") &
              (fx["term"] == "AKI") &
              (fx["comparison"].isin(["AKI_vs_Reference", "AKI_vs_nonAKI"]))].copy()

eff_bulk = be[(be["model"] == "unadjusted") & (be["analysis_level"] == "primary")].copy()

fg_prim = fg[fg["tier"] == "primary_ge50"].copy()
# Panel F shows AKI vs the control group only (GSE183276 AKI_vs_Reference,
# GSE210622 AKI_vs_nonAKI); the GSE183276 AKI_vs_CKD rows are excluded so CKD
# is never presented as a control.
fg_prim = fg_prim[fg_prim["comparison"].isin(["AKI_vs_Reference", "AKI_vs_nonAKI"])]
bg_prim = bg[bg["analysis_level"] == "primary"].copy()

# ------------------------------------------------------------------ figure --
fig = plt.figure(figsize=(11.7, 8.6))
gs = fig.add_gridspec(3, 4, height_ratios=[0.95, 1.05, 0.75], hspace=0.55,
                      wspace=0.42, left=0.07, right=0.97, top=0.90, bottom=0.10)
axA = fig.add_subplot(gs[0, 0:2])
axB = fig.add_subplot(gs[0, 2:4])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[1, 2:4])
axF = fig.add_subplot(gs[2, :])

# --------------------------------------------------------------- Panel A ---
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
# Stylised nephron: glomerulus, PT, descending/ascending loop (TAL), DCT, CD
axA.add_patch(patches.Circle((1.6, 8.6), 0.55, facecolor="#FDDBC7", edgecolor=AKI, lw=1.2))
axA.text(1.6, 8.6, "G", ha="center", va="center", fontsize=7)
axA.plot([1.6, 1.6], [7.9, 7.0], color="#666666", lw=2)
axA.text(1.05, 7.45, "PT", fontsize=7, color="#1B9E77", fontweight="bold")
axA.plot([1.6, 3.4], [7.0, 5.4], color="#666666", lw=2)
axA.plot([3.4, 3.4], [5.4, 2.6], color="#666666", lw=2)         # descending
axA.plot([3.4, 5.6], [2.6, 2.6], color=TAL_C, lw=3.5)          # loop bend = TAL
axA.plot([5.6, 5.6], [2.6, 5.4], color=TAL_C, lw=3.5)          # ascending = TAL
axA.plot([5.6, 7.4], [5.4, 6.4], color="#666666", lw=2)
axA.text(6.9, 6.7, "DCT", fontsize=7, color="#7570B3", fontweight="bold")
axA.plot([7.4, 8.6], [6.4, 6.0], color="#666666", lw=2)
axA.text(8.5, 5.6, "CD", fontsize=7, color="#66A61E", fontweight="bold")
axA.text(4.5, 1.95, "TAL (thick ascending limb)", fontsize=8, color=TAL_C,
         fontweight="bold", ha="center")
axA.text(4.5, 1.30, "Frozen programme: transport + differentiation identity",
         fontsize=7, ha="center", color=GREY)
groups = [("Transport", ["SLC12A1", "KCNJ10"]),
          ("Paracellular", ["CLDN16", "CLDN10"]),
          ("Differentiation", ["UMOD", "PTH1R"])]
yy = [1.05, 0.45, -0.15]
for (gname, genes), y in zip(groups, yy):
    axA.text(4.5, 4.7 + y, f"{gname}:  " + "  \u00b7  ".join(genes),
             fontsize=6.8, ha="center", family="monospace",
             bbox=dict(boxstyle="round,pad=0.25", fc="#F7F7F7", ec="#999999", lw=0.6))
panel_letter(axA, "A")

# --------------------------------------------------------------- Panel B ---
mk = pd.concat([pd.read_csv(os.path.join(EXT, f"atlas_{ds}_markers.csv"))
                for ds in ["GSE183276", "GSE210622"]])
mk6 = mk[mk["gene"].isin(TAL_GENES)].copy()
mk6["mean_expr_z"] = mk6.groupby(["cohort", "gene"])["mean_expr"].transform(
    lambda x: (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) > 0 else 0.0)
mk6["gene"] = pd.Categorical(mk6["gene"], categories=TAL_GENES, ordered=True)
mk6["broad_class"] = pd.Categorical(mk6["broad_class"],
                                    categories=["PT","TAL/LOH","DCT","CNT","CD"],
                                    ordered=True)
mk6 = mk6[mk6["broad_class"].notna() & mk6["gene"].notna()].copy()
class_order = ["PT", "TAL/LOH", "DCT", "CNT", "CD"]
cohort_order = ["GSE183276", "GSE210622"]
x_map = {(cohort, cls): i * 6 + j for i, cohort in enumerate(cohort_order)
         for j, cls in enumerate(class_order)}
mk6["x"] = [x_map[(c, str(k))] for c, k in zip(mk6["cohort"], mk6["broad_class"])]
sc = axB.scatter(mk6["x"], mk6["gene"],
                 s=(mk6["pct_expr"] * 220) + 2,
                 c=mk6["mean_expr_z"], cmap="RdBu_r", vmin=-2.0, vmax=2.0,
                 edgecolor="#666666", linewidth=0.3)
axB.axvline(5, color="#999999", lw=0.7)
axB.set_xticks([x_map[(c, k)] for c in cohort_order for k in class_order])
axB.set_xticklabels([f"{c.replace('GSE','')}\n{k}" for c in cohort_order for k in class_order],
                    rotation=45, ha="right", fontsize=5.2)
axB.axhline(2.5, color="#CCCCCC", lw=0.6)
axB.axhline(4.5, color="#CCCCCC", lw=0.6)
cb = fig.colorbar(sc, ax=axB, fraction=0.045, pad=0.02)
cb.set_label("mean expression (within-gene z)", fontsize=6)
cb.ax.tick_params(labelsize=5.5)
axB.set_title("Six-gene expression context across epithelial classes\n"
              "GSE183276 + GSE210622 (descriptive)", loc="left", fontsize=7.5)
axB.set_xlabel("Cohort-specific broad cell class")
axB.tick_params(axis="x", labelsize=6)
panel_letter(axB, "B")

# ------------------------------------------------------------- Panel C/D ---
rng = np.random.default_rng(20250730)
for ax, ds in [(axC, "GSE183276"), (axD, "GSE210622")]:
    sub = cohorts[ds]
    ctrl = sub.loc[sub["group"] == "control", "TAL_identity"]
    aki  = sub.loc[sub["group"] == "AKI", "TAL_identity"]
    boxes = ax.boxplot([ctrl, aki], positions=[1, 2], widths=0.52, patch_artist=True,
                       showfliers=False, zorder=2,
                       medianprops=dict(color="#222222", lw=1.1),
                       whiskerprops=dict(color="#555555", lw=0.9),
                       capprops=dict(color="#555555", lw=0.9),
                       boxprops=dict(lw=0.9))
    for patch, col in zip(boxes["boxes"], [CTRL, AKI]):
        patch.set_facecolor(CTRL_L if col == CTRL else AKI_L)
        patch.set_edgecolor(col)
    for vals, col, pos in [(ctrl, CTRL, 1), (aki, AKI, 2)]:
        jit = rng.uniform(-0.14, 0.14, size=len(vals))
        ax.scatter(pos + jit, vals, s=28, color=col, edgecolor="white",
                   linewidth=0.6, zorder=3, alpha=0.95)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"{CTRL_NAME[ds]}\n(n = {len(ctrl)})", f"AKI\n(n = {len(aki)})"])
    ax.set_xlim(0.4, 2.6)
    ax.set_ylim(-1.5, 1.5)
    ax.tick_params(axis="x", length=0)
    if ds == "GSE183276":
        ax.set_ylabel("Six-gene TAL detectability score\n(donor pseudobulk)")
    ax.set_title(FACET_TITLE[ds], loc="left", fontsize=6.8)
    if ds == "GSE210622":
        ax.text(0.5, -0.32, "mixed acquisition: 8 postmortem AKI\nvs 3 nephrectomy + 1 postmortem non-AKI",
                transform=ax.transAxes, fontsize=5.6, color=GREY, ha="left", va="top")
panel_letter(axC, "C")
panel_letter(axD, "D")

# --------------------------------------------------------------- Panel E ---
# Two evidence blocks: cell-resolved (donor) and tissue-level (sample/donor)
ec = eff_cell.sort_values("n_aki", ascending=False).reset_index(drop=True)
eb = eff_bulk.copy()
# row positions: cell-resolved block rows 2.0/1.0, tissue block rows 0.5/-0.5
rows_cell = np.arange(len(ec))[::-1] * 1.0 + 1.0   # [1.0, 2.0] for 2 rows -> top at 2.0
rows_bulk = np.arange(len(eb))[::-1] * 1.0 + 0.5 - 0.5  # [0.0, -1.0] for 2 rows -> top at 0.0
# make the cell block sit above the tissue block with a clear gap
rows_cell = rows_cell + 1.0  # [2.0, 3.0]
ylab, ypos = [], []
for y, r in zip(rows_cell, ec.itertuples(index=False)):
    axE.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]],
                 fmt="s", ms=6, color=AKI, mec="white", mew=0.8, lw=1.4,
                 capsize=3, capthick=1.2, zorder=4)
    est = f"{r.estimate:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]".replace("-", "\u2212")
    p = r.p_value
    pstr = f"P = {p:.1e}" if p < 0.001 else f"P = {p:.3f}"
    axE.text(0.30, y, f"{est}\n{pstr}", fontsize=6, va="center", ha="left",
             color="#222222", linespacing=1.3)
    ylab.append(r.dataset); ypos.append(y)
for y, r in zip(rows_bulk, eb.itertuples(index=False)):
    axE.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]],
                 fmt="o", ms=6, color=GOLD, mec="white", mew=0.8, lw=1.4,
                 capsize=3, capthick=1.2, zorder=4)
    est = f"{r.estimate:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]".replace("-", "\u2212")
    p = r.p_value
    pstr = f"P = {p:.1e}" if p < 0.001 else f"P = {p:.3f}"
    axE.text(0.30, y, f"{est}\n{pstr}", fontsize=6, va="center", ha="left",
             color="#222222", linespacing=1.3)
    ylab.append(r.dataset); ypos.append(y)
axE.axvline(0, color="#888888", ls="--", lw=0.9, zorder=1)
axE.axhspan(min(rows_cell) - 0.5, max(rows_cell) + 0.5, facecolor="#FBF0F0", zorder=0)
axE.axhspan(min(rows_bulk) - 0.5, max(rows_bulk) + 0.5, facecolor="#FBF6E8", zorder=0)
axE.set_yticks(ypos)
axE.set_yticklabels(ylab)
axE.set_ylim(min(rows_bulk) - 1.0, max(rows_cell) + 1.0)
axE.set_xlim(-2.6, 3.2)
axE.set_xticks([-2, -1, 0, 1, 2, 3])
axE.set_xlabel("AKI minus control, \u03b2 [95% CI]")
axE.tick_params(axis="y", length=0)
axE.text(0.02, max(rows_cell) + 0.62, "cell-resolved  \u2014  donor unit",
         transform=axE.get_yaxis_transform(), fontsize=6.4, color=GREY, ha="left")
axE.text(0.02, max(rows_bulk) - 0.55, "tissue-level  \u2014  sample/donor unit",
         transform=axE.get_yaxis_transform(), fontsize=6.4, color=GREY, ha="left")
axE.set_title("Cross-sectional effects, separated evidence blocks (no pooled P)",
              loc="left", fontsize=7.5)
axE.text(0.99, -0.18,
         "Displayed P values are model-based; permutation P is retained in the freeze receipt",
         transform=axE.transAxes, fontsize=5.6, color=GREY, ha="right", va="top")
panel_letter(axE, "E")

# --------------------------------------------------------------- Panel F ---
pivot = fg_prim.pivot_table(index="gene", columns="dataset", values="z_effect")
pivot = pivot.reindex(TAL_GENES)
pivot = pivot[["GSE183276", "GSE210622"]]
im = axF.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
axF.set_xticks(range(2))
axF.set_xticklabels(["GSE183276", "GSE210622"])
axF.set_yticks(range(6))
axF.set_yticklabels(TAL_GENES)
for i in range(6):
    for j in range(2):
        v = pivot.values[i, j]
        axF.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6.2,
                 color="white" if abs(v) > 1 else "#222222")
axF.set_title("Per-gene z-effect (AKI vs control), donor pseudobulk\n"
              "LODO negative: 100% both cohorts; directions 6/6 vs 5/6",
              loc="left", fontsize=7.5)
cb2 = fig.colorbar(im, ax=axF, fraction=0.025, pad=0.02)
cb2.set_label("z-effect", fontsize=6)
cb2.ax.tick_params(labelsize=5.5)
# Annotate bulk per-gene direction counts at the right
bulk_note = ""
for ds in ["GSE139061", "GSE30718"]:
    d = bg_prim[bg_prim["dataset"] == ds]
    neg = (d["direction"] == "negative").sum()
    bulk_note += f"{ds}: {neg}/{len(d)} genes negative  "
axF.text(0.5, -0.25, bulk_note, transform=axF.transAxes, fontsize=6.2, color=GREY,
         ha="center", va="bottom")
panel_letter(axF, "F")

fig.suptitle("Figure 2  |  Frozen TAL programme across independent cross-sectional cohorts",
             fontsize=11, fontweight="bold", y=0.975)
fig.text(0.5, 0.008,
         "Expression-based programme detectability is not a TAL abundance, TAL-specificity, or causal estimate.",
         ha="center", fontsize=6.2, color=GREY)
fig.savefig(os.path.join(OUT, "figure2_frozen.png"), dpi=300)
fig.savefig(os.path.join(OUT, "figure2_frozen.pdf"))
print("wrote figure2_frozen.png/.pdf to", OUT)
