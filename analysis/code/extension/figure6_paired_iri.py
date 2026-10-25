#!/usr/bin/env python3
"""
Figure 6 (extension) — paired transplant-IRI validation (GSE126805).

Panels (analysis_extension_plan.md WP6 / manuscript_v4 storyboard):
  A  sampling design schematic (pre/post-reperfusion, 3 months, 1 year;
     42 grafts, 41 primary pre/post pairs)
  B  paired pre/post trajectories (each line = one kidney allograft)
  C  paired effect summary (primary + exploratory contrasts, Holm P)
  D  six-gene paired effects + LOGO/LODO stability
  E  GSE43974 replication: DEFERRED (expression matrix not present locally)

Statistical unit: paired kidney allograft. No unpaired test, no donor pooling.
All numbers are read verbatim from the frozen results CSVs.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, patches

rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
})

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
META = os.path.join(BASE, "metadata", "gse126805_samples.csv")
RES = os.path.join(BASE, "results", "gse126805_results.csv")
OUT = os.path.join(BASE, "figures", "extension")

GREY = "#444444"; MAIN = "#1B6CA8"; EXP = "#7FBF7B"; GOLD = "#E0A800"

def panel_letter(ax, letter):
    ax.text(-0.05, 1.02, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="right", va="bottom")

samples = pd.read_csv(META)
# in_primary_pair arrives as the strings 'TRUE'/'FALSE'; normalise explicitly so
# the filter never silently drops all rows if dtype inference changes.
samples["in_primary_pair"] = samples["in_primary_pair"].astype(str).str.strip().str.upper() == "TRUE"
res = pd.read_csv(RES)

def val(label, field="estimate"):
    x = res[(res["label"] == label)][field]
    assert len(x) == 1, label
    return x.iloc[0]

# paired pre/post
paired = samples[(samples["in_primary_pair"]) &
                 (samples["timepoint"].isin(["baseline_pre", "baseline_post"]))]
pre = paired[paired["timepoint"] == "baseline_pre"][["kidney", "TAL_score"]].rename(
    columns={"TAL_score": "pre"})
post = paired[paired["timepoint"] == "baseline_post"][["kidney", "TAL_score"]].rename(
    columns={"TAL_score": "post"})
pairs = pre.merge(post, on="kidney", sort=True).reset_index(drop=True)
assert len(pairs) == 41

# per-gene deltas
pg = res[res["result_type"] == "per_gene"].copy()
pg = pg.sort_values("label").reset_index(drop=True)
logo = res[(res["result_type"] == "logo_detail")].copy()
logo = logo.sort_values("label").reset_index(drop=True)

# ------------------------------------------------------------------ figure --
fig = plt.figure(figsize=(11.7, 8.6))
gs = fig.add_gridspec(3, 4, height_ratios=[0.75, 1.15, 1.1], hspace=0.62, wspace=0.45,
                      left=0.08, right=0.97, top=0.93, bottom=0.07)
axA = fig.add_subplot(gs[0, 0:2])
axB = fig.add_subplot(gs[1, 0:2])
axC = fig.add_subplot(gs[1, 2:4])
axD = fig.add_subplot(gs[2, 0:3])
axE = fig.add_subplot(gs[2, 3])

# --------------------------------------------------------------- Panel A ---
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
axA.set_xticks([]); axA.set_yticks([])   # axis('off') still draws tick labels; clear them
tp_x = [1.5, 4.0, 6.5, 9.0]
tp_lab = ["Before\nreperfusion", "After\nreperfusion", "3 months", "1 year"]
for x, lab in zip(tp_x, tp_lab):
    axA.add_patch(patches.FancyBboxPatch((x - 1.0, 4.2), 2.0, 2.4,
                                         boxstyle="round,pad=0.12",
                                         facecolor="#EAF2F8", edgecolor=MAIN, lw=1.0))
    axA.text(x, 5.4, lab, ha="center", va="center", fontsize=7)
    if x != tp_x[-1]:
        axA.annotate("", xy=(x + 1.1, 5.4), xytext=(x + 1.0, 5.4),
                     arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.0))
axA.add_patch(patches.FancyBboxPatch((0.5, 1.2), 4.5, 1.6,
                                     boxstyle="round,pad=0.12",
                                     facecolor="#FDECEA", edgecolor="#B2182B", lw=1.2))
axA.text(2.75, 2.0, "PRIMARY paired contrast: after minus before", ha="center", va="center", fontsize=6.6,
         color="#B2182B", fontweight="bold")
axA.text(7.2, 1.55, "41 pre/post pairs\n42 grafts · 163 protocol biopsies",
         ha="center", va="center", fontsize=6.6, color=GREY)
axA.text(5, 9.3, "GSE126805 — transplant ischaemia-reperfusion protocol biopsies (bulk RNA-seq)",
         ha="center", fontsize=8, fontweight="bold")
axA.text(5, 8.4, "Statistical unit: paired kidney allograft (repeated measures, not independent replicates)",
         ha="center", fontsize=6.6, color=GREY)
panel_letter(axA, "A")

# --------------------------------------------------------------- Panel B ---
axB.set_title(f"Paired acute change (n = {len(pairs)} allografts)", loc="left", fontsize=8)
for i, r in pairs.iterrows():
    axB.plot([1, 2], [r["pre"], r["post"]], color="#B3B3B3", lw=0.7, zorder=1)
axB.scatter(np.ones(len(pairs)), pairs["pre"], s=14, color="#404040", zorder=2)
axB.scatter(np.full(len(pairs), 2), pairs["post"], s=14, color="#404040", zorder=2)
axB.plot([0.88, 1.12], [pairs["pre"].mean()] * 2, lw=3, color=MAIN, zorder=3)
axB.plot([1.88, 2.12], [pairs["post"].mean()] * 2, lw=3, color=MAIN, zorder=3)
axB.set_xticks([1, 2])
axB.set_xticklabels(["Before reperfusion", "After reperfusion"])
axB.set_ylabel("Frozen TAL programme score")
delta = val("baseline_post_vs_pre_delta")
axB.text(1.5, axB.get_ylim()[0], "",
         ha="center", fontsize=6.6)
axB.text(0.02, 0.02,
         f"Mean paired change = {delta:.3f} "
         f"[{val('baseline_post_vs_pre_delta','ci_low'):.3f}, "
         f"{val('baseline_post_vs_pre_delta','ci_high'):.3f}], "
         f"P = {val('baseline_post_vs_pre_delta','p_value'):.4f}",
         transform=axB.transAxes, fontsize=7, color=GREY, va="bottom")
panel_letter(axB, "B")

# --------------------------------------------------------------- Panel C ---
forest = pd.DataFrame({
    "label": ["Post - pre (primary)", "3 months - post (exploratory)",
              "1 year - post (exploratory)", "1 year - pre (exploratory)"],
    "est": [val("baseline_post_vs_pre_delta"),
            val("3mo_vs_post"), val("1y_vs_post"), val("1y_vs_pre")],
    "lo": [val("baseline_post_vs_pre_delta", "ci_low"),
           val("3mo_vs_post", "ci_low"), val("1y_vs_post", "ci_low"),
           val("1y_vs_pre", "ci_low")],
    "hi": [val("baseline_post_vs_pre_delta", "ci_high"),
           val("3mo_vs_post", "ci_high"), val("1y_vs_post", "ci_high"),
           val("1y_vs_pre", "ci_high")],
    "n": [val("baseline_post_vs_pre_delta", "n"), val("3mo_vs_post", "n"),
          val("1y_vs_post", "n"), val("1y_vs_pre", "n")],
})
p_labels = [f"P = {val('baseline_post_vs_pre_delta','p_value'):.4f}",
            f"Holm P = {val('3mo_vs_post_holm_P'):.4f}",
            f"Holm P = {val('1y_vs_post_holm_P'):.4f}",
            f"Holm P = {val('1y_vs_pre_holm_P'):.4f}"]
yy = np.arange(len(forest))[::-1]
for y, r, pl in zip(yy, forest.itertuples(index=False), p_labels):
    col = MAIN if r.label.startswith("Post") else EXP
    axC.errorbar(r.est, y, xerr=[[r.est - r.lo], [r.hi - r.est]],
                 fmt="o", ms=7, color=col, mec="white", mew=0.8,
                 lw=1.5, capsize=3, capthick=1.2)
    axC.text(1.05, y, f"{pl}; n={int(r.n)}", fontsize=6.2, va="center", ha="left", color=GREY)
axC.axvline(0, color="#888888", ls="--", lw=0.9)
axC.set_yticks(yy)
axC.set_yticklabels(forest["label"], fontsize=6.6)
axC.set_ylim(-0.7, 3.7)
axC.set_xlim(-1.0, 2.0)
axC.set_xlabel("Mean paired score change [95% CI]")
axC.set_title("Paired allograft contrasts", loc="left", fontsize=8)
axC.text(0.02, 0.02, "Only post - pre is primary; later contrasts are exploratory",
         transform=axC.transAxes, fontsize=6.2, color=GREY, va="bottom")
panel_letter(axC, "C")

# --------------------------------------------------------------- Panel D ---
genes = pg["label"].tolist()
gdelta = pg["estimate"].tolist()
ylab = [f"{g}  {d:+.3f}" for g, d in zip(genes, gdelta)]
bar_cols = [MAIN if d < 0 else "#C0504D" for d in gdelta]
axD.barh(np.arange(len(genes))[::-1], gdelta, color=bar_cols, alpha=0.85, height=0.6)
axD.set_yticks(np.arange(len(genes))[::-1])
axD.set_yticklabels(ylab, fontsize=6.8)
axD.axvline(0, color="#888888", lw=0.8)
axD.set_xlabel("Mean paired logCPM change (post minus pre)")
axD.set_title("Six-gene paired effects (allograft-level)",
              loc="left", fontsize=8)
logo_neg = val("logo_negative_fraction")
lodo_neg = val("lodo_negative_fraction")
pg_neg = (pg["estimate"] < 0).sum()
axD.text(0.02, 0.02,
         f"{pg_neg}/6 genes negative; LOGO negative {logo_neg:.0%}; "
         f"LODO negative {lodo_neg:.0%} ({int(round(val('lodo_negative_fraction')*41))}/41). "
         "UMOD is positive; no cell-intrinsic inference.",
         transform=axD.transAxes, fontsize=6.6, color=GREY, va="bottom")
panel_letter(axD, "D")

# --------------------------------------------------------------- Panel E ---
axE.set_xlim(0, 1); axE.set_ylim(0, 1); axE.axis("off")
axE.set_xticks([]); axE.set_yticks([])   # axis('off') still draws tick labels; clear them
axE.text(0.5, 0.7, "GSE43974 paired T2/T3", fontsize=8, fontweight="bold",
         ha="center")
axE.text(0.5, 0.45, "DEFERRED\nNo local GSE43974 expression matrix\n(archive removed at evidence freeze)",
         fontsize=6.6, ha="center", color=GREY)
axE.add_patch(patches.FancyBboxPatch((0.08, 0.28), 0.84, 0.5,
                                     boxstyle="round,pad=0.05",
                                     facecolor="#FBF3E4", edgecolor=GOLD, lw=1.0))
axE.text(0.5, 0.06, "Would require reproducible reacquisition of the public matrix\n"
                    "and author T2/T3 pairing table",
         fontsize=5.6, ha="center", color=GREY)
panel_letter(axE, "E")

fig.suptitle("Figure 6  |  Within-allograft tissue-level TAL programme change after transplant IRI",
             fontsize=11, fontweight="bold", y=0.965)
fig.savefig(os.path.join(OUT, "figure6_paired_iri.png"), dpi=300)
fig.savefig(os.path.join(OUT, "figure6_paired_iri.pdf"))
print("wrote figure6_paired_iri.png/.pdf to", OUT)
