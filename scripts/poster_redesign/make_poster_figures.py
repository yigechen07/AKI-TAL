#!/usr/bin/env python3
"""Journal-style poster figures (redesign 2026-11-09).

Reads only frozen result tables; writes SVG + PNG to
poster/redesign_2026-11-09/figures/. Run from the project root with the
original refresh venv (Matplotlib 3.11.2):
    deliverables/2026-10-21_refresh/.build/venv/bin/python -I \
        code/poster_redesign/make_poster_figures.py .
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
R = root / "results"
OUT = root / "poster" / "redesign_2026-11-09" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

INK, GRAY, LIGHT, RULE = "#141413", "#5E5D59", "#B9B7AE", "#E3E1D8"
CLAY, SKY = "#D97757", "#6A9BCC"
W = 13.1  # inches: inner width of one poster column card

plt.rcParams.update({
    "font.family": ["Avenir Next", "Helvetica Neue", "Arial Unicode MS"], "font.size": 19, "axes.labelsize": 20,
    "axes.titlesize": 21, "xtick.labelsize": 18, "ytick.labelsize": 18,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK, "xtick.color": GRAY,
    "ytick.color": GRAY, "text.color": INK, "axes.linewidth": 1.1,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "axes.unicode_minus": True, "svg.hashsalt": "aki-tal-poster",
    "figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none",
})


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def fmtp(p):
    if p < 0.001:
        m, e = f"{p:.1e}".split("e")
        return f"{m}×10" + str(int(e)).translate(SUP)
    return f"{p:.2g}"


def num(v, nd=2):
    return f"{v:.{nd}f}".replace("-", "−")


def panel(ax, letter, x=-0.02, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=26, fontweight="bold",
            va="bottom", ha="right" if x < 0 else "left")


fe = pd.read_csv(R / "formal_effects.csv")
be = pd.read_csv(R / "bulk_effects.csv")
ev = pd.read_csv(R / "formal_evidence.csv")
bs = pd.read_csv(R / "bulk_scores.csv")
g126 = pd.read_csv(R / "gse126805_results.csv").set_index("label")
smp = pd.read_csv(root / "metadata" / "gse126805_samples.csv")
cons = pd.read_csv(R / "extension" / "consensus_module_effects.csv")
conc = pd.read_csv(R / "extension" / "tal_pseudobulk_concordance.csv")


def fx(ds, comp, model):
    r = fe[(fe.dataset == ds) & (fe.comparison == comp) & (fe.model == model)]
    assert len(r) == 1
    return r.iloc[0]


def bx(ds, model):
    r = be[(be.dataset == ds) & (be.model == model) & (be.analysis_level == "primary")]
    assert len(r) == 1
    return r.iloc[0]


# --------------------------------------------------------------- Figure 2
strips = []
p = ev[ev.tier == "primary_ge50"]
for ds, mod, ctl, ctl_lab in [("GSE183276", "scRNA-seq", "Reference", "Reference"),
                              ("GSE210622", "snRNA-seq", "Reference", "Non-AKI")]:
    d = p[p.dataset == ds]
    strips.append((ds, mod, ctl_lab, d[d.group == ctl].TAL_identity.values,
                   d[d.group == "AKI"].TAL_identity.values))
q = bs[bs.analysis_role == "primary"]
for ds, mod, ctl, ctl_lab in [("GSE139061", "Bulk RNA-seq", "Reference", "Reference"),
                              ("GSE30718", "Microarray", "Stable", "Stable Tx")]:
    d = q[q.dataset == ds]
    strips.append((ds, mod, ctl_lab, d[d.group == ctl].TAL_score.values,
                   d[d.group == "AKI"].TAL_score.values))

rows = []
for ds, comp, lab, test in [("GSE183276", "AKI_vs_Reference", "scRNA-seq", "OLS t"),
                            ("GSE210622", "AKI_vs_nonAKI", "snRNA-seq", "OLS t")]:
    r = fx(ds, comp, "unadjusted")
    rows.append((ds, lab, f"{int(r.n_aki)} vs {int(r.n_ref)} donors", r.estimate, r.ci_low, r.ci_high, r.p_value, test))
for ds, lab, unit in [("GSE139061", "Bulk RNA-seq", "samples"), ("GSE30718", "Microarray", "donors")]:
    r = bx(ds, "unadjusted")
    rows.append((ds, lab, f"{int(r.n_aki)} vs {int(r.n_ref)} {unit}", r.estimate, r.ci_low, r.ci_high, r.p_value, "Welch t"))
r = g126.loc["baseline_post_vs_pre_delta"]
rows.append(("GSE126805", "Paired allografts", "41 kidneys, post − pre", r.estimate, r.ci_low, r.ci_high, r.p_value, "paired t"))

fig = plt.figure(figsize=(W, 8.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.72], hspace=0.6, wspace=0.42)
rng = np.random.default_rng(7)
for i, (ds, mod, ctl_lab, c, a) in enumerate(strips):
    ax = fig.add_subplot(gs[0, i])
    for x, vals, col in [(0, c, SKY), (1, a, CLAY)]:
        jit = rng.uniform(-0.17, 0.17, len(vals))
        ax.scatter(x + jit, vals, s=46, color=col, alpha=0.85, lw=0, zorder=3)
        ax.plot([x - 0.3, x + 0.3], [vals.mean()] * 2, color=INK, lw=2.6, zorder=4)
    ax.set_xticks([0, 1], [f"{ctl_lab}\n{len(c)}", f"AKI\n{len(a)}"])
    ax.set_xlim(-0.6, 1.6)
    ax.tick_params(axis="x", length=0, pad=6)
    ax.set_title(f"{ds}\n", fontweight="bold", fontsize=19, pad=4)
    ax.text(0.5, 1.02, mod, transform=ax.transAxes, ha="center", va="bottom", color=GRAY, fontsize=17)
    if i == 0:
        ax.set_ylabel("TAL score")
        panel(ax, "a", x=-0.42, y=1.13)

ax = fig.add_subplot(gs[1, :])
pos = ax.get_position()
# plot region = 30%..66% of the row; text columns sit outside it in axes fractions
ax.set_position([pos.x0 + 0.335 * pos.width, pos.y0, 0.29 * pos.width, pos.height])
fx_ = lambda f: (f - 0.335) / 0.29  # row fraction -> axes fraction
tr = ax.get_yaxis_transform()
n = len(rows)
for k, (ds, lab, nlab, est, lo, hi, pv, test) in enumerate(rows):
    y = n - 1 - k
    ax.plot([lo, hi], [y, y], color=INK, lw=2.4, solid_capstyle="round")
    ax.scatter([est], [y], s=190, color=CLAY, edgecolor=INK, lw=1.2, zorder=4, marker="D" if ds == "GSE126805" else "o")
    ax.text(fx_(0.0), y, ds, transform=tr, ha="left", va="center", fontweight="bold", fontsize=18)
    ax.text(fx_(0.162), y, lab, transform=tr, ha="left", va="center", color=GRAY, fontsize=16)
    ax.text(fx_(0.655), y, f"{num(est)} ({num(lo)} to {num(hi)})", transform=tr, ha="left", va="center", fontsize=16)
    ax.text(fx_(0.895), y, f"{fmtp(pv)}", transform=tr, ha="left", va="center", color=GRAY, fontsize=16)
ax.axvline(0, color=GRAY, lw=1.2, ls=(0, (4, 3)))
ax.set_xlim(-1.75, 0.15)
ax.set_xticks([-1.5, -1.0, -0.5, 0])
ax.spines["left"].set_visible(False)
ax.set_yticks([])
ax.set_ylim(-0.6, n - 0.3)
ax.set_xlabel("Difference in TAL score,\nAKI − control (95% CI)")
ax.text(fx_(0.655), n - 0.45, "Estimate (95% CI)", transform=tr, ha="left", va="bottom", color=GRAY, fontsize=15, fontstyle="italic")
ax.text(fx_(0.895), n - 0.45, "P", transform=tr, ha="left", va="bottom", color=GRAY, fontsize=15, fontstyle="italic")
ax.text(fx_(-0.03), n - 0.3, "b", transform=tr, fontsize=26, fontweight="bold", va="bottom", ha="left")
save(fig, "fig2_scores_forest")

# --------------------------------------------------------------- Figure 3
pre = smp[smp.timepoint == "baseline_pre"].set_index("kidney").TAL_score
post = smp[smp.timepoint == "baseline_post"].set_index("kidney").TAL_score
k = pre.index.intersection(post.index)
pre, post = pre[k], post[k]
dlt = (post - pre).sort_values()
assert len(dlt) == 41
fig = plt.figure(figsize=(W, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.35], wspace=0.32)
ax = fig.add_subplot(gs[0])
for kid in k:
    col = CLAY if post[kid] < pre[kid] else LIGHT
    ax.plot([0, 1], [pre[kid], post[kid]], color=col, lw=1.6, alpha=0.8, zorder=2 if col == CLAY else 1)
    ax.scatter([0, 1], [pre[kid], post[kid]], s=22, color=col, lw=0, zorder=3)
ax.plot([0, 1], [pre.mean(), post.mean()], color=INK, lw=4.0, zorder=5, solid_capstyle="round")
ax.set_xticks([0, 1], ["Before\nreperfusion", "After\nreperfusion"])
ax.tick_params(axis="x", length=0, pad=6)
ax.set_xlim(-0.25, 1.25)
ax.set_ylabel("TAL score")
panel(ax, "a", x=-0.3, y=1.02)
ax = fig.add_subplot(gs[1])
cols = [CLAY if v < 0 else LIGHT for v in dlt.values]
ax.bar(np.arange(len(dlt)), dlt.values, color=cols, width=0.78)
ax.axhline(0, color=GRAY, lw=1.1)
r = g126.loc["baseline_post_vs_pre_delta"]
ax.axhline(r.estimate, color=INK, lw=2.2, ls=(0, (5, 3)))
ax.text(len(dlt) - 0.5, r.estimate - 0.08, f"mean Δ = {num(r.estimate)}", ha="right", va="top", fontsize=17)
ax.text(0, dlt.max() * 0.95, f"{int((dlt < 0).sum())} of {len(dlt)} kidneys lower", ha="left", va="top",
        fontsize=18, fontweight="bold", color=CLAY)
ax.set_xticks([])
ax.spines["bottom"].set_visible(False)
ax.set_xlabel("Kidney allografts, ranked by change")
ax.set_ylabel("Δ TAL score (after − before)")
panel(ax, "b", x=-0.16, y=1.02)
save(fig, "fig3_paired_allografts")

# --------------------------------------------------------------- Figure 4
fig = plt.figure(figsize=(W, 5.2))
gs = fig.add_gridspec(1, 2, width_ratios=[0.75, 1.6], wspace=0.22)
ax = fig.add_subplot(gs[0])
lim = 6
x, y = conc.logFC_GSE183276.clip(-lim, lim), conc.logFC_GSE210622.clip(-lim, lim)
rep = conc.support_tier == "FDR<0.05 both"
ax.scatter(x[~rep], y[~rep], s=5, color=LIGHT, alpha=0.5, lw=0, rasterized=True)
ax.scatter(x[rep], y[rep], s=9, color=CLAY, alpha=0.75, lw=0, rasterized=True)
ax.axhline(0, color=GRAY, lw=0.9); ax.axvline(0, color=GRAY, lw=0.9)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xticks([-5, 0, 5]); ax.set_yticks([-5, 0, 5])
ax.set_aspect("equal")
ax.set_xlabel("log$_2$FC, GSE183276")
ax.set_ylabel("log$_2$FC, GSE210622")
rho = conc[["logFC_GSE183276", "logFC_GSE210622"]].corr("spearman").iloc[0, 1]
ax.set_title(f"{int(rep.sum()):,} replicated genes (orange)\nSpearman ρ = {rho:.2f}",
             fontsize=16, color=GRAY, loc="left", pad=10)
fig.text(0.0, 0.97, "a", fontsize=26, fontweight="bold", va="top")
ax = fig.add_subplot(gs[1])
crow = []
for ds, lab, nlab in [("GSE139061", "Bulk RNA-seq", "39 AKI vs 9 Ref."),
                      ("GSE30718", "Microarray", "26 AKI vs 11 stable"),
                      ("GSE126805", "Paired allografts", "41 kidneys, post − pre")]:
    rr = cons[(cons.cohort == ds) & (cons.model == "main_signed_z")].iloc[0]
    crow.append((ds, lab, nlab, rr.estimate, rr.ci_low, rr.ci_high, rr.p_holm))
for kk, (ds, lab, nlab, est, lo, hi, ph) in enumerate(crow):
    yy = 2 - kk
    ax.plot([lo, hi], [yy, yy], color=INK, lw=2.4, solid_capstyle="round")
    ax.scatter([est], [yy], s=190, color=CLAY, edgecolor=INK, lw=1.2, zorder=4, marker="D" if ds == "GSE126805" else "o")
    ax.text(-1.75, yy + 0.12, ds, ha="left", va="center", fontweight="bold", fontsize=18)
    ax.text(-1.75, yy - 0.22, lab, ha="left", va="center", color=GRAY, fontsize=15)
    ax.text(1.22, yy + 0.12, f"{num(est)} ({num(lo)} to {num(hi)})", ha="left", va="center", fontsize=16)
    ax.text(1.22, yy - 0.22, f"Holm P = {ph:.2g}", ha="left", va="center", color=GRAY, fontsize=15)
ax.axvline(0, color=GRAY, lw=1.2, ls=(0, (4, 3)))
ax.set_xlim(-1.8, 1.95); ax.set_ylim(-0.55, 2.45)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([0, 0.5, 1.0]); ax.spines["bottom"].set_bounds(-0.15, 1.15)
ax.set_xlabel("Response score, AKI − control (95% CI)")
fig.text(0.42, 0.97, "b", fontsize=26, fontweight="bold", va="top")
save(fig, "fig4_gene_extension")

# --------------------------------------------------------------- Figure 5
spec = []
for ds, comp in [("GSE183276", "AKI_vs_Reference"), ("GSE210622", "AKI_vs_nonAKI")]:
    u, a = fx(ds, comp, "unadjusted"), fx(ds, comp, "adjusted_PT_injury")
    spec.append((ds, (u.estimate, u.ci_low, u.ci_high), (a.estimate, a.ci_low, a.ci_high)))
for ds in ["GSE139061", "GSE30718"]:
    u, a = bx(ds, "unadjusted"), bx(ds, "adjusted_PT_injury")
    spec.append((ds, (u.estimate, u.ci_low, u.ci_high), (a.estimate, a.ci_low, a.ci_high)))
u, a = g126.loc["baseline_post_vs_pre_FE"], g126.loc["FE_baseline_post_PT_inj_adj"]
spec.append(("GSE126805", (u.estimate, u.ci_low, u.ci_high), (a.estimate, a.ci_low, a.ci_high)))
extra = []
for ds, comp, model, lab in [("GSE183276", "AKI_vs_Reference", "TAL_minus_PT", "TAL − PT, GSE183276"),
                             ("GSE210622", "AKI_vs_nonAKI", "TAL_minus_PT", "TAL − PT, GSE210622"),
                             ("GSE183276", "AKI_vs_CKD", "unadjusted", "AKI vs CKD, GSE183276")]:
    rr = fx(ds, comp, model)
    extra.append((lab, rr.estimate, rr.ci_low, rr.ci_high))

fig, ax = plt.subplots(figsize=(W, 5.4))
ytop = len(spec) + len(extra) + 0.6
yy = ytop
ticks = []
for ds, un, ad in spec:
    yy -= 1
    for off, (e, lo, hi), face, lab in [(0.17, un, CLAY, "Unadjusted"), (-0.17, ad, "white", "Adjusted for PT injury")]:
        ax.plot([lo, hi], [yy + off] * 2, color=INK if face == CLAY else GRAY, lw=2.2, solid_capstyle="round")
        ax.scatter([e], [yy + off], s=130, color=face, edgecolor=INK, lw=1.4, zorder=4, label=lab)
    ticks.append((yy, ds))
yy -= 0.55
ax.plot([-3.35, 1.1], [yy, yy], color=RULE, lw=1.4, clip_on=False)
for lab, e, lo, hi in extra:
    yy -= 1
    ax.plot([lo, hi], [yy, yy], color=GRAY, lw=2.2, solid_capstyle="round")
    ax.scatter([e], [yy], s=130, color=LIGHT, edgecolor=INK, lw=1.2, zorder=4)
    ticks.append((yy, lab))
for t, lab in ticks:
    ax.text(-3.35, t, lab, ha="left", va="center", fontsize=17,
            fontweight="bold" if lab.startswith("GSE") else "normal")
ax.axvline(0, color=GRAY, lw=1.2, ls=(0, (4, 3)))
ax.set_xlim(-3.4, 1.15); ax.set_ylim(yy - 0.7, ytop - 0.2)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([-1.5, -1.0, -0.5, 0, 0.5, 1.0]); ax.spines["bottom"].set_bounds(-2.15, 1.05)
ax.set_xlabel("Difference in TAL score (95% CI)")
h, l = ax.get_legend_handles_labels()
ax.legend(h[:2], l[:2], loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False,
          fontsize=16, handletextpad=0.3, columnspacing=1.5, borderaxespad=0.2)
save(fig, "fig5_specificity")
print("written:", sorted(p.name for p in OUT.iterdir()))
