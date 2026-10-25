#!/usr/bin/env python3
"""
Figure 5 (extension) — nephron-wide tubular identity context.

Panels (analysis_extension_plan.md WP5 / manuscript_v4 storyboard):
  A  segment-by-cohort effect forest (PT/TAL/DCT/CNT; per_segment_effect)
  B  donor-by-segment identity heatmap (cohort blocks kept separate)
  C  TAL vs PT identity donor scatter (each point = donor)
  D  TAL vs other-segment identity summary (partial r, cohort-specific + pooled)
  E  failed specificity tests (PT-adjusted, TAL-minus-PT, AKI-vs-CKD)
  F  discordant TAL-PT injury coupling (directions differ; labelled not pooled)

All effect numbers are read from frozen result tables; no new models are fit
here. Donors remain cohort-specific; no donor pooling across cohorts.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

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

AKI_PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BASE = os.path.dirname(AKI_PROJ)
AKI_RES = os.path.join(AKI_PROJ, "results")
AKI_EXT = os.path.join(AKI_PROJ, "results", "extension")
DCT_RES = os.path.join(BASE, "AKI_DCT_analysis", "results")
OUT = os.path.join(AKI_PROJ, "figures", "extension")

AKI = "#B2182B"; CTRL = "#2166AC"; GREY = "#444444"
SEG_COL = {"PT": "#1B9E77", "TAL": "#D95F02", "DCT": "#7570B3", "CNT": "#E7298A"}

def panel_letter(ax, letter):
    ax.text(-0.05, 1.02, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="right", va="bottom")

# ------------------------------------------------------------------ data ---
s3c_models = pd.read_csv(os.path.join(DCT_RES, "stage3c_model_results.csv"))
s3c_scores = pd.read_csv(os.path.join(DCT_RES, "stage3c_program_scores.csv"))
fe = pd.read_csv(os.path.join(AKI_RES, "formal_evidence.csv"))
fx = pd.read_csv(os.path.join(AKI_RES, "formal_effects.csv"))
coupling = pd.read_csv(os.path.join(AKI_RES, "tal_pt_coupling_results.csv"))

# segment effects. The frozen table stores estimate, se and t_value but leaves
# df and ci_* NA on per_segment_effect rows. The CIs are therefore recomputed
# here as estimate +/- t(0.975, df) * se using the documented residual df for
# the per-segment mixed models (GSE183276 df=20 for 22 donors, GSE210622
# df=10 for 12 donors), which is the same df used to derive the frozen p values.
seg_eff = s3c_models[s3c_models["model"] == "per_segment_effect"].copy()
seg_map = {"AKI_effect_PT": "PT", "AKI_effect_TAL": "TAL",
           "AKI_effect_DCT": "DCT", "AKI_effect_CNT": "CNT"}
seg_eff["segment"] = seg_eff["term"].map(seg_map)
DOC_DF = {"GSE183276": 20.0, "GSE210622": 10.0}
seg_eff["df_use"] = seg_eff["dataset"].map(DOC_DF)
seg_eff["ci_low"] = seg_eff["estimate"] - stats.t.ppf(0.975, seg_eff["df_use"]) * seg_eff["se"]
seg_eff["ci_high"] = seg_eff["estimate"] + stats.t.ppf(0.975, seg_eff["df_use"]) * seg_eff["se"]

# donor-level identity scores (four segments)
ident = s3c_scores[s3c_scores["program"].isin(
    ["PT_identity", "TAL_identity", "DCT_identity_transport", "CNT_identity"])].copy()
ident["segment"] = ident["program"].map(
    {"PT_identity": "PT", "TAL_identity": "TAL",
     "DCT_identity_transport": "DCT", "CNT_identity": "CNT"})
ident_wide = ident.pivot_table(index=["donor", "dataset", "diagnosis"],
                               columns="segment", values="score").reset_index()

# TAL vs PT (formal evidence donor-level; primary tier)
fe_prim = fe[fe["tier"] == "primary_ge50"].copy()
fe_prim = fe_prim[fe_prim["group"].isin(["AKI", "Reference"])].copy()

# coupling results
coup_ident = coupling[coupling["analysis"] == "TAL_identity_vs_PT_identity"]
coup_inj = coupling[coupling["analysis"] == "TAL_identity_vs_PT_injury"]
coup_other = coupling[coupling["analysis"] == "TAL_vs_OtherIdentityMean"]

# failed specificity (formal_effects): select by exact (model, comparison)
# pairs so AKI_vs_CKD is never conflated with AKI_vs_Reference.
spec = fx[(fx["tier"] == "primary_ge50")].copy()
spec_sel = [
    ("PT-injury adjusted", "adjusted_PT_injury", "AKI_vs_Reference", "GSE183276"),
    ("PT-injury adjusted", "adjusted_PT_injury", "AKI_vs_nonAKI", "GSE210622"),
    ("TAL minus PT identity", "TAL_minus_PT", "AKI_vs_Reference", "GSE183276"),
    ("TAL minus PT identity", "TAL_minus_PT", "AKI_vs_nonAKI", "GSE210622"),
    ("AKI vs CKD (unadjusted)", "unadjusted", "AKI_vs_CKD", "GSE183276"),
]
spec_rows = []
for lab, model, comp, ds in spec_sel:
    sub = spec[(spec["model"] == model) & (spec["comparison"] == comp) &
               (spec["dataset"] == ds) & (spec["term"] == "AKI")]
    if len(sub) == 0:
        continue
    r = sub.iloc[0]
    spec_rows.append({
        "label": lab, "dataset": ds, "comparison": comp,
        "estimate": r["estimate"], "ci_low": r["ci_low"], "ci_high": r["ci_high"],
        "p_value": r["p_value"], "n_aki": r["n_aki"], "n_ref": r["n_ref"],
    })
spec_df = pd.DataFrame(spec_rows)

# ------------------------------------------------------------------ figure --
fig = plt.figure(figsize=(11.7, 10.5))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.9], hspace=0.62, wspace=0.35,
                      left=0.09, right=0.97, top=0.91, bottom=0.06)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[2, 0])
axF = fig.add_subplot(gs[2, 1])

# --------------------------------------------------------------- Panel A ---
rows = []
for ds in ["GSE183276", "GSE210622"]:
    for seg in ["PT", "TAL", "DCT", "CNT"]:
        r = seg_eff[(seg_eff["dataset"] == ds) & (seg_eff["segment"] == seg)]
        rows.append((ds, seg, r["estimate"].iloc[0], r["ci_low"].iloc[0],
                     r["ci_high"].iloc[0], r["p_value"].iloc[0]))
segA = pd.DataFrame(rows, columns=["dataset", "segment", "est", "lo", "hi", "p"])
for ds in ["GSE183276", "GSE210622"]:
    sub = segA[segA["dataset"] == ds]
    y = np.arange(4)[::-1]
    if ds == "GSE210622":
        y = y + 4.5
    for yi, r in zip(y, sub.itertuples(index=False)):
        axA.errorbar(r.est, yi, xerr=[[r.est - r.lo], [r.hi - r.est]],
                     fmt="o", ms=6, color=SEG_COL[r.segment], mec="white", mew=0.7,
                     lw=1.4, capsize=3, capthick=1.2, zorder=4)
        axA.text(1.55, yi, f"{r.est:.2f}", fontsize=6, va="center", ha="left",
                 color=GREY)
axA.axvline(0, color="#888888", ls="--", lw=0.9)
axA.set_yticks(np.concatenate([np.arange(4)[::-1], np.arange(4)[::-1] + 4.5]))
axA.set_yticklabels(list(segA[segA["dataset"] == "GSE183276"]["segment"]) +
                    list(segA[segA["dataset"] == "GSE210622"]["segment"]))
axA.set_ylim(-1.0, 9.0)
axA.set_xlim(-2.6, 2.6)
axA.set_xlabel("AKI minus control, segment identity z-effect [95% CI]")
axA.axhspan(3.6, 8.6, facecolor="#F2F2F2", zorder=0)
axA.text(0.02, 4.15, "GSE210622", transform=axA.get_yaxis_transform(), fontsize=6.4, color=GREY)
axA.text(0.02, 0.25, "GSE183276", transform=axA.get_yaxis_transform(), fontsize=6.4, color=GREY)
axA.set_title("Segment identity effects (donor models: GSE183276 n=22; GSE210622 n=12)", loc="left", fontsize=8)
panel_letter(axA, "A")

# --------------------------------------------------------------- Panel B ---
hb = ident_wide.copy()
hb["diagnosis"] = hb["diagnosis"].replace({"non_AKI": "non-AKI"})
hb = hb.sort_values(["dataset", "diagnosis", "donor"])
hb["donor_f"] = [f"{dx}|{x}" for dx, x in zip(hb["diagnosis"], hb["donor"])]
hb["donor_f"] = pd.Categorical(hb["donor_f"], categories=hb["donor_f"].unique())
mat = hb[["donor_f", "PT", "TAL", "DCT", "CNT"]].set_index("donor_f").T
mat = mat.reindex(["PT", "TAL", "DCT", "CNT"])
# mat is 4 (segments) x 57 (donors); plot as-is so rows = segments and the
# 4 y-tick labels spread over the full axes height. origin='lower' so row 0
# (PT) is at the bottom and the cohort labels at y=-1.35 sit below the axes.
im = axB.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                origin="lower")
axB.set_yticks(range(4))
axB.set_yticklabels(["PT", "TAL", "DCT", "CNT"])
n183 = (hb["dataset"] == "GSE183276").sum()
axB.set_xticks(range(len(hb)))
axB.set_xticklabels(hb["donor_f"], rotation=90, fontsize=3.6)
axB.axvline(n183 - 0.5, color="black", lw=1)
axB.text(n183 / 2, -1.35, "GSE183276", ha="center", fontsize=7, color=GREY)
axB.text(n183 + (len(hb) - n183) / 2, -1.35, "GSE210622", ha="center", fontsize=7, color=GREY)
axB.set_title("Donor-by-segment identity (z-score; cohort blocks kept separate)",
              loc="left", fontsize=8)
cb = fig.colorbar(im, ax=axB, fraction=0.045, pad=0.02)
cb.set_label("identity z", fontsize=6)
cb.ax.tick_params(labelsize=5.5)
panel_letter(axB, "B")

# --------------------------------------------------------------- Panel C ---
cdf = fe_prim.copy()
cdf["group"] = cdf["group"].replace({"Reference": "control"})
for ds in ["GSE183276", "GSE210622"]:
    sub = cdf[cdf["dataset"] == ds]
    for grp, col in [("AKI", AKI), ("control", CTRL)]:
        s = sub[sub["group"] == grp]
        axC.scatter(s["PT_identity"], s["TAL_identity"], s=26, color=col,
                    edgecolor="white", linewidth=0.5, alpha=0.9,
                    label=f"{ds} {grp} (n={len(s)})")
axC.axline((0, 0), slope=1, color="#B3B3B3", ls=":", lw=0.8)
axC.set_xlabel("PT identity score (donor)")
axC.set_ylabel("TAL identity score (donor)")
axC.legend(fontsize=5.8, frameon=False, loc="lower right")
axC.set_title("TAL vs PT identity, donor-level", loc="left", fontsize=8)
# panel letter inside top-left to clear the y tick labels
axC.text(0.02, 0.97, "C", transform=axC.transAxes, fontsize=12,
         fontweight="bold", ha="left", va="top")

# --------------------------------------------------------------- Panel D ---
rowsD = []
for _, r in coup_ident.iterrows():
    rowsD.append(("TAL vs PT identity", r["dataset"], r["partial_r"],
                  r["ci_lower"], r["ci_upper"]))
for _, r in coup_other.iterrows():
    rowsD.append(("TAL vs other identity (PT+DCT+CNT)", r["dataset"], r["partial_r"],
                  r["ci_lower"], r["ci_upper"]))
ddf = pd.DataFrame(rowsD, columns=["analysis", "dataset", "r", "lo", "hi"])
for i, (lab, ds) in enumerate([("TAL vs PT identity", "GSE183276"),
                               ("TAL vs PT identity", "GSE210622"),
                               ("TAL vs PT identity", "pooled"),
                               ("TAL vs other identity (PT+DCT+CNT)", "GSE183276"),
                               ("TAL vs other identity (PT+DCT+CNT)", "GSE210622"),
                               ("TAL vs other identity (PT+DCT+CNT)", "pooled")]):
    r = ddf[(ddf["analysis"] == lab) & (ddf["dataset"] == ds)]
    if len(r) == 0:
        continue
    r = r.iloc[0]
    y = 5 - i
    col = "#7A4E00" if ds == "pooled" else GREY
    axD.errorbar(r["r"], y, xerr=[[r["r"] - r["lo"]], [r["hi"] - r["r"]]],
                 fmt="o" if ds == "pooled" else "s", ms=6 if ds == "pooled" else 5,
                 color=col, mec="white", mew=0.7, lw=1.3, capsize=3, capthick=1.1)
    axD.text(0.85, y, f"r = {r['r']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}]",
             fontsize=6, va="center", ha="left", color=GREY)
axD.axvline(0, color="#888888", ls="--", lw=0.9)
axD.set_yticks([5, 4, 3, 2, 1, 0])
axD.set_yticklabels(["TAL vs PT  GSE183276 (n=27)", "TAL vs PT  GSE210622 (n=12)", "TAL vs PT  fixed-effect",
                     "TAL vs other  GSE183276 (n=22)", "TAL vs other  GSE210622 (n=12)", "TAL vs other  fixed-effect"])
axD.set_xlim(-1.1, 1.4)
axD.set_xlabel("Partial correlation (diagnosis-adjusted)")
axD.set_title("TAL identity coordination (cohort estimates + fixed-effect summary)",
              loc="left", fontsize=8)
panel_letter(axD, "D")

# --------------------------------------------------------------- Panel E ---
rowsE = []
for _, r in spec_df.iterrows():
    rowsE.append((r["label"], r["dataset"], r["estimate"], r["ci_low"], r["ci_high"],
        r["p_value"], r["n_aki"], r["n_ref"]))
edf = pd.DataFrame(rowsE, columns=["label", "dataset", "est", "lo", "hi", "p", "n_aki", "n_ref"])
# order: PT-injury adj (183, 210), TAL-minus-PT (183, 210), AKI vs CKD (183)
order = [("PT-injury adjusted", "GSE183276"), ("PT-injury adjusted", "GSE210622"),
         ("TAL minus PT identity", "GSE183276"), ("TAL minus PT identity", "GSE210622"),
         ("AKI vs CKD (unadjusted)", "GSE183276")]
ylabs, ys = [], []
for i, (lab, ds) in enumerate(order):
    r = edf[(edf["label"] == lab) & (edf["dataset"] == ds)]
    if len(r) == 0:
        continue
    r = r.iloc[0]
    y = len(order) - 1 - i
    axE.errorbar(r["est"], y, xerr=[[r["est"] - r["lo"]], [r["hi"] - r["est"]]],
                 fmt="o", ms=5.5, color="#E0A800", mec="white", mew=0.7,
                 lw=1.3, capsize=3, capthick=1.1)
    p = r["p"]
    pstr = f"P={p:.3f}" if p >= 0.001 else f"P={p:.1e}"
    axE.text(1.9, y, pstr, fontsize=6, va="center", ha="left", color=GREY)
    short_lab = {"PT-injury adjusted": "PT-injury adj.",
                 "TAL minus PT identity": "TAL\u2212PT identity",
                 "AKI vs CKD (unadjusted)": "AKI vs CKD"}.get(lab, lab)
    ds_short = "183276" if ds == "GSE183276" else "210622"
    ylabs.append(f"{short_lab}\n{ds_short}: {int(r['n_aki'])}+{int(r['n_ref'])}"); ys.append(y)
axE.axvline(0, color="#888888", ls="--", lw=0.9)
axE.set_yticks(ys)
axE.set_yticklabels(ylabs, fontsize=6.4)
axE.set_ylim(min(ys) - 0.7, max(ys) + 0.7)
axE.set_xlim(-2.8, 2.8)
axE.set_xlabel("AKI effect [95% CI]")
axE.set_title("Failed / unsupported specificity tests (retained, not pooled)\n"
              "PT-adjusted result is not replicated; no PT independence or TAL selectivity",
              loc="left", fontsize=7.2)
panel_letter(axE, "E")

# --------------------------------------------------------------- Panel F ---
rowsF = []
for _, r in coup_inj.iterrows():
    rowsF.append((r["dataset"], r["partial_r"], r["ci_lower"], r["ci_upper"]))
fdf = pd.DataFrame(rowsF, columns=["dataset", "r", "lo", "hi"])
for i, (ds, lab) in enumerate([("GSE183276", "GSE183276"), ("GSE210622", "GSE210622")]):
    r = fdf[fdf["dataset"] == ds].iloc[0]
    y = 1 - i
    col = GREY
    axF.errorbar(r["r"], y, xerr=[[r["r"] - r["lo"]], [r["hi"] - r["r"]]],
                 fmt="s", ms=7, color=col, mec="white", mew=0.8, lw=1.5,
                 capsize=3, capthick=1.2)
    n_donor = int(coup_inj.loc[coup_inj["dataset"] == ds, "n_donors"].iloc[0])
    axF.text(1.10, y, f"r = {r['r']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}], n={n_donor}",
             fontsize=6.2, va="center", ha="left", color=GREY)
axF.axvline(0, color="#888888", ls="--", lw=0.9)
axF.set_yticks([1, 0])
axF.set_yticklabels(["GSE183276", "GSE210622"])
axF.set_ylim(-0.7, 1.7)
axF.set_xlim(-1.15, 2.0)
axF.set_xlabel("Partial correlation: TAL identity vs PT injury (| diagnosis)")
axF.set_title("Discordant TAL\u2013PT injury coupling (not pooled)",
              loc="left", fontsize=8)
axF.text(0.5, -0.18, "Directions conflict across cohorts; no TAL\u2013PT injury coupling claim",
         transform=axF.transAxes, fontsize=6, color=GREY, ha="center", va="top")
panel_letter(axF, "F")

fig.suptitle("Figure 5  |  TAL programme variation accompanies broader tubular identity remodelling",
             fontsize=11, fontweight="bold", y=0.965)
fig.savefig(os.path.join(OUT, "figure5_nephron_context.png"), dpi=300)
fig.savefig(os.path.join(OUT, "figure5_nephron_context.pdf"))
print("wrote figure5_nephron_context.png/.pdf to", OUT)
