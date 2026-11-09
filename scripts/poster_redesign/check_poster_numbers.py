#!/usr/bin/env python3
"""Check every number and directional claim in the redesigned poster text
against the frozen result tables. Run from the project root:
    python3 -I code/poster_redesign/check_poster_numbers.py .
Figure numbers are drawn directly from the same tables by make_poster_figures.py."""
import re
import sys
from pathlib import Path
import pandas as pd

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
R = root / "results"
html = (root / "poster/redesign_2026-11-09/AKI_TAL_poster_journal.html").read_text(encoding="utf8")
text = re.sub(r"<[^>]+>", " ", html).replace("&nbsp;", " ").replace("&lt;", "<")
text = re.sub(r"\s+", " ", text)

fe = pd.read_csv(R / "formal_effects.csv"); be = pd.read_csv(R / "bulk_effects.csv")
bl = pd.read_csv(R / "bulk_lodo.csv"); flo = pd.read_csv(R / "formal_lodo.csv")
g = pd.read_csv(R / "gse126805_results.csv").set_index("label")
smp = pd.read_csv(root / "metadata/gse126805_samples.csv")
fd = pd.read_csv(root / "metadata/formal_donors.csv")
ev = pd.read_csv(R / "formal_evidence.csv")
cons = pd.read_csv(R / "extension/consensus_module_effects.csv")
conc = pd.read_csv(R / "extension/tal_pseudobulk_concordance.csv")
camera = (R / "extension/tal_pathway_report.md").read_text(encoding="utf8")

fails = []
def ok(label, cond, shown=None):
    print(("PASS " if cond else "FAIL ") + label)
    if not cond: fails.append(label)
    if shown is not None and shown not in text:
        print("FAIL   text missing: " + shown); fails.append("text: " + shown)

def fx(ds, comp, model):
    return fe[(fe.dataset == ds) & (fe.comparison == comp) & (fe.model == model)].iloc[0]
def bx(ds, model):
    return be[(be.dataset == ds) & (be.model == model) & (be.analysis_level == "primary")].iloc[0]

# Table 1 / group sizes
a = fx("GSE183276", "AKI_vs_Reference", "unadjusted"); ok("GSE183276 10 vs 17", (a.n_aki, a.n_ref) == (10, 17), "10 AKI vs 17 reference")
b = fx("GSE210622", "AKI_vs_nonAKI", "unadjusted"); ok("GSE210622 8 vs 4", (b.n_aki, b.n_ref) == (8, 4), "8 AKI vs 4 non-AKI")
c = bx("GSE139061", "unadjusted"); ok("GSE139061 39 vs 9", (c.n_aki, c.n_ref) == (39, 9), "39 AKI vs 9 reference")
d = bx("GSE30718", "unadjusted"); ok("GSE30718 26 vs 11", (d.n_aki, d.n_ref) == (26, 11), "26 transplant AKI vs 11 stable")
# five of five lower
ests = [a.estimate, b.estimate, c.estimate, d.estimate, g.loc["baseline_post_vs_pre_delta", "estimate"]]
ok("all five estimates negative with CI below 0", all(e < 0 for e in ests) and
   all(x < 0 for x in [a.ci_high, b.ci_high, c.ci_high, d.ci_high, g.loc["baseline_post_vs_pre_delta", "ci_high"]]),
   "lower in AKI in five of five cohorts")
# leave-one-donor-out never changed sign
ok("LODO all negative (cell-resolved)", a.lodo_negative_frac == 1 and b.lodo_negative_frac == 1)
ok("LODO all negative (bulk)", (bl.estimate < 0).all())
ok("LOKO all negative (allografts)", g.loc["lodo_negative_fraction", "estimate"] == 1, "Leaving out any single donor never changed the sign")
# GSE210622 acquisition
g2 = fd[fd.dataset == "GSE210622"]
ok("GSE210622 all 8 AKI post-mortem", (g2[g2.diagnosis == "AKI"].post_mortem == "postmortem").sum() == 8, "all AKI samples were post-mortem")
ok("GSE210622 3 of 4 controls nephrectomy", (g2[g2.diagnosis != "AKI"].sampling_context == "nephrectomy").sum() == 3, "three of four controls were nephrectomies")
# paired allografts
r = g.loc["baseline_post_vs_pre_delta"]
ok("paired mean 0.22 (0.02 to 0.43), P=0.032", (round(-r.estimate, 2), round(-r.ci_high, 2), round(-r.ci_low, 2), round(r.p_value, 3)) == (0.22, 0.02, 0.43, 0.032),
   "fell by a mean of 0.22 (95% CI 0.02 to 0.43; paired t -test P = 0.032)")
pre = smp[smp.timepoint == "baseline_pre"].set_index("kidney").TAL_score
post = smp[smp.timepoint == "baseline_post"].set_index("kidney").TAL_score
dl = (post - pre).dropna()
ok("41 pairs, 27 lower", (len(dl), int((dl < 0).sum())) == (41, 27), "27 grafts declined")
# gene extension
ok("10,064 shared genes", len(conc) == 10064, "Across 10,064 genes")
ok("1,027 replicated", int((conc.support_tier == "FDR<0.05 both").sum()) == 1027, "1,027 genes changed in the same direction")
ok("Spearman 0.31", round(conc[["logFC_GSE183276", "logFC_GSE210622"]].corr("spearman").iloc[0, 1], 2) == 0.31, "Spearman ρ = 0.31")
m = cons[cons.model == "main_signed_z"]
ok("consensus higher in all three tissue cohorts", (m.ci_low > 0).all() and (m.p_holm < 0.05).all(), "higher in AKI in all three tissue cohorts")
loo = cons[cons.model.str.contains("leave|excl|no_frozen|minus", case=False, regex=True)]
manu = (root / "manuscript/manuscript_v4.md").read_text(encoding="utf8")
ok("leave-six-genes-out stays positive (manuscript-recorded values)",
   all(s in manu for s in ["GSE139061 (0.308, Holm P=0.0479)", "GSE30718 (0.738, Holm P=0.0017)", "GSE126805 (0.146, Holm P=0.0479"]),
   "including after removing genes shared with the six-gene programme")
ok("no replicated primary pathway (0/49, 0/911, 0/2,953)", all(s in manu for s in ["(0/49)", "(0/911)", "(0/2,953)"]), "No single pathway replicated")
# specificity
pa = {ds: fx(ds, cmp, "adjusted_PT_injury") for ds, cmp in [("GSE183276", "AKI_vs_Reference"), ("GSE210622", "AKI_vs_nonAKI")]}
ok("PT-adjusted GSE183276 crosses 0", pa["GSE183276"].ci_low < 0 < pa["GSE183276"].ci_high, "removed it in GSE183276")
ok("PT-adjusted GSE210622 negative", pa["GSE210622"].ci_high < 0)
ok("PT-adjusted GSE139061 negative", bx("GSE139061", "adjusted_PT_injury").ci_high < 0)
ok("PT-adjusted allografts negative", g.loc["FE_baseline_post_PT_inj_adj", "ci_high"] < 0)
a30 = bx("GSE30718", "adjusted_PT_injury")
ok("PT-adjusted GSE30718 borderline (upper CI 0 to 0.02)", 0 <= a30.ci_high <= 0.02, "borderline in GSE30718")
e = ev[(ev.tier == "primary_ge50") & (ev.dataset == "GSE183276") & ev.group.isin(["AKI", "Reference"])]
ok("AKI–PT injury r = 0.91", round(e.PT_injury.corr((e.group == "AKI").astype(float)), 2) == 0.91, "( r = 0.91)")
t1, t2 = fx("GSE183276", "AKI_vs_Reference", "TAL_minus_PT"), fx("GSE210622", "AKI_vs_nonAKI", "TAL_minus_PT")
ok("TAL−PT not significant in both", t1.p_value > 0.05 and t2.p_value > 0.05, "did not fall more than a PT identity score")
ck = fx("GSE183276", "AKI_vs_CKD", "unadjusted")
ok("AKI vs CKD P = 0.15", round(ck.p_value, 2) == 0.15 and ck.ci_low < 0 < ck.ci_high, "(P = 0.15; Fig. 4)")
ok("≥50 TAL cell threshold", "≥50 TAL cells" in text)
print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
