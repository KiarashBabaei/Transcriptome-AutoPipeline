"""
ml/gleason_lr.py — Machine-learning analysis of snoRNA-containing retained-intron
(RI) PSI values, with patient-aware (leakage-free) evaluation.

Input : results/RI_events_with_snoRNA.tsv  (SUPPA2 PSI, events x SRR runs)
        metadata/SraRunTable.csv           (SRA run table, PRJNA751490)
Output: results/ml/*.tsv and figures/ml_results.png

Usage : python3 ml/gleason_lr.py [psi.tsv] [SraRunTable.csv] [--n-perm 1000]

Design decisions (fixed BEFORE looking at results):
  * Bulk libraries only (the 12 single-cell runs are a different assay).
  * Each library (BioSample) was sequenced on 4 lanes = 4 SRR runs. Lanes are
    technical replicates, so PSI is averaged per library (n = 34, not 136).
  * Tumour and adjacent-normal libraries of the same patient carry the same
    Gleason score. Patient ID is not in the run table, so a proxy
    (age + Gleason string) is used and every evaluation keeps all libraries of a
    patient together (leave-one-patient-out, LOPO).
  * One pre-specified model for every task: L2 logistic regression, C = 0.1,
    class_weight="balanced". All preprocessing (missingness filter, median
    imputation, zero-variance filter, scaling) is fitted inside each training fold.
  * Metrics: ROC-AUC with 95% patient-bootstrap CI, balanced accuracy, MCC, and a
    permutation p-value (labels shuffled at patient level, whole LOPO re-run).
  * All analyses below are reported, whatever the result.

Analyses
  M0  Naive evaluation (runs as samples, random 80/20 x 20) — reference only,
      shows how much replicate leakage inflates performance.
  M1  Gleason Low (<= 6) vs High (>= 7) from event-level PSI           [LOPO]
  M2  Same task from snoRNA-level PSI (events averaged per snoRNA)     [LOPO]
  M3  Tumour vs adjacent normal from event-level PSI                   [LOPO]
  M4  Unsupervised: PCA + hierarchical clustering; agreement of clusters with
      tumour/normal labels (adjusted Rand index, permutation p).

Requirements: pip install pandas numpy scikit-learn scipy matplotlib joblib
"""
import argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.model_selection import LeaveOneGroupOut, StratifiedShuffleSplit
from sklearn.metrics import (roc_auc_score, roc_curve, balanced_accuracy_score,
                             matthews_corrcoef, confusion_matrix, adjusted_rand_score)

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("psi", nargs="?", default="results/RI_events_with_snoRNA.tsv")
ap.add_argument("meta", nargs="?", default="metadata/SraRunTable.csv")
ap.add_argument("--n-perm", type=int, default=1000)
ap.add_argument("--n-boot", type=int, default=2000)
args = ap.parse_args()
OUT, FIG = Path("results/ml"), Path("figures")
OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
SEED = 42


# ── model ────────────────────────────────────────────────────────────────────
class MissingFilter(BaseEstimator, TransformerMixin):
    """Drop features missing in > max_frac of the TRAINING samples."""
    def __init__(self, max_frac=0.5): self.max_frac = max_frac
    def fit(self, X, y=None): self.keep_ = np.isnan(X).mean(0) <= self.max_frac; return self
    def transform(self, X): return X[:, self.keep_]


def make_model():
    return Pipeline([("miss", MissingFilter(0.5)), ("imp", SimpleImputer(strategy="median")),
                     ("var", VarianceThreshold(0.0)), ("sc", StandardScaler()),
                     ("lr", LogisticRegression(C=0.1, class_weight="balanced", max_iter=5000))])


def gleason_total(s):
    s = str(s)
    return int(s.split("=")[1].split(",")[0]) if "=" in s else sum(map(int, s.split("+")))


# ── data ─────────────────────────────────────────────────────────────────────
psi_raw = pd.read_csv(args.psi, sep="\t")
sno_map = psi_raw[["eventID", "snoRNAName"]].drop_duplicates()
psi = psi_raw.drop_duplicates("eventID").set_index("eventID").drop(columns="snoRNAName")

meta = pd.read_csv(args.meta)
meta = meta[(meta.LibrarySource == "TRANSCRIPTOMIC") & meta.Run.isin(psi.columns)].copy()
meta["patient"] = np.where(meta.gleason.notna(), meta.age.astype(str) + "|" + meta.gleason.astype(str), meta.BioSample)

lib = meta.groupby("BioSample").agg(condition=("condition", "first"), gleason=("gleason", "first"),
                                    patient=("patient", "first"), runs=("Run", list))
X_ev = pd.DataFrame({b: np.nanmean(psi[r].to_numpy(), axis=1) for b, r in lib.runs.items()},
                    index=psi.index).T                                    # libraries x events
X_sno = (X_ev.T.join(sno_map.set_index("eventID")).groupby("snoRNAName").mean().T)  # libraries x snoRNAs

print(f"Bulk runs: {len(meta)} | libraries: {len(lib)} | events: {X_ev.shape[1]} | snoRNAs: {X_sno.shape[1]}")


# ── evaluation helpers ───────────────────────────────────────────────────────
def lopo_proba(X, y, g):
    p = np.full(len(y), np.nan)
    for tr, te in LeaveOneGroupOut().split(X, y, g):
        if len(np.unique(y[tr])) == 2:
            p[te] = make_model().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return p


def shuffle_by_patient(y, g, rng):
    pats = np.unique(g)
    lab = np.array([y[g == p][0] for p in pats])
    if all((y[g == p] == y[g == p][0]).all() for p in pats):          # patient-level label
        return np.array([dict(zip(pats, rng.permutation(lab)))[x] for x in g])
    yp = y.copy()                                                      # within-patient label (T vs N)
    for p in pats:
        idx = np.where(g == p)[0]
        if rng.random() < 0.5: yp[idx] = 1 - yp[idx]
    return yp


def one_permutation(X, y, g, seed):
    yp = shuffle_by_patient(y, g, np.random.default_rng(seed))
    return roc_auc_score(yp, lopo_proba(X, yp, g))


def evaluate(name, X, y, g):
    p = lopo_proba(X, y, g); yhat = (p >= 0.5).astype(int)
    auc = roc_auc_score(y, p)
    seeds = np.random.SeedSequence(SEED).spawn(args.n_perm)
    null = np.array(Parallel(n_jobs=-1)(delayed(one_permutation)(X, y, g, s) for s in seeds))
    rng = np.random.default_rng(SEED); pats = np.unique(g); boots = []
    for _ in range(args.n_boot):
        idx = np.concatenate([np.where(g == q)[0] for q in rng.choice(pats, len(pats))])
        if len(np.unique(y[idx])) == 2: boots.append(roc_auc_score(y[idx], p[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    cm = confusion_matrix(y, yhat)
    return dict(analysis=name, n_libraries=len(y), n_patients=len(pats), n_class0=int((y == 0).sum()),
                n_class1=int((y == 1).sum()), AUC=auc, AUC_CI_low=lo, AUC_CI_high=hi,
                balanced_accuracy=balanced_accuracy_score(y, yhat), MCC=matthews_corrcoef(y, yhat),
                accuracy=(yhat == y).mean(), majority_baseline_accuracy=max(y.mean(), 1 - y.mean()),
                permutation_p=(1 + (null >= auc).sum()) / (args.n_perm + 1),
                TN=cm[0, 0], FP=cm[0, 1], FN=cm[1, 0], TP=cm[1, 1]), p, null


# ── M0: naive run-level evaluation (reference only) ──────────────────────────
gl = meta[meta.gleason.notna()]
X_run = psi[gl.Run].T.to_numpy(); y_run = (gl.gleason.map(gleason_total) >= 7).astype(int).to_numpy()
m0 = [roc_auc_score(y_run[te], make_model().fit(X_run[tr], y_run[tr]).predict_proba(X_run[te])[:, 1])
      for tr, te in StratifiedShuffleSplit(20, test_size=0.2, random_state=SEED).split(X_run, y_run)]
print(f"M0 naive (runs as samples, random 80/20 x20): AUC {np.mean(m0):.3f} ± {np.std(m0):.3f}  <- inflated")

# ── M1–M3 ───────────────────────────────────────────────────────────────────
gl_lib = lib[lib.gleason.notna()]
y_gl = (gl_lib.gleason.map(gleason_total) >= 7).astype(int).to_numpy(); g_gl = gl_lib.patient.to_numpy()
tn_lib = lib[lib.condition.isin(["Tumor", "Normal"])]
y_tn = (tn_lib.condition == "Tumor").astype(int).to_numpy(); g_tn = tn_lib.patient.to_numpy()

results, curves, nulls = [], {}, {}
for name, X, y, g in [
        ("M1 Gleason <=6 vs >=7 | event PSI", X_ev.loc[gl_lib.index].to_numpy(), y_gl, g_gl),
        ("M2 Gleason <=6 vs >=7 | snoRNA PSI", X_sno.loc[gl_lib.index].to_numpy(), y_gl, g_gl),
        ("M3 Tumour vs Normal | event PSI", X_ev.loc[tn_lib.index].to_numpy(), y_tn, g_tn)]:
    r, p, null = evaluate(name, X, y, g)
    results.append(r); curves[name] = (y, p); nulls[name] = (null, r["AUC"])
    print(f"{name}: AUC {r['AUC']:.3f} [{r['AUC_CI_low']:.2f}, {r['AUC_CI_high']:.2f}]  "
          f"bal.acc {r['balanced_accuracy']:.3f}  MCC {r['MCC']:.3f}  perm p {r['permutation_p']:.3f}")

# ── M4: unsupervised ─────────────────────────────────────────────────────────
Xu = make_model()[:-1].fit_transform(X_ev.loc[tn_lib.index].to_numpy())
pcs = PCA(n_components=2, random_state=SEED).fit(Xu); Z = pcs.transform(Xu)
clusters = AgglomerativeClustering(n_clusters=2, linkage="ward").fit_predict(Xu)
ari = adjusted_rand_score(y_tn, clusters)
rng = np.random.default_rng(SEED)
ari_null = np.array([adjusted_rand_score(shuffle_by_patient(y_tn, g_tn, rng), clusters) for _ in range(args.n_perm)])
ari_p = (1 + (ari_null >= ari).sum()) / (args.n_perm + 1)
results.append(dict(analysis="M4 Unsupervised: 2 clusters vs Tumour/Normal (ARI)", n_libraries=len(y_tn),
                    n_patients=len(np.unique(g_tn)), ARI=ari, permutation_p=ari_p,
                    PC1_var=pcs.explained_variance_ratio_[0], PC2_var=pcs.explained_variance_ratio_[1]))
print(f"M4 unsupervised: ARI {ari:.3f}  perm p {ari_p:.3f}")

res = pd.DataFrame(results)
res.insert(1, "M0_naive_AUC_reference", f"{np.mean(m0):.3f} ± {np.std(m0):.3f}")
res.to_csv(OUT / "ML_summary.tsv", sep="\t", index=False)

# ── exploratory coefficients (model refit on all M1 data) ────────────────────
full = make_model().fit(X_ev.loc[gl_lib.index].to_numpy(), y_gl)
kept = X_ev.columns[full["miss"].keep_][full["var"].get_support()]
coef = pd.DataFrame({"eventID": kept, "coef_High_Gleason": full["lr"].coef_[0]}) \
    .merge(sno_map.groupby("eventID").snoRNAName.agg(";".join).reset_index(), on="eventID") \
    .sort_values("coef_High_Gleason", key=abs, ascending=False)
coef.to_csv(OUT / "M1_coefficients_exploratory.tsv", sep="\t", index=False)

# ── figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
for name, (y, p) in curves.items():
    f, t, _ = roc_curve(y, p); ax[0].plot(f, t, lw=2, label=f"{name.split('|')[0].strip()}  AUC={roc_auc_score(y, p):.2f}")
ax[0].plot([0, 1], [0, 1], "k--", lw=1, label="chance"); ax[0].legend(fontsize=8)
ax[0].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC, leave-one-patient-out (pooled)")
for name, (null, obs) in nulls.items():
    ax[1].hist(null, bins=30, alpha=0.45, label=name.split("|")[0].strip())
    ax[1].axvline(obs, lw=2, color=ax[1].patches[-1].get_facecolor()[:3])
ax[1].set(xlabel="AUC under shuffled labels", ylabel="count", title="Permutation null (lines = observed AUC)")
ax[1].legend(fontsize=8)
for lab, col, nm in [(0, "#2E86C1", "Normal"), (1, "#C0392B", "Tumour")]:
    ax[2].scatter(Z[y_tn == lab, 0], Z[y_tn == lab, 1], c=col, label=nm, s=40, alpha=0.8)
ax[2].set(xlabel=f"PC1 ({pcs.explained_variance_ratio_[0]:.0%})", ylabel=f"PC2 ({pcs.explained_variance_ratio_[1]:.0%})",
          title=f"PCA of RI-snoRNA PSI (ARI vs T/N = {ari:.2f}, p = {ari_p:.2f})")
ax[2].legend()
plt.tight_layout(); plt.savefig(FIG / "ml_results.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: {OUT}/ML_summary.tsv, {OUT}/M1_coefficients_exploratory.tsv, {FIG}/ml_results.png")
