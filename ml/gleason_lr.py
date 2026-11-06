import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import (
    accuracy_score, roc_auc_score, precision_score,
    recall_score, f1_score, confusion_matrix, roc_curve,
)

# ── Load ──────────────────────────────────────────────────────────────────────

psi_raw  = pd.read_csv("RI_events_with_snoRNA.tsv", sep="\t")
meta_raw = pd.read_csv("SraRunTable.csv")

# ── Parse Gleason → binary label ─────────────────────────────────────────────

def parse_gleason(g):
    g = str(g).strip()
    if "=" in g:
        try: return int(g.split("=")[-1].split(",")[0].strip())
        except: pass
    if "+" in g:
        try: return int(g.split("+")[0]) + int(g.split("+")[1].split(",")[0])
        except: pass
    try: return int(g)
    except: return None

meta = meta_raw[["Run", "gleason"]].dropna(subset=["gleason"]).copy()
meta["gleason_score"] = meta["gleason"].apply(parse_gleason)
meta = meta.dropna(subset=["gleason_score"]).copy()
meta["gleason_score"] = meta["gleason_score"].astype(int)
meta["label"] = (meta["gleason_score"] >= 7).astype(int)  # 0 = Low, 1 = High

# ── Prepare PSI matrix ────────────────────────────────────────────────────────

# Collapse duplicate eventIDs by averaging, then transpose to sample × event
psi = (
    psi_raw
    .groupby("eventID", as_index=False)
    .mean(numeric_only=True)
    .set_index("eventID")
    .T
    .reset_index()
    .rename(columns={"index": "Run"})
)

merged = psi.merge(meta[["Run", "label"]], on="Run", how="inner")

features = [c for c in merged.columns if c not in {"Run", "label"}]
X = merged[features].values.astype(float)
y = merged["label"].values

print(f"Dataset: {X.shape[0]} samples × {X.shape[1]} events | "
      f"Low={sum(y==0)}, High={sum(y==1)}")

# ── Repeated 20-split performance — the honest number to report ───────────────
# A single 80/20 split on 136 samples is too noisy (1 sample = 3.6% accuracy).
# We repeat 20 times with different random seeds and report mean ± std.
# Each iteration is fully leak-free: imputer, var filter, scaler all fit on
# training fold only.

print("\nRunning 20 repeated stratified splits...")

sss = StratifiedShuffleSplit(n_splits=20, test_size=0.2, random_state=0)
rep_aucs_lr, rep_accs_lr, rep_f1s_lr = [], [], []
rep_aucs_svm, rep_accs_svm, rep_f1s_svm = [], [], []

for tr_idx, va_idx in sss.split(X, y):
    X_tr_r, X_va_r = X[tr_idx], X[va_idx]
    y_tr_r, y_va_r = y[tr_idx], y[va_idx]

    # Preprocessing — fit on train only
    imp_r = SimpleImputer(strategy="median").fit(X_tr_r)
    X_tr_r = imp_r.transform(X_tr_r)
    X_va_r = imp_r.transform(X_va_r)

    vm_r = X_tr_r.var(axis=0) > 0
    X_tr_r, X_va_r = X_tr_r[:, vm_r], X_va_r[:, vm_r]

    sc_r = StandardScaler().fit(X_tr_r)
    X_tr_r = sc_r.transform(X_tr_r)
    X_va_r = sc_r.transform(X_va_r)

    # Logistic Regression
    m_r = LogisticRegressionCV(
        Cs=np.logspace(-3, 0, 20), cv=5, penalty="l2",
        solver="lbfgs", class_weight="balanced",
        max_iter=5000, random_state=42,
    ).fit(X_tr_r, y_tr_r)

    rep_aucs_lr.append(roc_auc_score(y_va_r, m_r.predict_proba(X_va_r)[:, 1]))
    rep_accs_lr.append(accuracy_score(y_va_r, m_r.predict(X_va_r)))
    rep_f1s_lr.append(f1_score(y_va_r, m_r.predict(X_va_r)))

    # Linear SVM
    svm_r = CalibratedClassifierCV(
        LinearSVC(C=0.1, class_weight="balanced", max_iter=5000, random_state=42)
    ).fit(X_tr_r, y_tr_r)

    rep_aucs_svm.append(roc_auc_score(y_va_r, svm_r.predict_proba(X_va_r)[:, 1]))
    rep_accs_svm.append(accuracy_score(y_va_r, svm_r.predict(X_va_r)))
    rep_f1s_svm.append(f1_score(y_va_r, svm_r.predict(X_va_r)))

print("\nRepeated 20-split performance (mean ± std) — USE THESE FOR PRESENTATION:")
print(f"  Logistic Regression:")
print(f"    AUC      : {np.mean(rep_aucs_lr):.3f} ± {np.std(rep_aucs_lr):.3f}")
print(f"    Accuracy : {np.mean(rep_accs_lr):.3f} ± {np.std(rep_accs_lr):.3f}")
print(f"    F1       : {np.mean(rep_f1s_lr):.3f} ± {np.std(rep_f1s_lr):.3f}")
print(f"  Linear SVM:")
print(f"    AUC      : {np.mean(rep_aucs_svm):.3f} ± {np.std(rep_aucs_svm):.3f}")
print(f"    Accuracy : {np.mean(rep_accs_svm):.3f} ± {np.std(rep_accs_svm):.3f}")
print(f"    F1       : {np.mean(rep_f1s_svm):.3f} ± {np.std(rep_f1s_svm):.3f}")

# ── Single fixed split for detailed analysis (confusion matrix, coefficients) ─

X_tr, X_va, y_tr, y_va = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# Preprocessing — fit on training only
imputer = SimpleImputer(strategy="median")
X_tr = imputer.fit_transform(X_tr)
X_va = imputer.transform(X_va)

var_mask = X_tr.var(axis=0) > 0
X_tr = X_tr[:, var_mask]
X_va = X_va[:, var_mask]
features = [f for f, keep in zip(features, var_mask) if keep]

scaler = StandardScaler()
X_tr = scaler.fit_transform(X_tr)
X_va = scaler.transform(X_va)

print(f"\nAfter filtering: {X_tr.shape[1]} features | "
      f"Train={X_tr.shape[0]}, Val={X_va.shape[0]}")

# ── Logistic Regression (fixed split) ────────────────────────────────────────

model = LogisticRegressionCV(
    Cs=np.logspace(-3, 0, 20),
    cv=5,
    penalty="l2",
    solver="lbfgs",
    class_weight="balanced",
    max_iter=5000,
    random_state=42,
)
model.fit(X_tr, y_tr)
print(f"Best regularisation C = {model.C_[0]:.4f}")

y_pred = model.predict(X_va)
y_prob = model.predict_proba(X_va)[:, 1]

print("\nLogistic Regression — single split validation:")
print(f"  Accuracy  : {accuracy_score(y_va, y_pred):.4f}")
print(f"  ROC-AUC   : {roc_auc_score(y_va, y_prob):.4f}")
print(f"  Precision : {precision_score(y_va, y_pred):.4f}")
print(f"  Recall    : {recall_score(y_va, y_pred):.4f}")
print(f"  F1        : {f1_score(y_va, y_pred):.4f}")

cm = confusion_matrix(y_va, y_pred)
print(f"\nConfusion matrix (rows=actual, cols=predicted):")
print(f"  TN={cm[0,0]}  FP={cm[0,1]}")
print(f"  FN={cm[1,0]}  TP={cm[1,1]}")

y_tr_pred = model.predict(X_tr)
y_tr_prob = model.predict_proba(X_tr)[:, 1]
print(f"\nTrain performance (for reference):")
print(f"  Accuracy  : {accuracy_score(y_tr, y_tr_pred):.4f}")
print(f"  ROC-AUC   : {roc_auc_score(y_tr, y_tr_prob):.4f}")
print(f"  F1        : {f1_score(y_tr, y_tr_pred):.4f}")
print(f"\nGeneralisation gap (train - val):")
print(f"  Δ Accuracy: {accuracy_score(y_tr,y_tr_pred)-accuracy_score(y_va,y_pred):+.4f}")
print(f"  Δ ROC-AUC : {roc_auc_score(y_tr,y_tr_prob)-roc_auc_score(y_va,y_prob):+.4f}")

# ── Linear SVM (fixed split) ──────────────────────────────────────────────────

svm = CalibratedClassifierCV(
    LinearSVC(C=0.1, class_weight="balanced", max_iter=5000, random_state=42)
)
svm.fit(X_tr, y_tr)

svm_pred = svm.predict(X_va)
svm_prob = svm.predict_proba(X_va)[:, 1]
svm_cm   = confusion_matrix(y_va, svm_pred)

print("\n── Linear SVM — single split validation:")
print(f"  Accuracy  : {accuracy_score(y_va, svm_pred):.4f}")
print(f"  ROC-AUC   : {roc_auc_score(y_va, svm_prob):.4f}")
print(f"  Precision : {precision_score(y_va, svm_pred):.4f}")
print(f"  Recall    : {recall_score(y_va, svm_pred):.4f}")
print(f"  F1        : {f1_score(y_va, svm_pred):.4f}")
print(f"\nConfusion matrix:")
print(f"  TN={svm_cm[0,0]}  FP={svm_cm[0,1]}")
print(f"  FN={svm_cm[1,0]}  TP={svm_cm[1,1]}")

# ── Top coefficients ──────────────────────────────────────────────────────────

SNORNA = {
    "ENSG00000224078", "ENSG00000255717", "ENSG00000203875",
    "ENSG00000142541", "ENSG00000234741", "ENSG00000197989",
    "ENSG00000196756", "ENSG00000163597", "ENSG00000174365",
}
GENE_SYMBOLS = {
    "ENSG00000204628": "GPC1",    "ENSG00000142534": "RPS11",
    "ENSG00000166012": "TAMM41",  "ENSG00000224078": "SNHG15",
    "ENSG00000105373": "METTL21A","ENSG00000255717": "SNHG29",
    "ENSG00000156976": "EIF4A2",  "ENSG00000142541": "RPL13A",
    "ENSG00000140988": "RPS16",   "ENSG00000008988": "RPS20",
    "ENSG00000203875": "SNHG12",  "ENSG00000163938": "GNL3",
    "ENSG00000215021": "YRDC",    "ENSG00000174444": "RPS3",
    "ENSG00000234741": "SNHG20",  "ENSG00000197989": "SNHG8",
    "ENSG00000196756": "SNHG17",  "ENSG00000163597": "SNHG16",
}

coef = model.coef_[0]
N_TOP = 15

def annotate(event_id):
    ensg = event_id.split(";")[0].split(".")[0]
    symbol = GENE_SYMBOLS.get(ensg, ensg)
    tag = " [snoRNA]" if ensg in SNORNA else ""
    return f"{symbol}{tag}"

top_pos = np.argsort(coef)[::-1][:N_TOP]
top_neg = np.argsort(coef)[:N_TOP]

print(f"\nTop {N_TOP} positive coefficients (associated with High Gleason):")
for i in top_pos:
    print(f"  {coef[i]:+.4f}  {annotate(features[i])}")

print(f"\nTop {N_TOP} negative coefficients (associated with Low Gleason):")
for i in top_neg:
    print(f"  {coef[i]:+.4f}  {annotate(features[i])}")

# ── Figure ────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(
    "Gleason Classification — LR & SVM (L2, leak-free, 20 repeated splits)",
    fontsize=13, y=1.02,
)

# Confusion matrix (from fixed split)
sns.heatmap(
    cm, annot=True, fmt="d", cmap="Blues", ax=axes[0],
    xticklabels=["Pred Low", "Pred High"],
    yticklabels=["Act Low", "Act High"],
    linewidths=0.5, cbar=False,
    annot_kws={"size": 14, "fontweight": "bold"},
)
axes[0].set_title("LR Confusion Matrix\n(single fixed split)")

# ROC curve (from fixed split)
fpr, tpr, _ = roc_curve(y_va, y_prob)
fpr_s, tpr_s, _ = roc_curve(y_va, svm_prob)
axes[1].plot(fpr, tpr, lw=2, label=f"LR  AUC = {roc_auc_score(y_va, y_prob):.3f}")
axes[1].plot(fpr_s, tpr_s, lw=2, linestyle="--", label=f"SVM AUC = {roc_auc_score(y_va, svm_prob):.3f}")
axes[1].plot([0, 1], [0, 1], "k--", lw=1)
axes[1].fill_between(fpr, tpr, alpha=0.1)
axes[1].set_xlabel("False Positive Rate")
axes[1].set_ylabel("True Positive Rate")
axes[1].set_title("ROC Curve (single fixed split)")
axes[1].legend()

# Repeated splits — AUC distribution boxplot
axes[2].boxplot(
    [rep_aucs_lr, rep_aucs_svm],
    tick_labels=["LR", "SVM"],
    patch_artist=True,
    boxprops=dict(facecolor="#aec6cf", color="black"),
    medianprops=dict(color="black", lw=2),
)
axes[2].set_ylabel("ROC-AUC")
axes[2].set_title("AUC Distribution\n(20 repeated splits)")
axes[2].set_ylim(0.5, 1.05)
axes[2].axhline(np.mean(rep_aucs_lr), color="#2196F3", lw=1.2, ls="--",
                label=f"LR mean={np.mean(rep_aucs_lr):.3f}")
axes[2].axhline(np.mean(rep_aucs_svm), color="#FF7043", lw=1.2, ls="--",
                label=f"SVM mean={np.mean(rep_aucs_svm):.3f}")
axes[2].legend(fontsize=8)
axes[2].grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("gleason_lr_results.png", dpi=150, bbox_inches="tight")
print("\nFigure saved: gleason_lr_
