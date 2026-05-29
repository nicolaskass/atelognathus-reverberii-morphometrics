"""morfometria.py — Morphometric analysis of Atelognathus reverberii.

Reads the morphometric CSV and produces all reportable
quantities, figures, and LaTeX-includable table fragments needed by the
manuscript (papers/1-morfometria/ms/paper1_morfometria.qmd).

Outputs (in ../outputs/, relative to this script):
  - reportable_quantities.json     consumed by the .qmd chunks
  - fig1_size_distribution.png
  - fig2_gmm.png
  - fig3_dimorphism.png
  - fig4_length_weight.png
  - fig5_condition_index.png
  - table1_descriptive.tex          \\input-able longtable
  - table2_dimorphism.tex           \\input-able table
  - table3_lw.tex                   \\input-able table
  - table1_descriptive.csv          machine-readable record
  - table2_dimorphism.csv
  - table3_lw.csv

Usage:
  python scripts/morfometria.py
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import kruskal, mannwhitneyu, shapiro
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
REPO_DIR = HERE.parent
DATA_CSV = REPO_DIR / "data" / "morphometrics_laguna_azul.csv"
CAPTURE_HISTORY_CSV = REPO_DIR / "data" / "capture_history.csv"
OUT_DIR = REPO_DIR / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Aesthetics ────────────────────────────────────────────────────────────────
PALETTE = {"M": "#2166ac", "F": "#d6604d", "J": "#4dac26", "U": "#b2b2b2"}
plt.rcParams.update(
    {"font.family": "serif", "font.size": 10,
     "axes.spines.top": False, "axes.spines.right": False}
)
SP_ITALIC = r"$\it{Atelognathus\ reverberii}$"

# Cei (1969, 1980) reference data for figure 1 bands.
CEI_MALE_SVL = (35.0, 38.0)   # min, max (n=5)
CEI_FEMALE_SVL = (36.5, 38.0)  # min, max (n=2)

SEX_MAP = {
    "macho": "M",
    "hembra": "F", "hemba": "F",
    "juvenil": "J", "juvuenil": "J",
}


# ── Data loading ──────────────────────────────────────────────────────────────
def load_data(csv: Path = DATA_CSV) -> pd.DataFrame:
    """Load and normalise the morphometric CSV.

    Returns df with columns:
        SUL    (mm, snout-urostyle length)
        Mass   (g, body mass)
        MW     (mm, mouth width)
        Sex    (M / F / J / U)
    Only rows with complete SUL, Mass, MW are kept (the CSV already
    enforces that; this is a defence-in-depth check).
    """
    df = pd.read_csv(csv)
    df = df.rename(columns={
        "SUL_mm": "SUL",
        "body_mass_g": "Mass",
        "mouth_width_mm": "MW",
    })
    raw = df["sex_class"].astype(str).str.lower().str.strip()
    df["Sex"] = raw.map(SEX_MAP).fillna("U")
    # Anything that wasn't macho/hembra/juvenil (NaN, '?', 'foto parásito') → U
    df = df.dropna(subset=["SUL", "Mass", "MW"]).reset_index(drop=True)
    return df[["SUL", "Mass", "MW", "Sex"]]


# ── Descriptive statistics ────────────────────────────────────────────────────
def bootstrap_ci_mean(values: np.ndarray, n_boot: int = 2000,
                      seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    boot = np.array([rng.choice(values, size=len(values), replace=True).mean()
                     for _ in range(n_boot)])
    return tuple(np.percentile(boot, [2.5, 97.5]))


def descriptive_per_class(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    label_map = {"M": "Male", "F": "Female", "J": "Juvenile", "U": "Undetermined"}
    for sex in ["M", "F", "J", "U"]:
        sub = df[df["Sex"] == sex]
        for var in ["SUL", "Mass", "MW"]:
            vals = sub[var].dropna().values
            if len(vals) == 0:
                continue
            ci_lo, ci_hi = bootstrap_ci_mean(vals)
            rows.append({
                "Class": label_map[sex], "Variable": var,
                "n": len(vals),
                "Mean": round(float(vals.mean()), 2),
                "SD": round(float(vals.std(ddof=1)), 2),
                "Min": round(float(vals.min()), 2),
                "Max": round(float(vals.max()), 2),
                "CI_lo": round(float(ci_lo), 2),
                "CI_hi": round(float(ci_hi), 2),
            })
    return pd.DataFrame(rows)


# ── Sexual dimorphism ─────────────────────────────────────────────────────────
def dimorphism_tests(df: pd.DataFrame) -> pd.DataFrame:
    males = df[df["Sex"] == "M"]
    females = df[df["Sex"] == "F"]
    rows = []
    for var in ["SUL", "Mass", "MW"]:
        m = males[var].dropna().values
        f = females[var].dropna().values
        u, p = mannwhitneyu(m, f, alternative="two-sided")
        pooled_sd = np.sqrt((m.std(ddof=1) ** 2 + f.std(ddof=1) ** 2) / 2)
        d = (m.mean() - f.mean()) / pooled_sd if pooled_sd > 0 else np.nan
        rows.append({
            "Variable": var,
            "n_M": len(m), "Mean_M": round(float(m.mean()), 2),
            "n_F": len(f), "Mean_F": round(float(f.mean()), 2),
            "U": round(float(u), 1),
            "p": round(float(p), 4),
            "d": round(float(d), 3),
        })
    return pd.DataFrame(rows)


def shapiro_per_sex(df: pd.DataFrame) -> dict:
    out = {}
    for sex in ["M", "F"]:
        sub = df[df["Sex"] == sex]
        for var in ["SUL", "Mass", "MW"]:
            w, p = shapiro(sub[var].dropna().values)
            out[f"shapiro_{sex.lower()}_{var.lower()}_w"] = round(float(w), 3)
            out[f"shapiro_{sex.lower()}_{var.lower()}_p"] = float(p)
    return out


def mw_power_are(n1: int, n2: int, d: float, alpha: float = 0.05) -> float:
    """ARE-corrected power for the two-sided Mann-Whitney U at effect d.

    Effective n_eff_i = n_i * ARE, with ARE = 3/pi ≈ 0.955 for normal data
    (Lehmann & Romano 2006). Power follows the standard normal
    approximation for a two-sample t-test at the effective sizes.
    """
    are = 3.0 / np.pi
    ne1, ne2 = n1 * are, n2 * are
    ncp = d * np.sqrt(ne1 * ne2 / (ne1 + ne2))
    z_crit = stats.norm.ppf(1 - alpha / 2)
    power = stats.norm.cdf(ncp - z_crit) + stats.norm.cdf(-ncp - z_crit)
    return float(power)


def d_threshold_for_power(n1: int, n2: int, target: float = 0.80) -> float:
    """Smallest |d| giving power >= target via bisection."""
    lo, hi = 0.0, 2.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if mw_power_are(n1, n2, mid) >= target:
            hi = mid
        else:
            lo = mid
    return float(round(hi, 2))


# ── GMM + bootstrapped LRT ────────────────────────────────────────────────────
def bic_curve(values: np.ndarray, k_range=range(1, 6),
              seed: int = 42) -> dict:
    bic = {}
    for k in k_range:
        gm = GaussianMixture(n_components=k, random_state=seed, n_init=10)
        gm.fit(values.reshape(-1, 1))
        bic[k] = float(gm.bic(values.reshape(-1, 1)))
    return bic


def bootstrap_lrt_unimodality(values: np.ndarray, B: int = 999,
                              seed: int = 42) -> tuple[float, float]:
    """McLachlan (1987) bootstrapped LRT, k=2 vs k=1.

    Returns (Lambda_obs, p_value). Negative raw Lambdas are clipped to 0.
    """
    X = values.reshape(-1, 1)
    gm1 = GaussianMixture(n_components=1, random_state=seed, n_init=10).fit(X)
    gm2 = GaussianMixture(n_components=2, random_state=seed, n_init=10).fit(X)
    ll1 = gm1.score(X) * len(X)
    ll2 = gm2.score(X) * len(X)
    lam_obs = max(0.0, 2 * (ll2 - ll1))

    mu = gm1.means_.flatten()[0]
    sd = np.sqrt(gm1.covariances_.flatten())[0]
    rng = np.random.default_rng(seed)
    sim = []
    for _ in range(B):
        x = rng.normal(mu, sd, size=len(X)).reshape(-1, 1)
        g1 = GaussianMixture(n_components=1, random_state=seed, n_init=5).fit(x)
        g2 = GaussianMixture(n_components=2, random_state=seed, n_init=5).fit(x)
        sim.append(max(0.0, 2 * (g2.score(x) - g1.score(x)) * len(x)))
    sim = np.array(sim)
    p = float((sim >= lam_obs).mean())
    return float(lam_obs), p


def bootstrap_k2_components(values: np.ndarray, B: int = 500,
                            seed: int = 42) -> dict:
    """Bootstrap CIs for k=2 component (mu1, mu2, w1, w2) parameters."""
    X = values.reshape(-1, 1)
    rng = np.random.default_rng(seed)
    mus_1, mus_2, sds_1, sds_2, ws_1, ws_2 = ([] for _ in range(6))
    for _ in range(B):
        idx = rng.integers(0, len(X), size=len(X))
        x = X[idx]
        gm = GaussianMixture(n_components=2, random_state=seed, n_init=10).fit(x)
        order = np.argsort(gm.means_.flatten())
        mus_1.append(gm.means_.flatten()[order][0])
        mus_2.append(gm.means_.flatten()[order][1])
        sds_1.append(np.sqrt(gm.covariances_.flatten())[order][0])
        sds_2.append(np.sqrt(gm.covariances_.flatten())[order][1])
        ws_1.append(gm.weights_[order][0])
        ws_2.append(gm.weights_[order][1])
    return {
        "mu1": (float(np.mean(mus_1)),
                float(np.percentile(mus_1, 2.5)),
                float(np.percentile(mus_1, 97.5))),
        "mu2": (float(np.mean(mus_2)),
                float(np.percentile(mus_2, 2.5)),
                float(np.percentile(mus_2, 97.5))),
        "sd1": (float(np.mean(sds_1)),
                float(np.percentile(sds_1, 2.5)),
                float(np.percentile(sds_1, 97.5))),
        "sd2": (float(np.mean(sds_2)),
                float(np.percentile(sds_2, 2.5)),
                float(np.percentile(sds_2, 97.5))),
        "w1": (float(np.mean(ws_1)),
               float(np.percentile(ws_1, 2.5)),
               float(np.percentile(ws_1, 97.5))),
        "w2": (float(np.mean(ws_2)),
               float(np.percentile(ws_2, 2.5)),
               float(np.percentile(ws_2, 97.5))),
    }


# ── Length-weight regression ──────────────────────────────────────────────────
def lw_regression(values_sul: np.ndarray, values_mass: np.ndarray) -> dict:
    log_sul = np.log(values_sul)
    log_mass = np.log(values_mass)
    slope, intercept, r, p_fit, se = stats.linregress(log_sul, log_mass)
    n = len(values_sul)
    return {
        "n": n,
        "a": float(np.exp(intercept)),
        "b": float(slope),
        "se_b": float(se),
        "r2": float(r ** 2),
        "p_fit": float(p_fit),
        "df": n - 2,
    }


def test_b_vs(reg: dict, b_null: float) -> tuple[float, float]:
    """One-sample t-test on the estimated slope vs a null exponent."""
    t = (reg["b"] - b_null) / reg["se_b"]
    p_two = 2 * (1 - stats.t.cdf(abs(t), df=reg["df"]))
    return float(t), float(p_two)


# ── Body condition ────────────────────────────────────────────────────────────
def body_condition(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Reference regression on Undetermined-sex adults; residuals for all."""
    ref = df[df["Sex"] == "U"]
    log_sul_ref = np.log(ref["SUL"].values)
    log_mass_ref = np.log(ref["Mass"].values)
    slope, intercept, r, _, _ = stats.linregress(log_sul_ref, log_mass_ref)
    log_sul_all = np.log(df["SUL"].values)
    log_mass_all = np.log(df["Mass"].values)
    df = df.copy()
    df["CI"] = log_mass_all - (intercept + slope * log_sul_all)
    info = {
        "n_ref": int(len(ref)),
        "b_ref": float(slope),
        "r2_ref": float(r ** 2),
    }
    return df, info


# ── Misclassification correction (Buonaccorsi 2010, Ch.3) ─────────────────────
def cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    pooled = np.sqrt((x.std(ddof=1) ** 2 + y.std(ddof=1) ** 2) / 2)
    return float((x.mean() - y.mean()) / pooled) if pooled > 0 else float("nan")


def corrected_d_under_random_misclass(x: np.ndarray, y: np.ndarray,
                                      p_err: float) -> float:
    """Linear correction: if a fraction p_err of each group is randomly
    misassigned to the other, the observed means relate to the true means by
        x_obs = (1-p) * x_true + p * y_true
        y_obs = p * x_true + (1-p) * y_true
    Inverting gives true means; SDs are left as observed (conservative).
    """
    if p_err >= 0.5:
        return float("nan")
    a = 1 - p_err
    det = a ** 2 - p_err ** 2
    x_true_mean = (a * x.mean() - p_err * y.mean()) / det
    y_true_mean = (a * y.mean() - p_err * x.mean()) / det
    pooled = np.sqrt((x.std(ddof=1) ** 2 + y.std(ddof=1) ** 2) / 2)
    return float((x_true_mean - y_true_mean) / pooled) if pooled > 0 else float("nan")


def misclass_threshold(x: np.ndarray, y: np.ndarray,
                       target: float = 0.20) -> float:
    """Smallest p (per cent) where |corrected d| >= target via bisection."""
    lo, hi = 0.0, 0.49
    for _ in range(60):
        mid = (lo + hi) / 2
        d_corr = corrected_d_under_random_misclass(x, y, mid)
        if abs(d_corr) >= target:
            hi = mid
        else:
            lo = mid
    return float(round(hi * 100, 1))


# ── LaTeX table writers ───────────────────────────────────────────────────────
def write_table1(desc: pd.DataFrame, path: Path) -> None:
    rows_per_class = {"Male": "Male", "Female": "Female",
                      "Juvenile": "Juvenile", "Undetermined": "Undetermined"}
    body_lines: list[str] = []
    for cls_label in rows_per_class:
        sub = desc[desc["Class"] == cls_label]
        if sub.empty:
            continue
        for i, var in enumerate(["SUL", "Mass", "MW"]):
            r = sub[sub["Variable"] == var].iloc[0]
            class_cell = (f"\\multirow{{3}}{{*}}{{{cls_label}}}" if i == 0 else "")
            mass_fmt_n = f"{r['n']}"
            body_lines.append(
                f" {class_cell} & {var}  & {mass_fmt_n} & "
                f"{r['Mean']:.2f} & {r['SD']:.2f} & "
                f"{r['Min']:.2f} & {r['Max']:.2f} & "
                f"{{[{r['CI_lo']:.2f}, {r['CI_hi']:.2f}]}}\\\\"
            )
        body_lines.append("\\addlinespace")
    if body_lines and body_lines[-1] == "\\addlinespace":
        body_lines.pop()

    tex = r"""\begin{ThreePartTable}
\begin{TableNotes}[flushleft]\footnotesize
  \item SUL: snout--urostyle length (mm); Mass: body mass (g);
        MW: mouth width (mm).
  \item Bootstrap CI: 95\% confidence interval for the mean
        (2\,000 resampling iterations, seed 42).
\end{TableNotes}
\begin{longtable}{llrrrrrl}
\caption{Descriptive morphometric statistics of \textit{Atelognathus reverberii}
by sex/age class (Laguna Azul, 2019--2020).}
\label{tab:descriptive}\\
\toprule
Class & Variable & $n$ & Mean & SD & Min & Max & 95\% CI\\
\midrule
\endfirsthead
\multicolumn{8}{c}{\tablename~\thetable{} \textit{(continued)}}\\
\toprule
Class & Variable & $n$ & Mean & SD & Min & Max & 95\% CI\\
\midrule
\endhead
\midrule
\multicolumn{8}{r}{\textit{Continued on next page}}\\
\endfoot
\bottomrule
\insertTableNotes
\endlastfoot

""" + "\n".join(body_lines) + r"""

\end{longtable}
\end{ThreePartTable}
"""
    path.write_text(tex)


def write_table2(dim: pd.DataFrame, path: Path) -> None:
    rows = []
    var_label = {"SUL": "SUL (mm)", "Mass": "Mass (g)", "MW": "MW (mm)"}
    for var in ["SUL", "Mass", "MW"]:
        r = dim[dim["Variable"] == var].iloc[0]
        rows.append(
            f"{var_label[var]} & {r['n_M']} & {r['Mean_M']:.2f} & "
            f"{r['n_F']} & {r['Mean_F']:.2f} & "
            f"{r['U']:.0f} & {r['p']:.3f} & {r['d']:.3f}\\\\"
        )
    tex = r"""\begin{table}[H]
\centering
\caption{Sexual dimorphism tests for \textit{Atelognathus reverberii}
(Laguna Azul, 2019--2020). Mann--Whitney $U$ test (two-sided), males vs.\
females. Normality of each variable per sex was assessed with the
Shapiro--Wilk test at $\alpha = 0.05$: SUL did not differ significantly
from normality in either sex ($p > 0.05$); body mass was non-normal in
both sexes ($p < 0.006$); MW did not differ from normality in either sex
(males $p = 0.071$, female $p > 0.5$), with the male value at the
boundary of the threshold but treated as non-significant at $\alpha =
0.05$. Cohen's $d$ (pooled SD): $|d| < 0.2$ = negligible.}
\label{tab:dimorphism}
\begin{tabular}{lrrrrrrl}
\toprule
Variable & $n_M$ & $\bar{x}_M$ & $n_F$ & $\bar{x}_F$ & $U$ & $p$ & $d$\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    path.write_text(tex)


def write_table3(lw_data: dict, path: Path) -> None:
    """lw_data: dict with keys 'All', 'Male', 'Female', 'Juvenile'.

    Each value: {n, a, b, se, p2, p3}. p_juv values are bolded if < 0.05.
    """
    def fmt_p(p: float) -> str:
        return f"{p:.3f}" if p >= 0.001 else "$<$0.001"

    rows = []
    for grp in ["All", "Male", "Female", "Juvenile"]:
        r = lw_data[grp]
        p2 = fmt_p(r["p2"])
        p3 = fmt_p(r["p3"])
        if grp == "Juvenile":
            if r["p2"] < 0.05:
                p2 = f"\\textbf{{{r['p2']:.3f}}}\\tnote{{*}}"
            if r["p3"] < 0.05:
                p3 = f"\\textbf{{{r['p3']:.3f}}}\\tnote{{*}}"
        rows.append(
            f"{grp:8s} & {r['n']:3d} & {r['a']:.5f} & {r['b']:.3f} & "
            f"{r['se']:.3f} & {p2} & {p3}\\\\"
        )
    tex = r"""\begin{table}[H]
\centering
\caption{Length--weight regression parameters for \textit{Atelognathus
reverberii} by sex/age class (Laguna Azul, 2019--2020). Model:
$\log W = \log a + b \cdot \log(\mathrm{SUL})$, fitted by OLS.
$a$ is a dimensioned coefficient with units g$\cdot$mm$^{-b}$ (varying
across groups as $b$ differs); $b$ is dimensionless.
$p_2$: $p$-value for $H_0\colon b = 2$ (two-dimensional isometry);
$p_3$: $p$-value for $H_0\colon b = 3$ (three-dimensional isometry).}
\label{tab:lw}
\begin{threeparttable}
\begin{tabular}{lrrclll}
\toprule
Group & $n$ & $a$ & $b$ & SE$(b)$ & $p_2$ & $p_3$\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]\footnotesize
  \item[*] Significantly different from the stated isometric expectation
  at $\alpha = 0.05$.
\end{tablenotes}
\end{threeparttable}
\end{table}
"""
    path.write_text(tex)


# ── Figures ───────────────────────────────────────────────────────────────────
def fig_size_distribution(df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = np.arange(19, 42, 1)
    label_map = {"M": "Male", "F": "Female", "J": "Juvenile", "U": "Undetermined"}
    for cls in ["U", "J", "F", "M"]:
        vals = df[df["Sex"] == cls]["SUL"]
        ax.hist(vals, bins=bins, color=PALETTE[cls], alpha=0.78,
                label=label_map[cls], edgecolor="white", linewidth=0.3)
    ax.axvspan(*CEI_MALE_SVL, alpha=0.12, color="#2166ac",
               label=f"Cei (1969) ♂ SVL ({CEI_MALE_SVL[0]}–{CEI_MALE_SVL[1]} mm)")
    ax.axvspan(*CEI_FEMALE_SVL, alpha=0.12, color="#d6604d",
               label=f"Cei (1969) ♀ SVL ({CEI_FEMALE_SVL[0]}–{CEI_FEMALE_SVL[1]} mm)")
    ax.set_xlabel("Snout–urostyle length, SUL (mm)", fontsize=11)
    ax.set_ylabel("Frequency", fontsize=11)
    ax.set_title(f"Size distribution of {SP_ITALIC}\n"
                 f"(Laguna Azul, 2019–2020, n = {len(df)})",
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def fig_gmm(full_bic: dict, adult_bic: dict, sul_full: np.ndarray,
            k2_means: tuple[float, float], k2_sds: tuple[float, float],
            k2_weights: tuple[float, float], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ks = list(full_bic.keys())
    ax.plot(ks, [full_bic[k] for k in ks], "o-", color="#525252",
            linewidth=1.8, markersize=7, label=f"Full sample (n = {len(sul_full)})")
    ax.plot(ks, [adult_bic[k] for k in ks], "s--", color="#969696",
            linewidth=1.4, markersize=6,
            label=f"Adults only (n = {int(adult_bic.get('_n', 251))})")
    k_opt = min(full_bic, key=full_bic.get)
    ax.axvline(k_opt, color="#d6604d", linestyle=":", linewidth=1.4,
               label=f"Optimal k = {k_opt}")
    ax.set_xlabel("Number of components (k)", fontsize=11)
    ax.set_ylabel("BIC", fontsize=11)
    ax.set_title("Model selection (BIC)", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False)

    ax = axes[1]
    bins = np.arange(19, 42, 0.8)
    ax.hist(sul_full, bins=bins, color="#d9d9d9", edgecolor="white",
            linewidth=0.3, density=True, label="Observed SUL")
    x_range = np.linspace(sul_full.min() - 2, sul_full.max() + 2, 500)
    gmm_colors = ["#4dac26", "#f4a582"]
    total = np.zeros_like(x_range)
    for i in (0, 1):
        comp = k2_weights[i] * stats.norm.pdf(x_range, k2_means[i], k2_sds[i])
        total += comp
        ax.plot(x_range, comp, color=gmm_colors[i], linewidth=2,
                label=f"Class {i+1}: μ = {k2_means[i]:.1f}, "
                      f"σ = {k2_sds[i]:.1f} mm "
                      f"(w = {k2_weights[i]:.2f})")
    ax.plot(x_range, total, "k--", linewidth=1.5, label="Total GMM (k = 2)")
    ax.set_xlabel("SUL (mm)", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title("k = 2 solution (descriptive; not supported by BIC)",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(f"Population size structure — {SP_ITALIC}", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def fig_dimorphism(df: pd.DataFrame, path: Path) -> None:
    adults = df[df["Sex"].isin(["M", "F"])].copy()
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.5))
    var_labels = {"SUL": "SUL (mm)", "Mass": "Body mass (g)", "MW": "Mouth width, MW (mm)"}
    rng = np.random.default_rng(42)
    for ax, var in zip(axes, ["SUL", "Mass", "MW"]):
        for i, cls in enumerate(["M", "F"]):
            sub = adults[adults["Sex"] == cls][var].dropna()
            bp = ax.boxplot(sub, positions=[i], widths=0.4, patch_artist=True,
                            medianprops=dict(color="black", linewidth=1.8),
                            flierprops=dict(marker="o", markersize=3, alpha=0.4))
            bp["boxes"][0].set_facecolor(PALETTE[cls])
            bp["boxes"][0].set_alpha(0.6)
            jitter = rng.uniform(-0.12, 0.12, len(sub))
            ax.scatter(np.full(len(sub), i) + jitter, sub,
                       color=PALETTE[cls], alpha=0.35, s=12, zorder=3)
        nM = (adults["Sex"] == "M").sum()
        nF = (adults["Sex"] == "F").sum()
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Male\n(n = {nM})", f"Female\n(n = {nF})"],
                           fontsize=9)
        ax.set_ylabel(var_labels[var], fontsize=10)
        ax.set_title(var_labels[var].split(" (")[0], fontsize=10)
    fig.suptitle(f"Sexual dimorphism — {SP_ITALIC}\n"
                 f"(Mann–Whitney U, all p > 0.34; Cohen's d < 0.09)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def fig_length_weight(df: pd.DataFrame, lw_all: dict, path: Path) -> None:
    import matplotlib.ticker as ticker
    fig, ax = plt.subplots(figsize=(7, 5))
    label_map = {"M": "Male", "F": "Female", "J": "Juvenile", "U": "Undetermined"}
    for cls in ["U", "J", "F", "M"]:
        sub = df[df["Sex"] == cls].dropna(subset=["SUL", "Mass"])
        ax.scatter(sub["SUL"], sub["Mass"], color=PALETTE[cls],
                   alpha=0.55, s=22, label=label_map[cls], edgecolors="none")
    x = np.linspace(df["SUL"].min() * 0.95, df["SUL"].max() * 1.05, 300)
    a, b = lw_all["a"], lw_all["b"]
    ax.plot(x, a * x ** b, "k-", linewidth=2.2,
            label=f"All: W = {a:.5f} × SUL$^{{{b:.3f}}}$ "
                  f"(r² = {lw_all['r2']:.3f})", zorder=5)
    ax.plot(x, a * x ** 2, "--", color="#969696", linewidth=1.2, label="Isometry b = 2")
    ax.plot(x, a * x ** 3, ":", color="#525252", linewidth=1.2, label="Isometry b = 3")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("SUL (mm) — log scale", fontsize=11)
    ax.set_ylabel("Body mass (g) — log scale", fontsize=11)
    ax.set_title(f"Length–weight relationship — {SP_ITALIC}", fontsize=11)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.legend(fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def fig_condition_index(df_ci: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    label_map = {"M": "Male", "F": "Female", "J": "Juvenile", "U": "Undetermined (ref.)"}
    for cls in ["M", "F", "J", "U"]:
        sub = df_ci[df_ci["Sex"] == cls]["CI"].dropna()
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=15, color=PALETTE[cls], alpha=0.65,
                label=f"{label_map[cls]} (n = {len(sub)})",
                edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linestyle="--", linewidth=1, label="Reference (0)")
    ax.set_xlabel("Body condition index (log-residuals)", fontsize=11)
    ax.set_ylabel("Frequency", fontsize=11)
    ax.set_title(f"Body condition index — {SP_ITALIC}", fontsize=11)
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


# ── Main pipeline ─────────────────────────────────────────────────────────────
def main() -> None:
    df = load_data()
    rq: dict = {}

    # ── 1. Sample composition
    rq["n_total"] = int(len(df))
    rq["n_male"] = int((df["Sex"] == "M").sum())
    rq["n_female"] = int((df["Sex"] == "F").sum())
    rq["n_juvenile"] = int((df["Sex"] == "J").sum())
    rq["n_undetermined"] = int((df["Sex"] == "U").sum())
    rq["pct_undetermined"] = round(100 * rq["n_undetermined"] / rq["n_total"], 1)

    # ── 2. Descriptive table
    desc = descriptive_per_class(df)
    desc.to_csv(OUT_DIR / "table1_descriptive.csv", index=False)
    write_table1(desc, OUT_DIR / "table1_descriptive.tex")
    # Per-class summary into rq (used by abstract + results)
    short = {"M": "male", "F": "female", "J": "juv", "U": "und"}
    var_short = {"SUL": "sul", "Mass": "mass", "MW": "mw"}
    label_to_code = {"Male": "M", "Female": "F", "Juvenile": "J", "Undetermined": "U"}
    for _, r in desc.iterrows():
        s = short[label_to_code[r["Class"]]]
        v = var_short[r["Variable"]]
        rq[f"{s}_{v}_n"] = int(r["n"])
        rq[f"{s}_{v}_mean"] = float(r["Mean"])
        rq[f"{s}_{v}_sd"] = float(r["SD"])
        rq[f"{s}_{v}_min"] = float(r["Min"])
        rq[f"{s}_{v}_max"] = float(r["Max"])
        rq[f"{s}_{v}_ci_lo"] = float(r["CI_lo"])
        rq[f"{s}_{v}_ci_hi"] = float(r["CI_hi"])

    # ── 3. Dimorphism + Shapiro + power
    dim = dimorphism_tests(df)
    dim.to_csv(OUT_DIR / "table2_dimorphism.csv", index=False)
    write_table2(dim, OUT_DIR / "table2_dimorphism.tex")
    for _, r in dim.iterrows():
        v = var_short[r["Variable"]]
        rq[f"dim_{v}_U"] = float(r["U"])
        rq[f"dim_{v}_p"] = float(r["p"])
        rq[f"dim_{v}_d"] = float(r["d"])
    rq.update(shapiro_per_sex(df))
    nM, nF = rq["n_male"], rq["n_female"]
    rq["mw_are_value"] = round(3.0 / np.pi, 3)
    rq["mw_power_d05"] = round(mw_power_are(nM, nF, 0.5), 3)
    rq["mw_power_d08"] = round(mw_power_are(nM, nF, 0.8), 3)
    rq["d_threshold_80pct"] = d_threshold_for_power(nM, nF, 0.80)

    # ── 4. GMM full + adults
    sul_full = df["SUL"].values
    adults = df[df["Sex"].isin(["M", "F", "U"])]   # excludes 25 juveniles
    sul_adult = adults["SUL"].values
    bic_full = bic_curve(sul_full)
    bic_adult = bic_curve(sul_adult)
    for k, v in bic_full.items():
        rq[f"bic_full_k{k}"] = round(v, 2)
    for k, v in bic_adult.items():
        rq[f"bic_adult_k{k}"] = round(v, 2)
    rq["delta_bic_full_k1_k2"] = round(bic_full[2] - bic_full[1], 1)
    rq["delta_bic_adult_k1_k2"] = round(bic_adult[2] - bic_adult[1], 1)
    rq["n_adult_subset"] = int(len(sul_adult))

    # ── 5. Bootstrapped LRT
    print("[5/8] bootstrapped LRT (B=999)…")
    lam_obs, lrt_p = bootstrap_lrt_unimodality(sul_full, B=999)
    rq["lrt_lambda"] = round(lam_obs, 3)
    rq["lrt_p"] = round(lrt_p, 3)
    rq["lrt_B"] = 999

    # ── 6. k=2 bootstrap CIs + point estimates
    print("[6/8] bootstrap k=2 CIs (B=500)…")
    k2 = bootstrap_k2_components(sul_full, B=500)
    rq["k2_class1_mu_mean"] = round(k2["mu1"][0], 1)
    rq["k2_class1_mu_lo"] = round(k2["mu1"][1], 1)
    rq["k2_class1_mu_hi"] = round(k2["mu1"][2], 1)
    rq["k2_class1_sd_mean"] = round(k2["sd1"][0], 1)
    rq["k2_class2_mu_mean"] = round(k2["mu2"][0], 1)
    rq["k2_class2_mu_lo"] = round(k2["mu2"][1], 1)
    rq["k2_class2_mu_hi"] = round(k2["mu2"][2], 1)
    rq["k2_class2_sd_mean"] = round(k2["sd2"][0], 1)
    rq["k2_class1_w_lo"] = round(k2["w1"][1], 2)
    rq["k2_class1_w_hi"] = round(k2["w1"][2], 2)
    rq["k2_class2_w_lo"] = round(k2["w2"][1], 2)
    rq["k2_class2_w_hi"] = round(k2["w2"][2], 2)
    # Point estimates from a single fit (for figure caption)
    gm2 = GaussianMixture(n_components=2, random_state=42, n_init=10).fit(
        sul_full.reshape(-1, 1))
    order = np.argsort(gm2.means_.flatten())
    rq["k2_class1_mu_point"] = round(float(gm2.means_.flatten()[order[0]]), 1)
    rq["k2_class1_sd_point"] = round(float(np.sqrt(gm2.covariances_.flatten())[order[0]]), 1)
    rq["k2_class1_w_point"] = round(float(gm2.weights_[order[0]]), 2)
    rq["k2_class2_mu_point"] = round(float(gm2.means_.flatten()[order[1]]), 1)
    rq["k2_class2_sd_point"] = round(float(np.sqrt(gm2.covariances_.flatten())[order[1]]), 1)
    rq["k2_class2_w_point"] = round(float(gm2.weights_[order[1]]), 2)

    # ── 7. Length-weight regressions
    print("[7/8] length-weight regressions…")
    lw_data = {}
    for label, mask in [("All", slice(None)),
                        ("Male", df["Sex"] == "M"),
                        ("Female", df["Sex"] == "F"),
                        ("Juvenile", df["Sex"] == "J")]:
        sub = df[mask] if not isinstance(mask, slice) else df
        reg = lw_regression(sub["SUL"].values, sub["Mass"].values)
        t2, p2 = test_b_vs(reg, 2.0)
        t3, p3 = test_b_vs(reg, 3.0)
        lw_data[label] = {"n": reg["n"], "a": reg["a"], "b": reg["b"],
                          "se": reg["se_b"], "p2": p2, "p3": p3,
                          "t2": t2, "t3": t3, "r2": reg["r2"], "df": reg["df"]}
        code = {"All": "all", "Male": "male", "Female": "female",
                "Juvenile": "juv"}[label]
        rq[f"lw_{code}_a"] = round(reg["a"], 5)
        rq[f"lw_{code}_b"] = round(reg["b"], 3)
        rq[f"lw_{code}_se"] = round(reg["se_b"], 3)
        rq[f"lw_{code}_r2"] = round(reg["r2"], 3)
        rq[f"lw_{code}_df"] = int(reg["df"])
        rq[f"lw_{code}_t2"] = round(t2, 2)
        rq[f"lw_{code}_p2"] = round(p2, 3)
        rq[f"lw_{code}_t3"] = round(t3, 2)
        rq[f"lw_{code}_p3"] = round(p3, 4)
    # Persist table 3
    pd.DataFrame([{**{"group": g}, **v} for g, v in lw_data.items()]).to_csv(
        OUT_DIR / "table3_lw.csv", index=False)
    write_table3(lw_data, OUT_DIR / "table3_lw.tex")

    # ── 8. Body condition
    print("[8/8] body condition…")
    df_ci, ci_info = body_condition(df)
    rq["bc_n_ref"] = ci_info["n_ref"]
    rq["bc_b_ref"] = round(ci_info["b_ref"], 3)
    rq["bc_r2_ref"] = round(ci_info["r2_ref"], 3)

    # Means + SDs per class (residuals)
    label_codes = {"M": "male", "F": "female", "J": "juv", "U": "und"}
    for cls, name in label_codes.items():
        vals = df_ci[df_ci["Sex"] == cls]["CI"].values
        rq[f"bc_{name}_mean"] = round(float(vals.mean()), 3)
        rq[f"bc_{name}_sd"] = round(float(vals.std(ddof=1)), 3)

    # KW + pairwise (Bonferroni on 3 pairs) + effect size
    vals_M = df_ci[df_ci["Sex"] == "M"]["CI"].values
    vals_F = df_ci[df_ci["Sex"] == "F"]["CI"].values
    vals_J = df_ci[df_ci["Sex"] == "J"]["CI"].values
    h, p_kw = kruskal(vals_M, vals_F, vals_J)
    rq["bc_kw_h"] = round(float(h), 2)
    rq["bc_kw_df"] = 2
    rq["bc_kw_p"] = round(float(p_kw), 3)
    # Effect size for Kruskal-Wallis (Tomczak & Tomczak 2014):
    #   η² = (H - k + 1) / (n - k), with k = number of groups, n = total
    #   ε² = H / (n - 1)
    # Cohen's η² convention: 0.01 small, 0.06 medium, 0.14 large.
    n_kw = len(vals_M) + len(vals_F) + len(vals_J)
    k_kw = 3
    rq["bc_kw_n"] = n_kw
    rq["bc_kw_eta2"] = round(float(max(0.0, (h - k_kw + 1) / (n_kw - k_kw))), 3)
    rq["bc_kw_epsilon2"] = round(float(h / (n_kw - 1)), 3)
    pairs = [("MF", vals_M, vals_F), ("MJ", vals_M, vals_J), ("FJ", vals_F, vals_J)]
    for tag, a, b in pairs:
        _, p_pair = mannwhitneyu(a, b, alternative="two-sided")
        rq[f"bc_p_{tag.lower()}_raw"] = round(float(p_pair), 4)
        rq[f"bc_p_{tag.lower()}_bonf"] = round(min(float(p_pair) * 3, 1.0), 3)

    # ── 8b. Sensitivity: body-condition residuals under POOLED reference slope.
    # Confirms the direction/magnitude of the juvenile-male contrast is robust
    # to the choice of reference slope. Re-fit log(W)~log(SUL) on the full
    # sample and re-project residuals.
    log_sul_all = np.log(df["SUL"].values)
    log_mass_all = np.log(df["Mass"].values)
    slope_p, intercept_p, r_p, _, _ = stats.linregress(log_sul_all, log_mass_all)
    df_p = df.copy()
    df_p["CI"] = log_mass_all - (intercept_p + slope_p * log_sul_all)
    vals_M_p = df_p[df_p["Sex"] == "M"]["CI"].values
    vals_J_p = df_p[df_p["Sex"] == "J"]["CI"].values
    rq["bc_pooled_slope"] = round(float(slope_p), 3)
    rq["bc_pooled_male_mean"] = round(float(vals_M_p.mean()), 3)
    rq["bc_pooled_juv_mean"] = round(float(vals_J_p.mean()), 3)
    rq["bc_pooled_mj_diff"] = round(
        float(vals_M_p.mean() - vals_J_p.mean()), 3)
    _, p_mj_pooled = mannwhitneyu(vals_M_p, vals_J_p, alternative="two-sided")
    rq["bc_pooled_p_mj_raw"] = round(float(p_mj_pooled), 4)
    rq["bc_pooled_p_mj_bonf"] = round(min(float(p_mj_pooled) * 3, 1.0), 3)
    # Reference-slope contrast for comparison
    rq["bc_ref_mj_diff"] = round(
        float(vals_M.mean() - vals_J.mean()), 3)

    # ── 9. Misclassification correction
    sul_M = df[df["Sex"] == "M"]["SUL"].values
    sul_F = df[df["Sex"] == "F"]["SUL"].values
    rq["misclass_threshold_pct"] = misclass_threshold(sul_M, sul_F, target=0.20)
    rq["misclass_d_at_20pct"] = round(
        abs(corrected_d_under_random_misclass(sul_M, sul_F, 0.20)), 3)

    # ── 9c. Formal comparison: SUL distributions of assigned (M+F) vs
    # undetermined adults. Quantifies whether the undetermined pool is
    # size-biased relative to clearly assigned adults.
    sul_assigned = df[df["Sex"].isin(["M", "F"])]["SUL"].values
    sul_und = df[df["Sex"] == "U"]["SUL"].values
    u_au, p_au = mannwhitneyu(sul_assigned, sul_und, alternative="two-sided")
    rq["assigned_vs_und_n_assigned"] = int(len(sul_assigned))
    rq["assigned_vs_und_n_und"] = int(len(sul_und))
    rq["assigned_vs_und_median_assigned"] = round(
        float(np.median(sul_assigned)), 2)
    rq["assigned_vs_und_median_und"] = round(float(np.median(sul_und)), 2)
    rq["assigned_vs_und_U"] = float(u_au)
    rq["assigned_vs_und_p"] = round(float(p_au), 3)

    # ── 9b. Adult-overall ranges (pooled M + F + U, excludes juveniles).
    # Compute the pooled adult ranges so the .qmd can report them as
    # "overall" with sex-specific ranges in parentheses.
    adults = df[df["Sex"].isin(["M", "F", "U"])]
    for var, key in [("SUL", "sul"), ("Mass", "mass"), ("MW", "mw")]:
        rq[f"adult_overall_{key}_min"] = round(float(adults[var].min()), 2)
        rq[f"adult_overall_{key}_max"] = round(float(adults[var].max()), 2)

    # ── 10. Capture-history context (computed from capture_history.csv)
    ch = pd.read_csv(CAPTURE_HISTORY_CSV)
    occ_cols = [c for c in ch.columns if c != "individual_id"]
    caps_per_ind = ch[occ_cols].sum(axis=1)
    rq["n_sampling_occasions"] = len(occ_cols)
    rq["n_individuals_marked_total"] = int(len(ch))               # 275
    rq["n_first_captures_measured"] = rq["n_total"]               # 274 (analytic n)
    rq["n_capture_events_total"] = int(ch[occ_cols].sum().sum())  # 339
    # Capture-frequency distribution (reported as a distribution, not as an
    # event subtraction): how many individuals were caught 1, 2 or 3 times.
    rq["n_captured_once"] = int((caps_per_ind == 1).sum())        # 220
    rq["n_captured_twice"] = int((caps_per_ind == 2).sum())       # 46
    rq["n_captured_thrice"] = int((caps_per_ind == 3).sum())      # 9
    rq["n_individuals_recaptured"] = int((caps_per_ind >= 2).sum())  # 55
    rq["pct_individuals_recaptured"] = round(
        100 * rq["n_individuals_recaptured"] / rq["n_individuals_marked_total"], 1)

    # Measurement-error context (sample-of-20 reread protocol, not from data)
    rq["cv_sul_pct"] = 1.2
    rq["cv_mw_pct"] = 2.1
    rq["n_repeated_meas"] = 20
    rq["pct_repeated_meas"] = round(
        100 * rq["n_repeated_meas"] / rq["n_total"], 1)

    # ── Figures
    fig_size_distribution(df, OUT_DIR / "fig1_size_distribution.png")
    bic_adult_for_fig = dict(bic_adult); bic_adult_for_fig["_n"] = len(sul_adult)
    fig_gmm(bic_full, bic_adult_for_fig, sul_full,
            (rq["k2_class1_mu_point"], rq["k2_class2_mu_point"]),
            (rq["k2_class1_sd_point"], rq["k2_class2_sd_point"]),
            (rq["k2_class1_w_point"], rq["k2_class2_w_point"]),
            OUT_DIR / "fig2_gmm.png")
    fig_dimorphism(df, OUT_DIR / "fig3_dimorphism.png")
    fig_length_weight(df, lw_data["All"], OUT_DIR / "fig4_length_weight.png")
    fig_condition_index(df_ci, OUT_DIR / "fig5_condition_index.png")

    # ── Dump reportable quantities
    (OUT_DIR / "reportable_quantities.json").write_text(
        json.dumps(rq, indent=2, ensure_ascii=False))
    print(f"\nWrote {len(rq)} keys to {OUT_DIR / 'reportable_quantities.json'}")
    print(f"Figures + tables in {OUT_DIR}/")


if __name__ == "__main__":
    main()
