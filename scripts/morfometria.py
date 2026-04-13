"""
morfometria.py
==============
Análisis morfométrico completo de Atelognathus reverberii.
Corresponde al Capítulo II / Paper 1 de la tesis doctoral.

Autores: Kass, N.A.; Kacoliris, F.P.
Uso:
    python morfometria.py --data ../data/Planilla_CMR.xlsx --outdir ../outputs/morfometria
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
from scipy import stats
from scipy.stats import shapiro, mannwhitneyu, kruskal
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).parent))
from data_loader import load_morpho, CEI_DATA

# ── Aesthetics ──────────────────────────────────────────────────────────────
PALETTE = {
    'M': '#2166ac',
    'F': '#d6604d',
    'J': '#4dac26',
    'Unknown': '#b2b2b2'
}
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})


# ── 1. Descriptive statistics ────────────────────────────────────────────────

def descriptive_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a summary table per sex/age class.
    Bootstrap 95% CI for mean (n=2000 iterations).
    """
    rows = []
    label_map = {'M': 'Male', 'F': 'Female', 'J': 'Juvenile', np.nan: 'Unknown'}

    for sex in ['M', 'F', 'J', np.nan]:
        if sex is np.nan:
            sub = df[df['Sexo'].isna()]
        else:
            sub = df[df['Sexo'] == sex]
        if len(sub) == 0:
            continue

        for var in ['LHU', 'Peso', 'AB']:
            vals = sub[var].dropna().values
            if len(vals) == 0:
                continue
            # Bootstrap CI
            rng = np.random.default_rng(42)
            boot_means = [
                rng.choice(vals, size=len(vals), replace=True).mean()
                for _ in range(2000)
            ]
            ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])
            rows.append({
                'Class': label_map[sex],
                'n': len(vals),
                'Variable': var,
                'Mean': round(vals.mean(), 2),
                'SD': round(vals.std(ddof=1), 2),
                'Min': round(vals.min(), 2),
                'Max': round(vals.max(), 2),
                'CI_lo': round(ci_lo, 2),
                'CI_hi': round(ci_hi, 2),
            })

    return pd.DataFrame(rows)


# ── 2. Sexual dimorphism tests ───────────────────────────────────────────────

def sexual_dimorphism_tests(df: pd.DataFrame) -> pd.DataFrame:
    """
    Mann-Whitney U test (males vs females) for each morphometric variable.
    Reports U, p-value, and Cohen's d effect size.
    """
    males = df[df['Sexo'] == 'M']
    females = df[df['Sexo'] == 'F']
    rows = []
    for var in ['LHU', 'Peso', 'AB']:
        m = males[var].dropna().values
        f = females[var].dropna().values
        u_stat, p_val = mannwhitneyu(m, f, alternative='two-sided')
        # Cohen's d
        pooled_sd = np.sqrt((m.std(ddof=1)**2 + f.std(ddof=1)**2) / 2)
        d = (m.mean() - f.mean()) / pooled_sd if pooled_sd > 0 else np.nan
        rows.append({
            'Variable': var,
            'n_males': len(m),
            'Mean_males': round(m.mean(), 2),
            'n_females': len(f),
            'Mean_females': round(f.mean(), 2),
            'U': round(u_stat, 1),
            'p_value': round(p_val, 4),
            'Cohen_d': round(d, 3),
            'Significant': p_val < 0.05
        })
    return pd.DataFrame(rows)


# ── 3. Gaussian Mixture Model ────────────────────────────────────────────────

def fit_gmm(df: pd.DataFrame, n_components_range=(1, 6)) -> dict:
    """
    Fits GMM on LHU using BIC to select optimal number of components.
    Returns best model, BIC scores, and class assignments.
    """
    X = df['LHU'].dropna().values.reshape(-1, 1)
    bic_scores = {}
    models = {}
    for k in range(*n_components_range):
        gm = GaussianMixture(n_components=k, random_state=42, n_init=10)
        gm.fit(X)
        bic_scores[k] = gm.bic(X)
        models[k] = gm

    best_k = min(bic_scores, key=bic_scores.get)
    best_model = models[best_k]
    labels = best_model.predict(X)

    # Sort components by mean (smallest = youngest)
    order = np.argsort(best_model.means_.flatten())
    label_remap = {old: new for new, old in enumerate(order)}
    labels_sorted = np.array([label_remap[l] for l in labels])

    return {
        'model': best_model,
        'best_k': best_k,
        'bic_scores': bic_scores,
        'X': X,
        'labels': labels_sorted,
        'means': best_model.means_.flatten()[order],
        'stds': np.sqrt(best_model.covariances_.flatten())[order],
        'weights': best_model.weights_[order],
    }


# ── 4. Length-weight relationship ────────────────────────────────────────────

def length_weight_regression(df: pd.DataFrame) -> dict:
    """
    Fits log(W) = log(a) + b * log(LHU) by sex class.
    Returns regression parameters and r² per group.
    """
    results = {}
    for sex in ['M', 'F', 'J', 'all']:
        if sex == 'all':
            sub = df.dropna(subset=['LHU', 'Peso'])
        else:
            sub = df[df['Sexo'] == sex].dropna(subset=['LHU', 'Peso'])
        if len(sub) < 5:
            continue
        log_lhu = np.log(sub['LHU'].values)
        log_w = np.log(sub['Peso'].values)
        slope, intercept, r, p, se = stats.linregress(log_lhu, log_w)
        results[sex] = {
            'n': len(sub),
            'a': round(np.exp(intercept), 5),
            'b': round(slope, 4),
            'r2': round(r**2, 4),
            'p': round(p, 6),
        }
    return results


# ── 5. Condition index ────────────────────────────────────────────────────────

def condition_index(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes scaled mass index (SMI) as residuals from log-log regression
    of Peso on LHU (all individuals pooled).
    """
    sub = df.dropna(subset=['LHU', 'Peso']).copy()
    log_lhu = np.log(sub['LHU'].values)
    log_w = np.log(sub['Peso'].values)
    slope, intercept, *_ = stats.linregress(log_lhu, log_w)
    predicted = intercept + slope * log_lhu
    sub['CI'] = log_w - predicted
    return sub


# ── 6. Comparison with Cei ────────────────────────────────────────────────────

def cei_comparison(df: pd.DataFrame) -> dict:
    """
    Compares present sample against Cei (1969, 1980) original description.
    Note: Cei used SVL (to cloaca); present study uses LHU (to urostyle).
    LHU < SVL; direct comparison interpreted with this caveat.
    """
    males = df[df['Sexo'] == 'M']['LHU'].dropna().values
    females = df[df['Sexo'] == 'F']['LHU'].dropna().values

    return {
        'males_present': {
            'n': len(males),
            'mean': round(males.mean(), 2),
            'sd': round(males.std(ddof=1), 2),
            'min': round(males.min(), 2),
            'max': round(males.max(), 2),
        },
        'females_present': {
            'n': len(females),
            'mean': round(females.mean(), 2),
            'sd': round(females.std(ddof=1), 2),
            'min': round(females.min(), 2),
            'max': round(females.max(), 2),
        },
        'cei_males': CEI_DATA['males'],
        'cei_females': CEI_DATA['females'],
        'note': CEI_DATA['note_methodology'],
    }


# ── Shared figure helpers ─────────────────────────────────────────────────────

SP_ITALIC = r'$\it{Atelognathus\ reverberii}$'   # italic species name for matplotlib
SP_SHORT  = r'$\it{A.\ reverberii}$'


# ── Figures ──────────────────────────────────────────────────────────────────

def fig_size_distribution(df: pd.DataFrame, outdir: Path) -> Path:
    """Fig 1: SUL distribution by class + Cei SVL reference bands."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = np.arange(19, 42, 1)
    df2 = df.copy()
    df2['Class'] = df2['Sexo'].fillna('Unknown')

    label_map = {'Unknown': 'Undetermined', 'J': 'Juvenile', 'F': 'Female', 'M': 'Male'}
    for label in ['Unknown', 'J', 'F', 'M']:
        vals = df2[df2['Class'] == label]['LHU']
        ax.hist(vals, bins=bins, color=PALETTE[label], alpha=0.78,
                label=label_map[label], edgecolor='white', linewidth=0.3)

    # Cei SVL reference bands
    ax.axvspan(CEI_DATA['males']['LHU_min'], CEI_DATA['males']['LHU_max'],
               alpha=0.12, color='#2166ac',
               label=f"Cei (1969) \u2642 SVL range "
                     f"({CEI_DATA['males']['LHU_min']}\u2013{CEI_DATA['males']['LHU_max']} mm)")
    ax.axvspan(CEI_DATA['females']['LHU_min'], CEI_DATA['females']['LHU_max'],
               alpha=0.12, color='#d6604d',
               label=f"Cei (1969) \u2640 SVL range "
                     f"({CEI_DATA['females']['LHU_min']}\u2013{CEI_DATA['females']['LHU_max']} mm)")

    ax.set_xlabel('Snout\u2013urostyle length, SUL (mm)', fontsize=11)
    ax.set_ylabel('Frequency', fontsize=11)
    ax.set_title(
        f'Size distribution of {SP_ITALIC}\n(Laguna Azul, 2019\u20132020, n\u202f=\u202f{len(df)})',
        fontsize=11
    )
    ax.legend(fontsize=8, frameon=False)

    fpath = outdir / 'fig1_size_distribution.png'
    fig.tight_layout()
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    return fpath


def fig_gmm(gmm_result: dict, df: pd.DataFrame, outdir: Path,
            adults_bic: dict = None) -> Path:
    """Fig 2: BIC curve (left) + k=2 overlay on SUL histogram (right).
    Left panel shows BIC for k=1..5 for full sample and adults-only.
    Right panel always shows the k=2 solution overlaid for descriptive
    comparison, regardless of the BIC-selected optimal k.
    """
    X = gmm_result['X'].flatten()
    k_opt = gmm_result['best_k']
    gmm_colors = ['#4dac26', '#f4a582', '#2166ac', '#b2182b', '#762a83']

    # Always compute k=2 for the right-panel overlay
    from sklearn.mixture import GaussianMixture as _GMM
    gm2 = _GMM(n_components=2, random_state=42, n_init=10).fit(
        X.reshape(-1, 1))
    order2 = np.argsort(gm2.means_.flatten())
    means2 = gm2.means_.flatten()[order2]
    stds2  = np.sqrt(gm2.covariances_.flatten())[order2]
    weights2 = gm2.weights_[order2]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # ── Left: BIC curves ──────────────────────────────────────────────────
    ax = axes[0]
    ks = list(gmm_result['bic_scores'].keys())
    bics = [gmm_result['bic_scores'][k_] for k_ in ks]
    ax.plot(ks, bics, 'o-', color='#525252', linewidth=1.8, markersize=7,
            label='Full sample (n\u202f=\u202f276)')
    if adults_bic:
        ks_a = list(adults_bic.keys())
        bics_a = [adults_bic[k_] for k_ in ks_a]
        ax.plot(ks_a, bics_a, 's--', color='#969696', linewidth=1.4, markersize=6,
                label='Adults only (n\u202f=\u202f251)')
    ax.axvline(k_opt, color='#d6604d', linestyle=':', linewidth=1.4,
               label=f'Optimal k\u202f=\u202f{k_opt}')
    ax.set_xlabel('Number of components (k)', fontsize=11)
    ax.set_ylabel('BIC', fontsize=11)
    ax.set_title('Model selection (BIC)', fontsize=11)
    ax.legend(fontsize=8.5, frameon=False)

    # ── Right: k=2 overlay (descriptive) ─────────────────────────────────
    ax = axes[1]
    bins = np.arange(19, 42, 0.8)
    ax.hist(X, bins=bins, color='#d9d9d9', edgecolor='white',
            linewidth=0.3, density=True, label='Observed SUL')

    x_range = np.linspace(X.min() - 2, X.max() + 2, 500)
    total_pdf = np.zeros_like(x_range)
    for i in range(2):
        comp_pdf = weights2[i] * stats.norm.pdf(x_range, means2[i], stds2[i])
        total_pdf += comp_pdf
        ax.plot(x_range, comp_pdf, color=gmm_colors[i], linewidth=2,
                label=f'Class {i+1}: \u03bc\u202f=\u202f{means2[i]:.1f}, '
                      f'\u03c3\u202f=\u202f{stds2[i]:.1f}\u202fmm '
                      f'(w\u202f=\u202f{weights2[i]:.2f})')
    ax.plot(x_range, total_pdf, 'k--', linewidth=1.5, label='Total GMM (k\u202f=\u202f2)')
    ax.set_xlabel('SUL (mm)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('k\u202f=\u202f2 solution (descriptive; not supported by BIC)',
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle(
        f'Population size structure \u2014 {SP_ITALIC}',
        fontsize=11
    )
    fpath = outdir / 'fig2_gmm.png'
    fig.tight_layout()
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    return fpath


def fig_dimorphism(df: pd.DataFrame, outdir: Path) -> Path:
    """Fig 3: Boxplots + strip plots for sexual dimorphism (SUL, mass, MW)."""
    adults = df[df['Sexo'].isin(['M', 'F'])].copy()
    adults['Class'] = adults['Sexo'].map({'M': 'Male', 'F': 'Female'})

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.5))
    var_labels = {
        'LHU': 'SUL (mm)',
        'Peso': 'Body mass (g)',
        'AB': 'Mouth width, MW (mm)'
    }

    for ax, var in zip(axes, ['LHU', 'Peso', 'AB']):
        for i, (cls, col) in enumerate([('Male', PALETTE['M']), ('Female', PALETTE['F'])]):
            sub = adults[adults['Class'] == cls][var].dropna()
            bp = ax.boxplot(sub, positions=[i], widths=0.4, patch_artist=True,
                            medianprops=dict(color='black', linewidth=1.8),
                            flierprops=dict(marker='o', markersize=3, alpha=0.4))
            bp['boxes'][0].set_facecolor(col)
            bp['boxes'][0].set_alpha(0.6)
            jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(sub))
            ax.scatter(np.full(len(sub), i) + jitter, sub,
                       color=col, alpha=0.35, s=12, zorder=3)

        nM = len(adults[adults['Class'] == 'Male'][var].dropna())
        nF = len(adults[adults['Class'] == 'Female'][var].dropna())
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'Male\n(n\u202f=\u202f{nM})', f'Female\n(n\u202f=\u202f{nF})'], fontsize=9)
        ax.set_ylabel(var_labels[var], fontsize=10)
        ax.set_title(var_labels[var].split(' (')[0], fontsize=10)

    fig.suptitle(
        f'Sexual dimorphism \u2014 {SP_ITALIC}\n'
        f'(Mann\u2013Whitney U, all p\u202f>\u202f0.34; Cohen\u2019s d\u202f<\u202f0.09)',
        fontsize=10
    )
    fpath = outdir / 'fig3_dimorphism.png'
    fig.tight_layout()
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    return fpath


def fig_length_weight(df: pd.DataFrame, lw_results: dict, outdir: Path) -> Path:
    """Fig 4: Length-weight relationship on log-log scale with isometric references."""
    fig, ax = plt.subplots(figsize=(7, 5))
    df2 = df.copy()
    df2['Class'] = df2['Sexo'].fillna('Unknown')

    label_map = {'M': 'Male', 'F': 'Female', 'J': 'Juvenile', 'Unknown': 'Undetermined'}
    for cls in ['Unknown', 'J', 'F', 'M']:
        sub = df2[df2['Class'] == cls].dropna(subset=['LHU', 'Peso'])
        ax.scatter(sub['LHU'], sub['Peso'], color=PALETTE[cls],
                   alpha=0.55, s=22, label=label_map[cls], edgecolors='none')

    x_range = np.linspace(df['LHU'].min() * 0.95, df['LHU'].max() * 1.05, 300)

    # Fitted regression (all individuals)
    if 'all' in lw_results:
        r = lw_results['all']
        ax.plot(x_range, r['a'] * x_range ** r['b'], 'k-', linewidth=2.2,
                label=f"All: W = {r['a']} \u00d7 SUL$^{{{r['b']}}}$ "
                      f"(r\u00b2\u202f=\u202f{r['r2']})",
                zorder=5)

    # Isometric references anchored to same intercept
    if 'all' in lw_results:
        a_ref = lw_results['all']['a']
        ax.plot(x_range, a_ref * x_range ** 2, '--', color='#969696',
                linewidth=1.2, label='Isometry b\u202f=\u202f2')
        ax.plot(x_range, a_ref * x_range ** 3, ':', color='#525252',
                linewidth=1.2, label='Isometry b\u202f=\u202f3')

    # Log-log scale (natural scale for L-W relationships)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('SUL (mm) — log scale', fontsize=11)
    ax.set_ylabel('Body mass (g) — log scale', fontsize=11)
    ax.set_title(
        f'Length\u2013weight relationship \u2014 {SP_ITALIC}',
        fontsize=11
    )

    # Clean tick labels on log axes
    import matplotlib.ticker as ticker
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

    ax.legend(fontsize=8.5, frameon=False)

    fpath = outdir / 'fig4_length_weight.png'
    fig.tight_layout()
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    return fpath


def fig_condition_index(df_ci: pd.DataFrame, outdir: Path) -> Path:
    """Fig 5: Body condition index distributions by class."""
    df_ci2 = df_ci.copy()
    df_ci2['Class'] = df_ci2['Sexo'].fillna('Unknown')
    label_map = {'M': 'Male', 'F': 'Female', 'J': 'Juvenile', 'Unknown': 'Undetermined (ref.)'}

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for cls in ['M', 'F', 'J', 'Unknown']:
        sub = df_ci2[df_ci2['Class'] == cls]['CI'].dropna()
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=15, color=PALETTE[cls], alpha=0.65,
                label=f'{label_map[cls]} (n\u202f=\u202f{len(sub)})',
                edgecolor='white', linewidth=0.3)

    ax.axvline(0, color='black', linestyle='--', linewidth=1,
               label='Reference (0)')
    ax.set_xlabel('Body condition index (log-residuals)', fontsize=11)
    ax.set_ylabel('Frequency', fontsize=11)
    ax.set_title(
        f'Body condition index \u2014 {SP_ITALIC}',
        fontsize=11
    )
    ax.legend(fontsize=9, frameon=False)

    fpath = outdir / 'fig5_condition_index.png'
    fig.tight_layout()
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    return fpath



# ── Main ─────────────────────────────────────────────────────────────────────

def run_all(data_path: str | Path, outdir: str | Path) -> dict:
    """
    Runs the complete morphometric analysis pipeline.
    Returns a dict with all results (suitable for Streamlit or CLI).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_morpho(data_path)

    desc = descriptive_stats(df)
    dimorphism = sexual_dimorphism_tests(df)
    gmm = fit_gmm(df, n_components_range=(1, 6))

    # Adults-only GMM for BIC comparison panel
    adults = df[df['Sexo'].isin(['M', 'F']) | df['Sexo'].isna()]
    gmm_adults = fit_gmm(adults, n_components_range=(1, 6))

    lw = length_weight_regression(df)
    df_ci = condition_index(df)
    cei_comp = cei_comparison(df)

    # Save tables
    desc.to_csv(outdir / 'table1_descriptive.csv', index=False)
    dimorphism.to_csv(outdir / 'table2_dimorphism.csv', index=False)

    # Figures
    f1 = fig_size_distribution(df, outdir)
    f2 = fig_gmm(gmm, df, outdir, adults_bic=gmm_adults['bic_scores'])
    f3 = fig_dimorphism(df, outdir)
    f4 = fig_length_weight(df, lw, outdir)
    f5 = fig_condition_index(df_ci, outdir)

    return {
        'df': df,
        'descriptive': desc,
        'dimorphism': dimorphism,
        'gmm': gmm,
        'length_weight': lw,
        'df_ci': df_ci,
        'cei_comparison': cei_comp,
        'figures': [f1, f2, f3, f4, f5],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Morphometric analysis of A. reverberii')
    parser.add_argument('--data', default='data/morphometrics_laguna_azul.csv')
    parser.add_argument('--outdir', default='outputs')
    args = parser.parse_args()

    print("Running morphometric analysis...")
    results = run_all(args.data, args.outdir)

    print("\n── Descriptive Statistics ──")
    print(results['descriptive'].to_string(index=False))

    print("\n── Sexual Dimorphism Tests ──")
    print(results['dimorphism'].to_string(index=False))

    print("\n── GMM: optimal k ──")
    print(f"  k = {results['gmm']['best_k']}")
    for i, (m, s, w) in enumerate(zip(
        results['gmm']['means'],
        results['gmm']['stds'],
        results['gmm']['weights']
    )):
        print(f"  Class {i+1}: μ={m:.2f} mm, σ={s:.2f} mm, weight={w:.3f}")

    print("\n── Length-Weight Regression ──")
    for grp, r in results['length_weight'].items():
        print(f"  {grp}: W = {r['a']} × LHU^{r['b']}  (r²={r['r2']}, n={r['n']})")

    print("\n── Comparison with Cei (1969) ──")
    cc = results['cei_comparison']
    print(f"  Males   present: n={cc['males_present']['n']}, mean LHU={cc['males_present']['mean']} mm")
    print(f"  Males   Cei1969: n={cc['cei_males']['n']}, mean SVL={cc['cei_males']['LHU_mean']} mm")
    print(f"  Females present: n={cc['females_present']['n']}, mean LHU={cc['females_present']['mean']} mm")
    print(f"  Females Cei1969: n={cc['cei_females']['n']}, mean SVL={cc['cei_females']['LHU_mean']} mm")

    print(f"\nFigures saved to: {Path(args.outdir).resolve()}")
