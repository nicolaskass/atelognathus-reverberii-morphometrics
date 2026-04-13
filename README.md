# Morphometric analysis of *Atelognathus reverberii* — Kass et al.

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)


Reproducibility repository for:

> Kass, N.A., Kass, C.A., Tettamanti, G., Kacoliris, F.P., and Williams, J.D.
> **Morphometric characterization and population structure of
> *Atelognathus reverberii* (Cei, 1969) based on a large field sample
> from the Somuncurá Plateau, Patagonia.**
> *Ichthyology & Herpetology*, sent for review.

All figures, tables, and statistics reported in the paper are produced
by running `scripts/morfometria.py` on the dataset in `data/`.
Results are bit-for-bit identical to the published values.

---

## Repository structure

```
atelognathus-reverberii-morphometrics/
│
├── data/
│   ├── morphometrics_laguna_azul.csv   ← clean dataset (276 individuals)
│   └── README_data.md                  ← full column descriptions
│
├── scripts/
│   ├── data_loader.py                  ← data loading and cleaning
│   └── morfometria.py                  ← all analyses and figures
│
├── figures/                            ← publication figures (300 DPI PNG)
│   ├── fig1_size_distribution.png
│   ├── fig2_gmm.png
│   ├── fig3_dimorphism.png
│   ├── fig4_length_weight.png
│   └── fig5_condition_index.png
│
├── outputs/                            ← generated tables (CSV)
│
├── requirements.txt
└── README.md
```

---

## Installation

### Requirements

- Python 3.9 or later
- pip

### Step-by-step setup

**1. Clone the repository**

```bash
git clone https://github.com/nicolaskass/atelognathus-reverberii-morphometrics
cd atelognathus-reverberii-morphometrics
```

**2. Create a virtual environment** (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate        # Linux / macOS
# .venv\Scripts\activate         # Windows
```

**3. Install dependencies**

```bash
pip install -r requirements.txt
```

Dependencies and why each is needed:

| Package | Version | Purpose |
|---------|---------|---------|
| `pandas` | ≥2.0 | Data loading, cleaning, and tabular operations |
| `numpy` | ≥1.24 | Numerical computations, bootstrap resampling |
| `matplotlib` | ≥3.7 | All publication figures |
| `scipy` | ≥1.11 | Statistical tests (Shapiro–Wilk, Mann–Whitney U, Kruskal–Wallis, linear regression) |
| `scikit-learn` | ≥1.3 | Gaussian Mixture Models for size-class analysis |
| `openpyxl` | ≥3.1 | Reading Excel files (only needed if using the raw `.xlsx` input) |

---

## Running the analysis

### Full pipeline (all figures and tables)

```bash
python scripts/morfometria.py \
    --data data/morphometrics_laguna_azul.csv \
    --outdir outputs/
```

This produces in `outputs/`:
- `table1_descriptive.csv` — Table 1 of the paper
- `table2_dimorphism.csv` — Table 2 of the paper
- `fig1_size_distribution.png` through `fig5_condition_index.png`

And prints a full summary to stdout including all statistics reported in the paper.

### Parameters explained

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--data` | `data/morphometrics_laguna_azul.csv` | Path to the input dataset. Must be a CSV with columns `ID`, `SUL_mm`, `mass_g`, `MW_mm`, `sex_class`. See `data/README_data.md` for full column specifications. |
| `--outdir` | `outputs/` | Directory where figures and tables are saved. Created automatically if it does not exist. |

### Running individual analyses

You can also import the analysis functions directly from Python:

```python
import sys
sys.path.insert(0, 'scripts')
from data_loader import load_morpho
from morfometria import (
    descriptive_stats,
    sexual_dimorphism_tests,
    fit_gmm,
    length_weight_regression,
    condition_index,
)

# Load data
df = load_morpho('data/morphometrics_laguna_azul.csv')

# Descriptive statistics with 95% bootstrap CI
desc = descriptive_stats(df)

# Sexual dimorphism: Mann-Whitney U + Cohen's d
dim = sexual_dimorphism_tests(df)

# Gaussian Mixture Models (BIC selection, k=1..5)
gmm = fit_gmm(df, n_components_range=(1, 6))
print(f"Optimal k: {gmm['best_k']}")  # → 1 (unimodal)

# Length-weight regression (log-log OLS)
lw = length_weight_regression(df)

# Body condition index (residuals from reference regression)
df_ci = condition_index(df)
```

---

## Analysis details

### Descriptive statistics

Means, standard deviations, and ranges per sex/age class.
Bootstrap 95% confidence intervals for the mean are computed with
2,000 resampling iterations (fixed seed 42) using `numpy.random.default_rng`.

**Why bootstrap CIs?** Standard normal-theory CIs assume normality.
With non-normal variables (body mass was significantly non-normal in both
sexes; Shapiro–Wilk p < 0.006) and unequal group sizes, bootstrap CIs
are more reliable and require no distributional assumptions.

### Sexual dimorphism (Mann–Whitney U test)

The two-sided Mann–Whitney U test is used rather than a t-test because
body mass was significantly non-normal in both sexes. To maintain a
consistent framework across all three variables (SUL, mass, MW), the
same non-parametric test is applied to all.

Effect size is Cohen's d (difference of means / pooled SD), interpreted as:
- |d| < 0.2 : negligible
- 0.2–0.5   : small
- 0.5–0.8   : medium
- ≥ 0.8     : large

**Power analysis:** Computed via the normal-theory asymptotic relative
efficiency (ARE) correction for the Mann–Whitney U test (Lehmann & Romano,
2006: ARE = 3/π ≈ 0.955 for normal data). Effective sample sizes are
n_eff = n × ARE. Power for d = 0.5: 74.7% (marginally below the 80%
conventional threshold; power exceeds 80% for d ≥ 0.53). Power for
d = 0.8: 98.8%.

### Gaussian Mixture Models (size-class structure)

GMMs are fitted to the SUL distribution using
`sklearn.mixture.GaussianMixture` with 10 random initialisations per model.
The optimal number of components k is selected by minimising BIC over
k ∈ {1, 2, 3, 4, 5}.

**Why start at k = 1?** Starting the BIC search at k = 1 allows the null
hypothesis of unimodality to be evaluated directly. Starting at k = 2 (a
common mistake) presupposes at least two components and biases model
selection toward spurious multimodality.

**Why not use the standard LRT with χ² tables?** Under H₀ (k = 1), the
mixing weight parameter of the alternative (k = 2) is unidentified,
violating the regularity conditions for the asymptotic χ² approximation.
This is a well-known irregular testing problem (McLachlan, 1987,
*Applied Statistics*). We use a bootstrapped LRT instead: 999 datasets
are simulated under the fitted k = 1 model and Λ = max(0, 2[ℓ₂ − ℓ₁])
is computed for each. The bootstrapped p-value is the proportion of
simulated Λ values ≥ the observed value.

**Result:** BIC selects k = 1 for both the full sample
(ΔBIC k=2 vs k=1 = +18.3) and the adult-only subset (ΔBIC = +14.9).
Bootstrapped LRT: Λ = 0, p = 1.00 (B = 999). The SUL distribution
is unimodal. The k = 2 solution is shown in Figure 2 (right panel)
for descriptive comparison only.

### Length–weight relationship

Model: log(W) = log(a) + b · log(SUL), fitted by OLS on log-transformed data.
The exponent b is tested against:
- b = 2 (two-dimensional isometry, expected for dorsoventrally flattened species)
- b = 3 (three-dimensional isometry, the standard expectation for compact bodies)

using a one-sample t-test on the estimated slope with SE(b) from the regression.

**Why two reference exponents?** *A. reverberii* has a dorsoventrally
flattened body form consistent with sheltering under flat rocks and in shallow
clay lagoons. Two-dimensional isometry (b ≈ 2) predicts that mass scales
with the square of linear dimensions, as would be expected for an animal
that grows primarily in the plane of the substrate rather than volumetrically.

### Body condition index

Condition = residuals from log(mass) ~ log(SUL) regression.

**Why use the undetermined-sex adults as the reference group?**
Fitting the baseline regression to all individuals and then comparing
residuals among those same individuals is circular: it forces the
overall mean residual to zero and can bias group comparisons. Using
the undetermined-sex adults (n = 129) as a sex-neutral reference
avoids this circularity (Green, 2001, *Ecology*).

Differences among classes are tested with Kruskal–Wallis H (df = 2),
followed by pairwise Mann–Whitney U with Bonferroni correction
(3 pairs; family-wise α = 0.05).

---

## Reproducing exact paper values

All random operations use fixed seeds. Running `morfometria.py` on the
provided dataset will reproduce every value in the paper exactly:

| Statistic | Paper | Script output |
|-----------|-------|---------------|
| Male SUL mean | 28.76 mm | ✓ |
| Female SUL mean | 28.64 mm | ✓ |
| Juvenile SUL mean | 24.57 mm | ✓ |
| Mann-Whitney U (SUL) | 1917 | ✓ |
| GMM optimal k | 1 | ✓ |
| BIC k=1 | 1396.97 | ✓ |
| ΔBIC (k=2 vs k=1) | +18.3 | ✓ |
| L-W exponent b (all) | 1.969 | ✓ |
| L-W exponent b (juvenile) | 2.530 | ✓ |
| Kruskal-Wallis H (condition) | 6.25 | ✓ |

---

## Using your own data

`morfometria.py` accepts any CSV with the following columns:

| Required column | Description |
|-----------------|-------------|
| `SUL_mm` | Body length in mm (snout to urostyle or vent) |
| `mass_g` | Body mass in g |

Optional columns that unlock additional analyses:

| Optional column | Description |
|-----------------|-------------|
| `MW_mm` | Mouth width in mm |
| `sex_class` | `male`, `female`, `juvenile`, or blank/NaN |
| `ID` | Individual identifier |
| `date` | Capture date |

If your length variable has a different column name (e.g., `SVL_mm`),
rename it before loading or edit the `MEASURE_COL` variable at the
top of `data_loader.py`.

---

## Citing this work

If you use this code or dataset, please cite:

```bibtex
@article{kass2026,
  author  = {Kass, Nicol{\'a}s Ariel and Kass, Camila Alejandra and
             Tettamanti, Germ{\'a}n and Kacoliris, Federico Pablo and
             Williams, Jorge Daniel},
  title   = {Morphometric characterization and population structure of
             \textit{Atelognathus reverberii} ({Cei}, 1969) based on a
             large field sample from the Somuncur{\'a} Plateau, {Patagonia}},
  journal = {Ichthyology \& Herpetology},
  year    = {2026},
  note    = {Sent}
}
```

---

## Contact

Nicolás Ariel Kass · nkass@fcnym.unlp.edu.ar  
Sección Herpetología, División Zoología Vertebrados  
Facultad de Ciencias Naturales y Museo, UNLP · La Plata, Argentina
