# Datasets — *Atelognathus reverberii* at Laguna Azul, 2019–2020

Field data collected by N.A.K. across 14 sampling occasions in three
primary sessions: October 2019 (oct_08 through oct_11), November 2019
(nov_05, nov_24, nov_26, nov_27, nov_28), and February 2020 (feb_01
through feb_05).

Each individual was marked at first capture with a unique visible implant
elastomer (VIE) colour code applied to one or more of the four limbs.
Three files are provided:

| File | Rows | What it contains |
|---|---|---|
| `morphometrics_laguna_azul.csv` | 274 | First-capture morphometric measurements; the analytic dataset for the paper |
| `capture_history.csv` | 275 | Wide binary matrix (individual × occasion); detection events per marked individual |
| `recaptures_with_measurements.csv` | 27 | Subset of recapture events with complete morphometric measurements |

The 275 individuals in `capture_history.csv` are the full marked
cohort; 274 of those yielded complete first-capture measurements (mass,
SUL, MW all recorded) and form the analytic sample in
`morphometrics_laguna_azul.csv`. The one excluded individual had
incomplete measurements at first capture (e.g., missing mass or MW
during field processing).

---

## File: `morphometrics_laguna_azul.csv`

First-capture morphometric measurements (n = 274 of 275 marked
individuals).

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `individual_id` | integer | — | Unique individual identifier (assigned at first capture) |
| `date` | string (DD-MM-YYYY) | — | Date of first capture |
| `body_mass_g` | float | g | Body mass measured immediately after capture (digital balance, resolution 0.1 g) |
| `SUL_mm` | float | mm | Snout–urostyle length: straight-line distance along the dorsal midline from the tip of the snout to the posterior end of the urostyle, with the animal in dorsal recumbency and hind limbs extended posteriorly. **Note:** this differs from the classical snout-vent length (SVL, measured to the cloaca) used by Cei (1969); SUL < SVL by a species-specific offset not yet calibrated. |
| `mouth_width_mm` | float | mm | Mouth width: transverse distance between the commissures of the mouth, used as a proxy for head width and gape (digital calliper, resolution 0.1 mm). Not validated against direct cranial measurements in Batrachylidae; reported as a body-form descriptor. |
| `sex_class` | string | — | Sex/age class assigned in the field by external characters (Cei 1969): `macho` (male) = nuptial pads present on pollex and second finger; `hembra` (female) = no nuptial pads, with relatively broader body; `juvenil` (juvenile) = no clear secondary sexual character or remnant larval tail present. Empty cells, `?`, or `Foto Parásito` are treated as undetermined by the analysis script |
| `stomach_content_present` | boolean | — | `True` if stomach lavage yielded prey items, `False` otherwise |
| `marca` | string | — | Free-text VIE mark description from the field datasheet |
| `mark_MAI`, `mark_MAD`, `mark_MPI`, `mark_MPD` | string | — | VIE colour code per limb (MAI = left anterior, MAD = right anterior, MPI = left posterior, MPD = right posterior); `Y` = present, blank = no mark on that limb |
| `notes` | string | — | Free-text field notes |

### Sample composition (after sex-label normalisation)

| Class | n | Note |
|---|---|---|
| Male | 75 | `macho` (also `Macho`, `MACHO`) |
| Female | 47 | `hembra` (also `Hembra`, `HEMBRA`, `hemba`) |
| Juvenile | 25 | `juvenil` (also `Juvenil`, `juvuenil`) |
| Undetermined | 127 | NaN, `?`, or `Foto Parásito` |
| **Total** | **274** | |

The analysis script in `../scripts/morfometria.py` performs the
case-insensitive label normalisation at load time.

---

## File: `capture_history.csv`

Wide binary capture-history matrix for the full marked cohort
(*n* = 275 individuals × 14 sampling occasions).

| Column | Type | Description |
|---|---|---|
| `individual_id` | integer | Unique identifier (matches `morphometrics_laguna_azul.csv` where the individual yielded complete first-capture measurements) |
| `oct_08`, `oct_09`, `oct_10`, `oct_11` | binary | October 2019 secondary occasions |
| `nov_05`, `nov_24`, `nov_26`, `nov_27`, `nov_28` | binary | November 2019 secondary occasions |
| `feb_01`, `feb_02`, `feb_03`, `feb_04`, `feb_05` | binary | February 2020 secondary occasions |

Cell value: 1 if the individual was detected on that occasion, 0 otherwise.
Row sum ≥ 2 indicates the individual was recaptured at least once. Across
the 339 detection events in this matrix, individuals were captured one,
two or three times: 220 individuals were captured once, 46 twice, and 9
three times, so 55 of 275 individuals (20.0%) were recaptured at least
once.

This file is the input to the companion capture–recapture analysis
(Kass et al., in preparation). The present paper uses it only to derive
the marked-cohort counts cited in the Methods.

---

## File: `recaptures_with_measurements.csv`

Subset of the recapture events for which complete morphometric
measurements were taken (n = 27 events). Not used in the present paper
(which restricts analysis to first-capture events to ensure
independence). Included here for completeness and for the companion
capture–recapture analysis.

Columns are the same as `morphometrics_laguna_azul.csv` minus the
`mark_MAI`/`mark_MAD`/`mark_MPI`/`mark_MPD` per-limb columns; the
`marca` text column carries the matched VIE pattern.

---

## Notes on data provenance

- All measurements were taken by a single observer (N.A.K.), eliminating
  inter-observer error within this dataset. Future monitoring surveys
  using this protocol should include an inter-observer calibration phase
  before comparing condition indices across teams.
- Intra-observer measurement error (CV) was < 1.2% for SUL and < 2.1%
  for MW, assessed from 20 re-measured individuals at first capture.
- Sex assignment used external morphological characters per Cei (1969).
  This method is unvalidated against internal examination for this
  species; the paper reports a Buonaccorsi (2010) linear-correction
  sensitivity analysis (including an asymmetric scenario) to bound the
  effect of plausible field-classification error on the dimorphism test.
- Field work was conducted under permit issued by the Ministerio de
  Ambiente y Desarrollo Sostenible, Río Negro Province, Argentina.

---

## Citation

Kass, N.A., Kass, C.A., Tettamanti, G., Kacoliris, F.P., and Williams,
J.D. (2026). Field-based morphometric characterization and size
distribution of a Laguna Azul population of *Atelognathus reverberii*
(Cei, 1969), Somuncurá Plateau, Patagonia. *Ichthyology & Herpetology*,
in review.

Data archived at Zenodo: [DOI to be assigned upon acceptance].
