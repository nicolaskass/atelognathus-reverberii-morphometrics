# Dataset: Morphometric data of *Atelognathus reverberii* — Laguna Azul, 2019–2020

## File: `morphometrics_laguna_azul.csv`

First-capture morphometric measurements of 276 uniquely marked individuals
of *Atelognathus reverberii* (Cei, 1969) collected at Laguna Azul,
Meseta de Somuncurá, Río Negro Province, Argentina.

Field work was conducted across 14 sampling occasions in three primary
sessions: October 2019 (4 secondary occasions), November 2019 (5 occasions),
and February 2020 (5 occasions). Each individual was marked at first capture
with a unique visible implant elastomer (VIE) colour code applied to one or
more limbs. This dataset contains only first-capture events; recaptured
individuals are excluded to ensure independence of observations.

---

## Column descriptions

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `ID` | integer | — | Unique individual identifier (assigned at first capture) |
| `date` | date (YYYY-MM-DD) | — | Date of first capture |
| `SUL_mm` | float | mm | Snout–urostyle length: straight-line distance along the dorsal midline from the tip of the snout to the posterior end of the urostyle, with the animal in dorsal recumbency and hind limbs extended posteriorly. **Note:** this differs from the classical snout-vent length (SVL, measured to the cloaca) used by Cei (1969); SUL < SVL by a species-specific offset not yet calibrated. |
| `mass_g` | float | g | Body mass measured immediately after capture (digital balance, resolution 0.1 g) |
| `MW_mm` | float | mm | Mouth width: transverse distance between the commissures of the mouth, used as a proxy for head width and gape (digital calliper, resolution 0.1 mm) |
| `sex_class` | string | — | Sex/age class assigned in the field based on external characters (Cei 1969): `male` = nuptial pads present on pollex and second finger; `female` = nuptial pads absent with relatively broader body; `juvenile` = no secondary sexual characters or remnant larval tail present; `NaN` = could not be determined in the field |
| `stomach_content_present` | boolean | — | `True` if stomach lavage yielded prey items (non-zero Eppendorf sample), `False` if lavage was negative or not performed |

---

## Sample composition

| Class | n |
|-------|---|
| Male | 75 |
| Female | 47 |
| Juvenile | 25 |
| Undetermined | 129 |
| **Total** | **276** |

---

## Notes

- All measurements were taken by a single observer (N.A.K.), eliminating
  inter-observer error within this dataset. Future monitoring surveys using
  this protocol should include an inter-observer calibration phase before
  comparing condition indices across teams.
- Sex assignment was made using external morphological characters per
  Cei (1969). This method is unvalidated against internal examination for
  this species; see the companion paper (Kass et al., *Ichthyology &
  Herpetology*, in review) for a full discussion of classification
  uncertainty and sensitivity analysis.
- Intra-observer measurement error (CV) was <1.2% for SUL and <2.1%
  for MW, assessed from 20 re-measured individuals.

---

## Citation

Kass, N.A., Kass, C.A., Tettamanti, G., Kacoliris, F.P., and Williams, J.D.
Morphometric characterization and population structure of *Atelognathus
reverberii* (Cei, 1969) based on a large field sample from the Somuncurá
Plateau, Patagonia. *Ichthyology & Herpetology*, in review.

Data archived at Zenodo: [DOI TO BE ASSIGNED UPON ACCEPTANCE]
