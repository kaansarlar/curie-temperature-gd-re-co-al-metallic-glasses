# Curie temperature in Gd-containing RE–Co–Al metallic glasses

Companion data and reproducible calculations for **Composition Based Estimation of Curie Temperature in Gd Containing RE Co Al Metallic Glasses**.

**Status:** working research manuscript, version 1. Not a published or accepted article. Numerical reproduction is complete; author declarations and a uniform primary-source eligibility audit remain to be finalized. No Zenodo deposit or DOI is provided.

## Scientific scope

The study compares an affine de Gennes–Co relation with a composition-dependent effective-scale relation for a restricted Gd-containing RE–Co–Al amorphous family. Phase eligibility must be established independently. Fe/Ni substitutions, RE-free systems, annealed/crystallized states and universal Tc prediction are not validated applications.

## Evidence layers

| Evaluation | MAE (K) | Maximum error (K) |
|---|---:|---:|
| Primary 15 records, affine model | 1.848 | 7.820 |
| Expanded 28 records, fixed RE–Co family | 1.747 | 5.951 |
| Expanded 28 records, complete family selection | 2.360 | 11.765 |
| Six new recipes from one source, frozen affine | 2.037 | 4.651 |
| Same six recipes, frozen RE–Co | 3.521 | 6.894 |

The 15- and 28-record scores are different datasets. The expanded model reuses earlier observations for development; historical holdout labels do not remain untouched tests. The six later recipes are outside the original joint composition support and from one publication; they are retrospective raw extrapolation checks. None was used for fitting or tuning. DOI/composition groups are not independently certified laboratories. The provisional author-dependency screen is not identity certification.

## Reproduce

Use Python 3.12 and the pinned requirements:

```bash
python -m pip install -r requirements.txt
python scripts/reproduce.py
python scripts/make_figures.py
```

The numerical audit independently reconstructs 99 archived outer models, their inner regularization/family choices, and frozen arithmetic. Results are written under `results/`. The figure script rebuilds the two original manuscript figures from archived held-out predictions. Source-article figures and PDFs are not redistributed.

## Files

- `manuscript/`: anonymous English draft in the supplied journal format and submission-status notes.
- `data/`: development fold predictions, six frozen checks, provenance flags and metrics.
- `models/`: frozen coefficients, outer fold models and inner decisions.
- `scripts/`: self-contained numerical reconstruction and figure generation.
- `results/`: reproducibility outputs.
- `figures/`: original PNG/SVG scientific figures.
- `references/`: IEEE reference mapping with original-source DOIs.
- `provenance/`: numerical logs, conditional author-dependency screen and runtime versions.

## Use and interpretation

The de Gennes factors are free-trivalent-ion proxies. They are not measured alloy moments or microscopic exchange constants. All elemental fractions refer to the whole alloy; Co25 means 0.25. Neither equation predicts glass formation, a calibrated uncertainty interval, or an annealing correction. The original study and frozen coefficients are retained even where the later check favors the simpler affine relation.

Code is MIT licensed. The original organization of curated data, documentation and figures follows the retained CC BY 4.0 terms; factual measurements must be attributed to their experimental sources. The manuscript draft and journal branding remain subject to the author's and journal's rights. Author declarations, source-method consistency and identity-resolved source-network sensitivity require review before journal submission.
