# PyGeoMorph v2.0
**Geometric Shell Morphometry for Marine Invertebrates**

[![Build](https://img.shields.io/badge/build-GitHub_Actions-0a1822)](.github/workflows/build.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-143247)](https://www.python.org/)
[![Platforms](https://img.shields.io/badge/platforms-Win%20%7C%20Linux%20%7C%20macOS--Intel%20%7C%20macOS--ARM-d6553d)](#)

> 🇮🇹 *Versione italiana disponibile in [README.it.md](README.it.md).*

> **Author**: Francesco Paolo Patti  
> **Affiliation**: Stazione Zoologica Anton Dohrn, Naples, Italy  
> **Contact**: francesco.patti@szn.it  
> **Year**: 2026

---

PyGeoMorph is an all-in-one Python suite for **geometric morphometric analysis**
of marine invertebrates (gastropods, bivalves, brachiopods, echinoids). The
application replaces, in a single open-source tool, the classical workflow
**tpsDig2 + tpsRelw + SPSS + SYSTAT** widely used in malacological studies —
including Criscione & Patti (2010) on Mediterranean *Rissoa* species, which
provides the reference template of the software and the conceptual starting
point of this 2026 modernisation.

---

## 1. Scientific framework

### 1.1 Continuity with Criscione & Patti (2010)

The original workflow of **Criscione, F. & Patti, F.P. (2010)** "*Similar shells
are not necessarily a reliable guide to phylogeny: Rissoa guerinii Recluz, 1843,
and Rissoa lia (Monterosato, 1884) (Caenogastropoda: Rissoidae): a case study*",
*The Nautilus* **124**(3): 117–128, was based on:

- 19 landmarks placed manually with tpsDig2 on the shell profile (apex, three
  last whorls, aperture, right and left profile);
- centroid size (CS), uniform components U1/U2, and relative warps (RW1...RW29)
  computed with tpsRelw v.1.45;
- ANOVA / ANCOVA, discriminant analysis, and TPS plotting in SPSS v.15 and
  SYSTAT v.12.

This workflow is fully preserved in PyGeoMorph as the built-in **"rissoa"
template** (19 landmarks, links and slider definitions for the semi-landmarks
reproduced from the 2010 paper), so that anyone publishing along that
methodology can directly reproduce previous results and compare them on new
samples.

### 1.2 Innovations over 2010

In the last 15 years geometric morphometrics has changed deeply. PyGeoMorph
updates the 2010 workflow to current standards (geomorph 4.x R-package, June
2024 / March 2026; Mitteroecker & Schaefer 2022; Adams, Collyer &
Kaliontzopoulou 2024):

| Aspect                   | 2010 pipeline (Criscione & Patti)        | PyGeoMorph v2.0 (2026)                                       |
|--------------------------|------------------------------------------|--------------------------------------------------------------|
| Software                 | tpsDig2 + tpsRelw + SPSS + SYSTAT        | Single all-in-one Python application                         |
| Landmark placement       | Manual, subjective                       | Automatic (curvature + active contours) + manual correction  |
| Landmark type            | Fixed only                               | **Fixed + sliding semi-landmarks** (BE-min sliding)          |
| Superimposition          | GPA on fixed LM                          | GPA with sliding semi-landmarks along shell whorls           |
| Shape variables          | Uniform + Relative Warps (RW1...RWn)     | Tangent-space PCA shape variables                            |
| Statistics               | Parametric ANOVA / ANCOVA                | **Procrustes ANOVA via RRPP** (999+ permutations)            |
| Group comparisons        | Classical discriminant analysis          | RRPP + pairwise Mahalanobis + MANOVA                         |
| Allometry                | Step-wise regression                     | Multivariate regression + permutation + **CAC** + **HOS**    |
| Outline analysis         | (not considered)                         | **Elliptic Fourier Analysis** (EFA), complementary           |
| Phylogeny                | (separate analysis)                      | **Kmult** phylogenetic signal, integrated                    |
| Modularity / Integration | (not available)                          | **Covariance Ratio test** (Adams 2016)                       |
| Export                   | .tps file, Excel sheet                   | Multi-block CSV + Rohlf TPS + JSON                           |
| Executables              | (legacy Windows-only)                    | **Win + Linux + macOS Intel + macOS Apple Silicon**          |

### 1.3 Methods implemented (with references)

- **Generalized Procrustes Analysis** with unit-centroid-size constraint on the
  consensus (Gower 1975; Goodall 1991; Rohlf & Slice 1990).
- **Sliding semi-landmarks** via chord-tangent BE-min sliding (Bookstein 1997;
  Gunz, Mitteroecker & Bookstein 2005).
- **Procrustes ANOVA via RRPP** (Randomization of Residuals in a Permutation
  Procedure; Collyer & Adams 2018; Adams & Collyer 2018, *Methods Ecol Evol*).
- **Common Allometric Component** (CAC; Mitteroecker, Gunz, Bernhard, Schaefer
  & Bookstein 2004, *J Hum Evol* 46:679–698).
- **Homogeneity of Slopes** test for group-specific allometric trajectories
  (Adams & Collyer 2009).
- **Pairwise Mahalanobis** distances with pooled within-group covariance.
- **Wilks' Λ MANOVA** with Bartlett's approximation.
- **Disparity** as Procrustes MSE, global and per group (Foote 1993; Zelditch
  et al. 2012).
- **Thin-Plate Spline** with kernel U(r) = r² log r and bending energy
  computation (Bookstein 1989, *IEEE PAMI* 11:567–585).
- **Normalised Elliptic Fourier Analysis** (Kuhl & Giardina 1982; Bonhomme,
  Picq, Gaucherel & Claude 2014, *Momocs*, *J Stat Softw* 56).
- **Kmult phylogenetic signal** (Adams 2014, *Syst Biol* 63:685–697) — optional,
  requires a user-supplied phylogenetic distance matrix.
- **Covariance Ratio test** for modularity / integration (Adams 2016, *Methods
  Ecol Evol* 7:565–572).
- **Combined automatic landmark detection**: equidistant outline + curvature
  peaks (deep-learning-style suggestion), to be used as an initial bootstrap
  before manual correction.

---

## 2. Software features

- **Image input**: JPG, PNG, TIF, PDF (via poppler / pdf2image), multi-file
  drag & drop.
- **Built-in templates**: Rissoidae (19 LM, after Criscione & Patti 2010),
  bivalve outline (15 LM), echinoid test (12 LM), generic outline (15 sliding
  LM).
- **Automatic landmark detection** in 3 modes: curvature peaks (default),
  equidistant outline, Shi-Tomasi / Harris corners.
- **Landmark editor**: add / move / delete, zoom, overlay numbered / dots /
  template links / EFA outline / slider highlighting.
- **Full statistical pipeline** triggered with one click: GPA → PCA → TPS →
  EFA → RRPP → MANOVA → Mahalanobis → Allometry + CAC + HOS → Disparity →
  (Kmult).
- **All statistics visible**: every result is shown both as a numeric card
  (Wilks Λ, F, Z, p, R²) and as an interactive Plotly figure.
- **Automatic export**:
  - **CSV** multi-block with raw landmarks, aligned shape coordinates, PC
    scores, eigenvalues, EFA scores, and Procrustes ANOVA table.
  - **TPS** in Rohlf tpsDig2 format (compatible with other morphometric
    software).
- **Cross-platform**: native executables for Windows, Linux, macOS Intel, and
  **macOS Apple Silicon**.

---

## 3. Installation (development, from source)

```bash
git clone <repo>
cd pygeomorph
python -m venv venv
source venv/bin/activate           # Windows: venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Open the browser at **http://127.0.0.1:5000**.

### System dependencies for PDF support

PyGeoMorph uses `pdf2image` to read `.pdf` files as photographic images. It
requires `poppler` installed at system level:

| OS              | Command                                                          |
|-----------------|------------------------------------------------------------------|
| Ubuntu / Debian | `sudo apt install poppler-utils`                                 |
| macOS           | `brew install poppler`                                           |
| Windows         | download binaries from [oschwartz10612/poppler-windows](https://github.com/oschwartz10612/poppler-windows/releases/) and add `bin/` to `PATH` |

---

## 4. Automatic executables via GitHub Actions

The `.github/workflows/build.yml` workflow automatically produces four native
binaries on every `v*` tag push:

- `PyGeoMorph-windows-x64.exe`
- `PyGeoMorph-linux-x64`
- `PyGeoMorph-macos-x86_64` (Mac Intel)
- `PyGeoMorph-macos-arm64` (Mac Apple Silicon: M1, M2, M3, M4)

### 4.1 Why two separate macOS binaries?

Building a macOS universal2 executable (Intel + ARM in a single file) with
PyInstaller fails when any binary wheel in the dependency tree is single-arch.
NumPy, SciPy, OpenCV, and Pillow on PyPI are normally single-arch: the
universal2 bundling fails arch validation with errors like `Bad CPU type in
executable` or `IncompatibleBinaryArchError`.

The **standard 2024-2026 solution** is no longer to attempt a universal2 build,
but to **build separately** on the two GitHub runners:

- `macos-13` → Intel x86_64 runner → produces `PyGeoMorph-macos-x86_64`
- `macos-14` → Apple Silicon arm64 runner (the new `macos-latest` default since
  late 2024) → produces `PyGeoMorph-macos-arm64`

This is the same pattern adopted by scientific Python projects such as *napari*
and *scikit-image-tools*.

### 4.2 To publish a release

```bash
git tag v2.0.0
git push origin v2.0.0
```

GitHub Actions will build the four executables in ~10 minutes and publish them
automatically as assets of a new GitHub Release.

---

## 5. Typical workflow

1. **Sidebar → Acquisition**: choose a template (e.g. *Rissoidae 19 LM*),
   group / taxon, number of landmarks, detection method.
2. **Upload images** (drag & drop, multi-file): for each specimen the software
   automatically extracts the outline and suggests N landmarks at curvature
   peaks.
3. **Digitize tab**: review and, if needed, correct the landmarks (Add / Move
   / Delete).
4. **Sidebar → Pipeline**:
   - set the number of RRPP permutations (default 999),
   - keep **Use sliding semi-landmarks** active if the template defines sliders,
   - enable **Compute Elliptic Fourier** for complementary outline analysis,
   - click **Run Full Pipeline**.
5. Inspect the results in the tabs:
   - **Procrustes (GPA)** — aligned configurations + consensus + per-specimen
     Procrustes distance plot;
   - **PCA / Shape Space** — interactive morphospace + scree + PC scores
     table;
   - **TPS Deformation** — deformation grid consensus → PC1 extreme, with
     bending energy;
   - **EFA Outline** — morphospace based on elliptic Fourier descriptors;
   - **Statistics · RRPP** — Procrustes ANOVA, MANOVA, Mahalanobis heatmap,
     disparity, Kmult (if a phylogenetic distance matrix has been loaded);
   - **Allometry · CAC** — R², permutation p-value, HOS, CAC vs log(CS)
     trajectory;
   - **Raw Data** — full registry + raw coordinate table.
6. **Export CSV** or **Export TPS** from the sidebar for archival and analysis
   in external software.

---

## 6. Citation

If you use PyGeoMorph in a publication, please cite:

> Patti, F.P. (2026). *PyGeoMorph: a modern all-in-one Python pipeline for
> geometric morphometric analysis of marine invertebrates.*
> Stazione Zoologica Anton Dohrn, Naples, Italy.
> https://doi.org/10.5281/zenodo.20117484

…and the original paper that provides the conceptual basis for the
malacological template:

> Criscione, F., & Patti, F.P. (2010). Similar shells are not necessarily a
> reliable guide to phylogeny: *Rissoa guerinii* Recluz, 1843, and *Rissoa lia*
> (Monterosato, 1884) (Caenogastropoda: Rissoidae): a case study.
> *The Nautilus*, **124**(3): 117–128.

---

## 7. References

- Adams, D.C. (2014). A generalized K statistic for estimating phylogenetic
  signal from shape and other high-dimensional multivariate data.
  *Syst Biol* 63: 685–697.
- Adams, D.C. (2016). Evaluating modularity in morphometric data: challenges
  with the RV coefficient and a new test measure.
  *Methods Ecol Evol* 7: 565–572.
- Adams, D.C., Collyer, M.L., & Kaliontzopoulou, A. (2024). *Geomorph: Software
  for geometric morphometric analyses*. R package version 4.0.8.
  https://CRAN.R-project.org/package=geomorph
- Bonhomme, V., Picq, S., Gaucherel, C., & Claude, J. (2014). Momocs: Outline
  analysis using R. *J Stat Softw* 56: 1–24.
- Bookstein, F.L. (1989). Principal warps: thin-plate splines and the
  decomposition of deformations. *IEEE PAMI* 11: 567–585.
- Bookstein, F.L. (1997). Landmark methods for forms without landmarks:
  morphometrics of group differences in outline shape.
  *Med Image Anal* 1: 225–243.
- Collyer, M.L., & Adams, D.C. (2018). RRPP: An R package for fitting linear
  models to high-dimensional data using residual randomization.
  *Methods Ecol Evol* 9: 1772–1779.
- Criscione, F., & Patti, F.P. (2010). Similar shells are not necessarily a
  reliable guide to phylogeny: *Rissoa guerinii* Recluz, 1843, and *Rissoa lia*
  (Monterosato, 1884) (Caenogastropoda: Rissoidae): a case study.
  *The Nautilus* 124: 117–128.
- Gunz, P., Mitteroecker, P., & Bookstein, F.L. (2005). Semilandmarks in three
  dimensions. In: D.E. Slice (ed.), *Modern Morphometrics in Physical
  Anthropology*, Kluwer.
- Kuhl, F.P., & Giardina, C.R. (1982). Elliptic Fourier features of a closed
  contour. *Comput Graph Image Proc* 18: 236–258.
- Mitteroecker, P., Gunz, P., Bernhard, M., Schaefer, K., & Bookstein, F.L.
  (2004). Comparison of cranial ontogenetic trajectories among great apes and
  humans. *J Hum Evol* 46: 679–698.
- Mitteroecker, P., & Schaefer, K. (2022). Thirty years of geometric
  morphometrics: Achievements, challenges, and the ongoing quest for biological
  meaningfulness. *Am J Biol Anthropol* 178: 181–210.
- Rohlf, F.J., & Slice, D.E. (1990). Extensions of the Procrustes method for
  the optimal superimposition of landmarks. *Syst Zool* 39: 40–59.

---

## 8. License

MIT

---

## 9. Contact

For bug reports, suggestions, collaborations, or requests for new templates
specific to other taxa of marine invertebrates:

**Francesco Paolo Patti**  
Stazione Zoologica Anton Dohrn, Naples, Italy  
✉ francesco.patti@szn.it
