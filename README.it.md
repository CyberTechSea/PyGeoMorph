# PyGeoMorph v2.0
**Geometric Shell Morphometry for Marine Invertebrates**

[![Build](https://img.shields.io/badge/build-GitHub_Actions-0a1822)](.github/workflows/build.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-143247)](https://www.python.org/)
[![Platforms](https://img.shields.io/badge/platforms-Win%20%7C%20Linux%20%7C%20macOS--Intel%20%7C%20macOS--ARM-d6553d)](#)

> 🇬🇧 *English version available at [README.md](README.md).*

> **Author**: Francesco Paolo Patti  
> **Affiliation**: Stazione Zoologica Anton Dohrn, Naples, Italy  
> **Contact**: francesco.patti@szn.it  
> **Year**: 2026

---

PyGeoMorph è una suite all-in-one in Python per la **morfometria geometrica** di
invertebrati marini (gasteropodi, bivalvi, brachiopodi, echinoidi). L'applicazione
sostituisce in un unico software open-source il workflow storico
**tpsDig2 + tpsRelw + SPSS + SYSTAT** utilizzato in molti studi malacologici — incluso
quello di Criscione & Patti (2010) sui *Rissoa* del Mediterraneo, che costituisce
il template di riferimento del software e il punto di partenza concettuale di
questa modernizzazione 2026.

---

## 1. Inquadramento scientifico

### 1.1 Continuità con Criscione & Patti (2010)

Il workflow originale di **Criscione, F. & Patti, F.P. (2010)** "*Similar shells are
not necessarily a reliable guide to phylogeny: Rissoa guerinii Recluz, 1843, and
Rissoa lia (Monterosato, 1884) (Caenogastropoda: Rissoidae): a case study*",
*The Nautilus* **124**(3): 117–128, basato su:

- 19 landmark posizionati manualmente con tpsDig2 sul profilo della conchiglia
  (apex, tre ultime spire, apertura, profilo destro/sinistro);
- calcolo di centroid size (CS), uniform components U1/U2 e relative warps (RW1...RW29)
  con tpsRelw v.1.45;
- ANOVA / ANCOVA, discriminant analysis e TPS plotting in SPSS v.15 e SYSTAT v.12;

…è interamente preservato in PyGeoMorph come **template "rissoa"** built-in
(19 landmark, links e definizione di slider dei semi-landmark riprodotti dal lavoro
del 2010), in modo che chi pubblica seguendo quella metodologia possa riprodurre
direttamente i risultati pregressi e confrontarli su nuovi campioni.

### 1.2 Innovazioni rispetto al 2010

Negli ultimi 15 anni la morfometria geometrica è cambiata profondamente. PyGeoMorph
aggiorna il workflow del 2010 ai metodi attualmente standard (geomorph 4.x R-package,
giugno 2024 / marzo 2026; Mitteroecker & Schaefer 2022; Adams, Collyer &
Kaliontzopoulou 2024):

| Aspetto                  | Pipeline 2010 (Criscione & Patti)        | PyGeoMorph v2.0 (2026)                                     |
|--------------------------|------------------------------------------|------------------------------------------------------------|
| Software                 | tpsDig2 + tpsRelw + SPSS + SYSTAT        | Singola applicazione Python all-in-one                     |
| Posizionamento landmark  | Manuale, soggettivo                      | Automatico (curvature + active contours) + correzione manuale |
| Tipo di landmark         | Solo fissi                               | **Fissi + sliding semi-landmark** (BE-min sliding)         |
| Allineamento             | GPA su LM fissi                          | GPA con sliding semi-landmark sui contorni della conchiglia |
| Variabili di forma       | Uniform + Relative Warps (RW1...RWn)     | Tangent-space PCA shape variables                          |
| Statistica               | ANOVA / ANCOVA parametriche              | **Procrustes ANOVA con RRPP** (999+ permutazioni)          |
| Confronto tra gruppi     | Discriminant analysis classica           | RRPP + Mahalanobis pairwise + MANOVA                       |
| Allometria               | Regressione step-wise                    | Multivariate regression + permutazione + **CAC** + **HOS** |
| Outline                  | (non considerato)                        | **Elliptic Fourier Analysis** (EFA) complementare          |
| Filogenesi               | (analisi separata)                       | **Kmult** phylogenetic signal integrato                    |
| Modularity / Integration | (non disponibile)                        | **Covariance Ratio test** (Adams 2016)                     |
| Esportazione             | File .tps, foglio Excel                  | CSV multi-blocco + TPS Rohlf + JSON                        |
| Eseguibili               | (solo desktop Windows storici)           | **Win + Linux + macOS Intel + macOS Apple Silicon**        |

### 1.3 Metodi implementati (con riferimenti)

- **Generalized Procrustes Analysis** con vincolo di centroid size unitario sul
  consensus (Gower 1975; Goodall 1991; Rohlf & Slice 1990).
- **Sliding semi-landmark** mediante chord-tangent BE-min sliding
  (Bookstein 1997; Gunz, Mitteroecker & Bookstein 2005).
- **Procrustes ANOVA RRPP** (Randomization of Residuals in a Permutation Procedure;
  Collyer & Adams 2018; Adams & Collyer 2018, *Methods Ecol Evol*).
- **Common Allometric Component** (CAC; Mitteroecker, Gunz, Bernhard, Schaefer &
  Bookstein 2004, *J Hum Evol* 46:679–698).
- **Homogeneity of Slopes** test per traiettorie allometriche specifiche
  (Adams & Collyer 2009).
- **Mahalanobis pairwise** distances con covarianza pooled fra gruppi.
- **Wilks' Λ MANOVA** con approssimazione di Bartlett.
- **Disparity** come MSE Procrustes globale e per gruppo (Foote 1993; Zelditch et al. 2012).
- **Thin-Plate Spline** con kernel U(r) = r² log r e calcolo della bending energy
  (Bookstein 1989, *IEEE PAMI* 11:567–585).
- **Elliptic Fourier Analysis** normalizzata (Kuhl & Giardina 1982; Bonhomme, Picq,
  Gaucherel & Claude 2014, *Momocs*, *J Stat Softw* 56).
- **Kmult phylogenetic signal** (Adams 2014, *Syst Biol* 63:685–697) — opzionale,
  richiede matrice di distanza filogenetica fornita dall'utente.
- **Covariance Ratio test** per modularity / integration (Adams 2016, *Methods Ecol
  Evol* 7:565–572).
- **Auto-rilevamento landmark** combinato: outline equidistante + curvature peaks
  (style "deep-learning suggestion", utilizzabile come bootstrap iniziale prima della
  correzione manuale).

---

## 2. Funzionalità del software

- **Caricamento immagini**: JPG, PNG, TIF, PDF (tramite poppler / pdf2image),
  drag & drop multiplo.
- **Template integrati**: Rissoidae (19 LM, after Criscione & Patti 2010), bivalve
  outline (15 LM), echinoid test (12 LM), generic outline (15 sliding LM).
- **Auto-rilevamento landmark** in 3 modalità: curvature peaks (default), equidistant
  outline, Shi-Tomasi / Harris corners.
- **Editor di landmark**: aggiungi / sposta / elimina, zoom, overlay numerati / dots
  / template links / outline EFA / highlight slider.
- **Pipeline statistica completa** lanciata con un click: GPA → PCA → TPS → EFA →
  RRPP → MANOVA → Mahalanobis → Allometria + CAC + HOS → Disparity → (Kmult).
- **Output statistici visibili**: ogni risultato è mostrato sia come card numerica
  (Wilks Λ, F, Z, p, R²) sia come grafico Plotly interattivo.
- **Esportazione automatica**:
  - **CSV** multi-blocco con landmark grezzi, shape coordinates allineate, PC scores,
    autovalori, EFA scores, e tabella Procrustes ANOVA.
  - **TPS** in formato Rohlf tpsDig2 (compatibile con altri software morfometrici).
- **Cross-platform**: eseguibili nativi Windows, Linux, macOS Intel,
  **macOS Apple Silicon**.

---

## 3. Installazione (sviluppo, da sorgente)

```bash
git clone <repo>
cd pygeomorph
python -m venv venv
source venv/bin/activate           # Windows: venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Apri il browser a **http://127.0.0.1:5000**.

### Dipendenze di sistema per il supporto PDF

PyGeoMorph usa `pdf2image` per leggere `.pdf` come immagini di calibro fotografico.
Richiede `poppler` installato a livello di sistema:

| OS              | Comando                                                          |
|-----------------|------------------------------------------------------------------|
| Ubuntu / Debian | `sudo apt install poppler-utils`                                 |
| macOS           | `brew install poppler`                                           |
| Windows         | scarica i binari da [oschwartz10612/poppler-windows](https://github.com/oschwartz10612/poppler-windows/releases/) e aggiungi `bin/` al `PATH` |

---

## 4. Eseguibili automatici da GitHub Actions

Il workflow `.github/workflows/build.yml` produce automaticamente quattro binari
nativi ad ogni push di tag `v*`:

- `PyGeoMorph-windows-x64.exe`
- `PyGeoMorph-linux-x64`
- `PyGeoMorph-macos-x86_64` (Mac Intel)
- `PyGeoMorph-macos-arm64` (Mac Apple Silicon: M1, M2, M3, M4)

### 4.1 Perché due binari macOS separati?

Costruire un eseguibile macOS universal2 (Intel + ARM in un solo file) con
PyInstaller fallisce quando una qualunque wheel binaria del progetto è single-arch.
NumPy, SciPy, OpenCV e Pillow distribuite su PyPI sono normalmente single-arch:
il bundling universal2 fallisce alla validazione delle architetture con errori
del tipo `Bad CPU type in executable` o `IncompatibleBinaryArchError`.

La **soluzione standard 2024-2026** non è più tentare un build universal2, ma
**buildare separatamente** sui due runner GitHub:

- `macos-13` → runner Intel x86_64 → produce `PyGeoMorph-macos-x86_64`
- `macos-14` → runner Apple Silicon arm64 (nuovo default `macos-latest` da fine
  2024) → produce `PyGeoMorph-macos-arm64`

Questo è lo stesso pattern adottato da progetti scientifici Python come *napari*
e *scikit-image-tools*.

### 4.2 Per pubblicare una release

```bash
git tag v2.0.0
git push origin v2.0.0
```

GitHub Actions builda i 4 eseguibili in ~10 minuti e li pubblica automaticamente
come asset di una nuova GitHub Release.

---

## 5. Workflow tipico

1. **Sidebar → Acquisition**: scegli template (es. *Rissoidae 19 LM*),
   gruppo / taxon, numero di landmark, metodo di rilevamento.
2. **Carica immagini** (drag & drop, multi-file): per ogni specimen
   il software estrae automaticamente l'outline e suggerisce N landmark sulle
   massime di curvatura.
3. **Tab Digitize**: revisiona ed eventualmente correggi i landmark
   (Add / Move / Delete).
4. **Sidebar → Pipeline**:
   - imposta il numero di permutazioni RRPP (default 999),
   - mantieni **Use sliding semi-landmarks** attivo se il template ha slider
     definiti,
   - attiva **Compute Elliptic Fourier** per analisi outline complementare,
   - clicca **Run Full Pipeline**.
5. Esplora i risultati nei tab:
   - **Procrustes (GPA)** — configurazioni allineate + consensus + Procrustes
     distance plot;
   - **PCA / Shape Space** — morfospazio interattivo + scree + tabella PC scores;
   - **TPS Deformation** — griglia di deformazione consensus → estremo PC1, con
     bending energy;
   - **EFA Outline** — morfospazio basato sui descrittori di Fourier ellittici;
   - **Statistics · RRPP** — Procrustes ANOVA, MANOVA, Mahalanobis heatmap,
     disparity, Kmult (se è caricata una matrice di distanze filogenetiche);
   - **Allometry · CAC** — R², p (permutazione), HOS, traiettoria CAC vs log(CS);
   - **Raw Data** — registry completo + tabella coordinate raw.
6. **Export CSV** o **Export TPS** dalla sidebar per l'archiviazione e l'analisi
   in software esterni.

---

## 6. Citazione

Se utilizzi PyGeoMorph in una pubblicazione, ti chiediamo di citare:

> Patti, F.P. (2026). *PyGeoMorph: a modern all-in-one Python pipeline for
> geometric morphometric analysis of marine invertebrates.*
> Stazione Zoologica Anton Dohrn, Naples, Italy.

…e il lavoro originale che costituisce la base concettuale del template
malacologico:

> Criscione, F., & Patti, F.P. (2010). Similar shells are not necessarily a
> reliable guide to phylogeny: *Rissoa guerinii* Recluz, 1843, and *Rissoa lia*
> (Monterosato, 1884) (Caenogastropoda: Rissoidae): a case study.
> *The Nautilus*, **124**(3): 117–128.

---

## 7. Riferimenti bibliografici

- Adams, D.C. (2014). A generalized K statistic for estimating phylogenetic signal
  from shape and other high-dimensional multivariate data. *Syst Biol* 63: 685–697.
- Adams, D.C. (2016). Evaluating modularity in morphometric data: challenges with
  the RV coefficient and a new test measure. *Methods Ecol Evol* 7: 565–572.
- Adams, D.C., Collyer, M.L., & Kaliontzopoulou, A. (2024). *Geomorph: Software for
  geometric morphometric analyses*. R package version 4.0.8.
  https://CRAN.R-project.org/package=geomorph
- Bonhomme, V., Picq, S., Gaucherel, C., & Claude, J. (2014). Momocs: Outline
  analysis using R. *J Stat Softw* 56: 1–24.
- Bookstein, F.L. (1989). Principal warps: thin-plate splines and the
  decomposition of deformations. *IEEE PAMI* 11: 567–585.
- Bookstein, F.L. (1997). Landmark methods for forms without landmarks:
  morphometrics of group differences in outline shape. *Med Image Anal* 1: 225–243.
- Collyer, M.L., & Adams, D.C. (2018). RRPP: An R package for fitting linear
  models to high-dimensional data using residual randomization.
  *Methods Ecol Evol* 9: 1772–1779.
- Criscione, F., & Patti, F.P. (2010). Similar shells are not necessarily a
  reliable guide to phylogeny: *Rissoa guerinii* Recluz, 1843, and *Rissoa lia*
  (Monterosato, 1884) (Caenogastropoda: Rissoidae): a case study.
  *The Nautilus* 124: 117–128.
- Gunz, P., Mitteroecker, P., & Bookstein, F.L. (2005). Semilandmarks in three
  dimensions. In: D.E. Slice (ed.), *Modern Morphometrics in Physical Anthropology*,
  Kluwer.
- Kuhl, F.P., & Giardina, C.R. (1982). Elliptic Fourier features of a closed
  contour. *Comput Graph Image Proc* 18: 236–258.
- Mitteroecker, P., Gunz, P., Bernhard, M., Schaefer, K., & Bookstein, F.L.
  (2004). Comparison of cranial ontogenetic trajectories among great apes and
  humans. *J Hum Evol* 46: 679–698.
- Mitteroecker, P., & Schaefer, K. (2022). Thirty years of geometric
  morphometrics: Achievements, challenges, and the ongoing quest for biological
  meaningfulness. *Am J Biol Anthropol* 178: 181–210.
- Rohlf, F.J., & Slice, D.E. (1990). Extensions of the Procrustes method for the
  optimal superimposition of landmarks. *Syst Zool* 39: 40–59.

---

## 8. Licenza

MIT

---

## 9. Contatti

Per segnalazioni, suggerimenti, collaborazioni o richieste di nuovi template
specifici per altri taxa di invertebrati marini:

**Francesco Paolo Patti**  
Stazione Zoologica Anton Dohrn, Napoli, Italia  
✉ francesco.patti@szn.it