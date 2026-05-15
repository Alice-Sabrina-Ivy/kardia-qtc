# KardiaMobile 6L QTc Measurement Tool

A Python script that measures the QT interval — and computes QTc — from
a [KardiaMobile 6L](https://www.kardia.com/) PDF export. Parses the
PDF's internal vector path data directly to recover the original 300 Hz
ECG samples, then measures QT on all leads with usable T waves.

## ⚠️ Disclaimer

**This is not a medical device.** It is an experimental, unvalidated
script for personal informational use. Its output should not be used to
make medical decisions without clinician review. If you are concerned
about QTc prolongation, get a clinical 12-lead ECG.

The detection parameters and PDF layout constants were tuned on
recordings from mid-2026 Kardia software. Results may differ for other
subjects, body types, or PDF formats. Use at your own risk.

## What it does

KardiaMobile produces a 6-lead PDF export with rhythm strips. Getting
QT from that PDF normally requires either eyeball measurement on a
printout or AliveCor's paid clinician-review service. This script:

1. Opens the PDF and decompresses its content streams
2. Extracts ECG samples directly from the PDF's vector path operators
   — recovering the original 300 Hz waveform without any image
   processing or rasterization
3. Partitions samples into all 6 leads by Y-coordinate
4. Detects R peaks in Leads II, III, and aVF
5. For each beat, measures **Q onset**, **T peak**, and **T end** (via
   the tangent method — the international clinical standard)
6. Filters out beats with unreliable T-wave detection
7. Reports mean QT, per-lead breakdown, and QTc by both **Bazett** and
   **Fridericia** formulas — with a standard error estimate
8. Saves a diagnostic plot showing all 4 rhythm-strip pages plus
   per-lead QT distributions, so you can verify the measurements
   visually

## Who it's for

- People with a KardiaMobile 6L who want programmatic QTc readings
  without a clinician-review subscription
- People on a medication with known QT-prolonging risk who want more
  frequent surveillance between clinic visits — common examples include
  certain antibiotics (azithromycin, levofloxacin, moxifloxacin),
  antidepressants (citalopram, escitalopram at high dose), antifungals
  (fluconazole, voriconazole), antipsychotics (haloperidol, ziprasidone,
  quetiapine), antiemetics and prokinetics (ondansetron, domperidone,
  droperidol), antiarrhythmics (sotalol, dofetilide, amiodarone), and
  methadone
- People with diagnosed long QT syndrome doing supplemental home
  tracking
- Anyone doing quantified-self work who wants their own data
- Developers wanting a starting point for ECG signal processing on
  Kardia exports

This is supplemental data. It does not replace clinical care, and it
is not appropriate for any kind of acute or emergent situation.

## Quick start

```bash
# All deps are pip-installable — no system packages required
pip install pypdf numpy scipy matplotlib

# Run
python measure_qtc.py path/to/kardia-export.pdf
```

Output prints to stdout: HR, mean QT, per-lead breakdown, QTc Bazett
(with standard error), QTc Fridericia, and an interpretation flag. A
diagnostic plot is written to `/tmp/qtc_work/qtc_diagnostic.png`
(override with `--workdir`).

**Always look at the diagnostic plot.** If the per-lead histograms
disagree wildly, or the Lead II strip shows wrong R-peak detection
(red triangles in obviously wrong places), don't trust the number.

## Taking consistent readings

For trending, consistency matters more than absolute accuracy:

- **Same time of day** (morning is usually cleanest — fewer confounds
  from caffeine, meals, activity)
- **Seated and calm for 5+ minutes** before recording
- **No caffeine, stimulants, or strenuous activity** in the prior hour
- **Standard Kardia 6L hand position** — both thumbs on the top
  electrodes, third electrode on left ankle or knee
- **30-second recording** (Kardia default)

Save each PDF with a date in the filename. Run the script on each and
keep a simple log: date, HR, QTc Bazett ± SE.

## Interpreting the output

| QTc Bazett | Flag         | General interpretation                          |
|------------|--------------|-------------------------------------------------|
| <450 ms    | Normal       | Within normal range                             |
| 450–460 ms | High-normal  | Borderline; pay attention to trend              |
| 460–470 ms | Borderline   | Discuss with prescriber                         |
| >470 ms    | Prolonged    | Clinical evaluation warranted                   |

These thresholds are commonly cited for adult females; adult male
upper limits are typically ~10 ms shorter. Your prescriber may use
different thresholds depending on context (e.g., specific medications,
congenital long QT).

**Change from baseline often matters more than the absolute value.** A
commonly used criterion when starting a QT-prolonging drug: investigate
further if QTc increases by >60 ms above personal baseline, even if
the absolute value is still in the "normal" range.

With the vector method, per-reading precision is around ±2 ms SE on
QTc Bazett, which means a real change of ~15–20 ms from baseline is
detectable — much tighter than the clinical 60 ms threshold. Take
3–5 readings over 1–2 weeks for a reliable personal baseline before
starting any new medication, then compare future readings to that
average.

## Why vector parsing rather than image processing

A more obvious approach would be to rasterize the PDF and detect the
ECG trace from pixels. That works (an earlier version of this script
did exactly that), but the vector approach is substantially better
because the Kardia PDF stores the ECG as actual sample-by-sample
vector path operators, not as a rendered image. Parsing them directly
recovers the original ADC samples.

| Metric              | Raster approach    | Vector approach    |
|---------------------|--------------------|--------------------|
| Pages used          | 1 of 4             | 4 of 4             |
| Leads used          | Lead II only       | II, III, aVF       |
| Beats per reading   | 4–6                | 30–40              |
| QT precision        | ±10–15 ms/beat     | ±2–5 ms/beat       |
| QTc Bazett SE       | ~±10 ms            | ~±2 ms             |
| Smoothing required  | yes (Savgol)       | no                 |
| Gridline noise      | filtered post-hoc  | absent (by color)  |

The vector approach reports systematically ~10 ms lower QTc than
rasterization because pixel-based trace detection inflates T-end
estimates slightly through trace thickness and gridline edge effects.

## How it works (technical)

The Kardia PDF is produced by Skia and stores the ECG as PDF vector
path operators inside its content stream, not as a rasterized image:

```
x y m          # moveto - start a path at (x, y)
x y l          # lineto - draw line to (x, y)
S              # stroke - render the path
```

Each `l` (lineto) endpoint is one ECG sample. The script's pipeline:

1. **Open the PDF with pypdf.** For each of the 4 rhythm-strip pages
   (PDF pages 2–5), decompress the content stream.

2. **Identify ECG paths by stroke color.** The content stream contains
   gridlines (light purple, `0.8196 0.8196 0.9608 RG`), text (various
   grays), and ECG samples (pure black, `0 0 0 RG`). The ECG color is
   set once just before the ECG paths, so everything after it is ECG.

3. **Extract sample points.** Each `l` operator endpoint is one ECG
   sample. Coordinates are in PDF points (1/72 inch). About 14,000
   samples per page across 6 leads = ~2,400 samples per lead per page.

4. **Partition by lead via Y-coordinate.** The 6 leads occupy strips
   at Y centers 172, 260, 352, 441, 529, 621 PDF points (top-to-bottom:
   I, II, III, aVR, aVL, aVF). Each strip's samples fall within ±40
   points of its center, so partitioning is unambiguous.

5. **Convert to physical units.** At 25 mm/s paper speed and 10 mm/mV
   amplitude:
   - 1 PDF point horizontally = 14.11 ms (so sample spacing 0.2362
     points = 3.33 ms = 300 Hz, AliveCor's native rate)
   - 1 PDF point vertically = 0.0353 mV

6. **Detect R peaks** in Leads II, III, aVF using `scipy.signal.find_peaks`
   with minimum height 0.5 mV (0.3 mV in Lead III, which has smaller R
   waves).

7. **For each beat, measure Q onset, T peak, T end:**
   - **Q onset:** walk backward from R until slope drops below
     0.02 mV/sample (6 mV/s) — the rapid-upstroke threshold
   - **T peak:** maximum amplitude in the 150–450 ms window after R
   - **T end (tangent method):** find the steepest descending slope
     in the 0–150 ms after T peak; extrapolate that tangent line to
     the post-T baseline; intersection is T end

8. **Filter contaminated beats** using density-based clustering. Beats
   within ±20 ms of the densest cluster center are kept; outliers are
   dropped. This is robust to T-wave detection failures in any single
   beat or lead.

9. **Compute QTc** with Bazett (QT/√RR) and Fridericia (QT/RR^⅓)
   formulas. RR is computed from ALL Lead-II R-peak intervals (not just
   the filtered beats), so HR isn't biased by which beats had measurable
   QT.

## Verifying the leads are identified correctly (do this first!)

If you're running this on your own Kardia exports for the first time,
the single most important sanity check is to **look at the diagnostic
plot** after a run (default `/tmp/qtc_work/qtc_diagnostic.png`). The
top 4 panels show Lead II from each rhythm page. You should see:

- Clean P-QRS-T morphology with QRS spikes pointing **up** and small
  positive T waves
- Red downward triangles (▽) marking R peaks — these should land on
  every R peak and only on R peaks
- Roughly 8–10 beats per page strip

If R peaks are wrong, leads look noisy, or the script errors with "no
ECG color block found," AliveCor may have changed their PDF format and
the constants in the script (lead Y-centers, color regex) will need
adjusting.

## Gotchas (and why the code is the way it is)

These came up during development; the code handles them, but if you
modify it, watch out:

### 1. Color identification is critical

The ECG is "black" (`0 0 0 RG`), but the PDF stream sets stroke color
multiple times for text rendering before getting to the ECG. The trick
is that the LAST color set before path drawing begins en masse is the
ECG color, and that color block contains thousands of m/l operators
(vs. a few for gridlines or text underlines). If AliveCor ever uses a
different color for the ECG, the regex `\b0\s+0\s+0\s+RG\b` will need
updating.

### 2. Page 2 has a calibration pulse

The first rhythm page (page 2 of the PDF) starts with a brief 1 mV
square-wave calibration pulse. This shows up in the extracted data
and would fool R-peak detection if treated as a beat. The script
drops the first and last beat of each page automatically to handle
this and other edge effects.

### 3. Y-center positions are page-layout-specific

The Y centers (172, 260, 352, 441, 529, 621) are based on observed
KardiaMobile 6L PDF layout as of mid-2026. If AliveCor changes the
layout, these constants will need updating.

### 4. Lead III T waves can be very small

Lead III often has T-wave amplitude near the 0.05 mV detection floor,
so many beats in Lead III get skipped (commonly 3 measurable per
reading vs. ~20 in Lead II). This is correct behavior — better to
skip a beat than report an unreliable measurement. The multi-lead
strategy means Leads II and aVF carry the bulk of the measurements.

### 5. Voltage baseline is the per-strip median

Each lead strip's vertical baseline ("0 mV") is computed as the median
Y of all that lead's samples. This works because most samples on an
ECG are at baseline (between beats); the QRS and T waves are a small
fraction of total samples and don't shift the median much. If the
baseline drifts during the recording, a more sophisticated approach
(rolling median) would be needed.

### 6. The Bazett standard error is approximate

The `qtc_bazett_se` value is `(QT_std / sqrt(N)) / sqrt(RR_s)`, which
is the SE of the QT mean propagated through the Bazett correction
**assuming RR is known exactly**. This is a slight underestimate of
the true SE (RR has its own uncertainty), but the underestimate is
small because RR uncertainty contributes much less than QT uncertainty
at typical heart rates. Treat the reported SE as a lower bound on
uncertainty.

## Limitations

- **Not a substitute for clinical ECG.** If QTc is reported as
  elevated, or trending upward by >30 ms from baseline, get a clinical
  12-lead ECG with cardiologist review before making drug decisions.

- **Lead II, III, aVF only.** Some QT prolongation patterns (notched T
  waves, prominent U waves) are better appreciated on precordial
  leads, which Kardia 6L recordings do not include.

- **Tangent method systematically underestimates** compared to
  threshold-based "T end at baseline crossing." This is by design —
  tangent is the international clinical standard. It reports
  ~10–20 ms shorter than threshold methods, which is fine for trend
  tracking as long as the same method is used consistently.

- **Validated on one subject's recordings only.** Detection parameters
  were tuned empirically. Other people's T-wave amplitudes, QRS
  morphology, or rhythm irregularities may require parameter tuning.
  Always inspect the diagnostic plot.

- **Tied to current Kardia PDF layout.** The vector-parsing approach
  depends on stable PDF structure (color codes, lead Y-centers,
  content-stream organization). Any AliveCor format change could
  require code updates.

## Files

- `measure_qtc.py` — main script
- `README.md` — this file

## Contributing

PRs welcome, particularly for:

- Other Kardia models (single-lead KardiaMobile, future 12-lead devices)
- AliveCor PDF format updates (different layout versions, color codes)
- Better T-end detection (alternative tangent variants, ML approaches)
- Per-lead weighted averaging based on T-wave SNR
- Automated lead-position auto-detection (instead of hardcoded
  Y-centers)

## License

[MIT](LICENSE) — do whatever you want with this, just keep the copyright
notice. No warranty, express or implied.
