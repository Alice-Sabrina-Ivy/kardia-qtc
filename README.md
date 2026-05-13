# KardiaMobile 6L QTc Measurement Tool

A Python script that measures the QT interval — and computes QTc — from a
[KardiaMobile 6L](https://www.kardia.com/) PDF export. Built for at-home
QTc tracking without a clinician-review subscription.

## ⚠️ Disclaimer

**This is not a medical device.** It is an experimental, unvalidated script
for personal informational use. Its output should not be used to make
medical decisions without clinician review. If you are concerned about
QTc prolongation, get a clinical 12-lead ECG.

The detection parameters and crop coordinates were tuned on one person's
recordings made with mid-2026 Kardia software. Results may differ for
other subjects, body types, or PDF formats. Use at your own risk.

## What it does

KardiaMobile produces a 6-lead PDF export with rhythm strips. Getting QT
from that PDF normally requires either eyeball measurement on a printout
or AliveCor's paid clinician-review service. This script:

1. Rasterizes the rhythm-strip page of the PDF at 300 DPI
2. Crops Lead II at known pixel coordinates
3. Extracts the waveform from the bitmap by tracking the trace per column
4. Detects R peaks, then for each beat measures **Q onset**, **T peak**,
   and **T end** (via the tangent method — the international clinical
   standard for QT measurement)
5. Filters out beats contaminated by gridline-crossing artifacts
6. Reports mean QT and computes QTc by both **Bazett** and **Fridericia**
   formulas
7. Saves a diagnostic plot showing the full strip plus four representative
   beats with Q, T-peak, and T-end marks — so you can verify the
   measurements visually

## Who it's for

- People with a KardiaMobile 6L who want programmatic QTc readings without
  a clinician-review subscription
- People on a medication with known QT-prolonging risk who want more
  frequent surveillance between clinic visits — common examples include
  certain antibiotics (azithromycin, levofloxacin, moxifloxacin),
  antidepressants (citalopram, escitalopram at high dose), antifungals
  (fluconazole, voriconazole), antipsychotics (haloperidol, ziprasidone,
  quetiapine), antiemetics and prokinetics (ondansetron, domperidone,
  droperidol), antiarrhythmics (sotalol, dofetilide, amiodarone), and
  methadone
- People with diagnosed long QT syndrome doing supplemental home tracking
- Anyone doing quantified-self work who wants their own data
- Developers wanting a starting point for ECG signal processing on Kardia
  exports

This is supplemental data. It does not replace clinical care, and it is
not appropriate for any kind of acute or emergent situation.

## Quick start

```bash
# Python deps (one time)
pip install pillow numpy scipy matplotlib

# System dep for pdftoppm:
#   macOS:           brew install poppler
#   Ubuntu/Debian:   sudo apt install poppler-utils
#   Fedora:          sudo dnf install poppler-utils
#   Windows:         install poppler binaries, ensure pdftoppm is on PATH

# Run
python measure_qtc.py path/to/kardia-export.pdf
```

Output prints to stdout: HR, mean QT, QTc Bazett, QTc Fridericia, and an
interpretation flag. A diagnostic plot is written to
`/tmp/qtc_work/qtc_diagnostic.png` (override with `--workdir`).

**Always look at the diagnostic plot.** If the measurement marks (green Q,
magenta T peak, blue T end) are obviously wrong on the per-beat panels,
don't trust the number.

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
keep a simple log: date, HR, QTc Bazett.

## Interpreting the output

| QTc Bazett | Flag         | General interpretation                          |
|------------|--------------|-------------------------------------------------|
| <450 ms    | Normal       | Within normal range                             |
| 450–460 ms | High-normal  | Borderline; pay attention to trend              |
| 460–470 ms | Borderline   | Discuss with prescriber                         |
| >470 ms    | Prolonged    | Clinical evaluation warranted                   |

These thresholds are commonly cited for adult females; adult male upper
limits are typically ~10 ms shorter. Your prescriber may use different
thresholds depending on context (e.g., specific medications, congenital
long QT).

**Change from baseline often matters more than the absolute value.** A
commonly used criterion when starting a QT-prolonging drug: investigate
further if QTc increases by >60 ms above personal baseline, even if the
absolute value is still in the "normal" range. This is why a multi-reading
baseline matters — take 3–5 readings over 1–2 weeks before starting any
new medication, average them, and compare future readings to that average.

## How it works (technical)

1. **Rasterize** the rhythm-strip page from the PDF at 300 DPI via
   `pdftoppm`. Pages 2–5 of a Kardia 6L PDF all contain the same 6 leads,
   just different time windows. Default is page 4 (use `--page` to change).

2. **Crop Lead II** at pixel coordinates `(50, 850, 2500, 1250)`. Lead II
   is preferred for QT measurement because P, QRS, and T waves are
   typically prominent and the heart's electrical axis usually aligns
   with it.

3. **Extract waveform** by finding the topmost cluster of dark pixels per
   column. The "first cluster from top" heuristic ignores most gridline
   noise; the rest is filtered downstream.

4. **Smooth** with a Savitzky-Golay filter (window 11, polyorder 3) to
   enable stable derivative computation. T waves can be very low amplitude
   in Lead II (often ~10 px ≈ 0.085 mV), and without heavy smoothing the
   gradient is dominated by quantization noise.

5. **Detect R peaks** with `scipy.signal.find_peaks` — minimum amplitude
   50 px, minimum spacing 150 px (~500 ms; allows HR up to ~120 BPM).

6. **For each beat, measure QT**:
   - **Q onset**: walk backward from the R peak until slope drops below
     0.8 px/px (no longer rapid upstroke).
   - **T peak**: maximum value in the 200–450 ms window after segment
     start (the segment starts ~135 ms before R, so the T-peak search
     window is effectively ~65–315 ms after R).
   - **T end via tangent method**: find the steepest descending slope in
     the 0–150 ms after T peak; extrapolate that tangent line to baseline;
     intersection is T end.

7. **Filter contaminated beats** using density-based cluster finding (see
   "Gotchas" below for why this is needed).

8. **Compute QTc** with both Bazett (QT/√RR) and Fridericia (QT/RR^⅓).
   Bazett over-corrects at high HR and under-corrects at low HR;
   Fridericia is more stable across HRs but Bazett remains the clinical
   standard.

## Verifying the crop is correct (do this first!)

If you're using this on your own Kardia exports for the first time, the
single most important sanity check is: **look at `lead_ii.jpg` in the
work directory after a run.** It should show:

- Only Lead II (no other leads visible)
- A clean black trace on white background
- QRS spikes pointing up
- Small positive T waves following each QRS
- Roughly 8–10 beats visible across the strip

If you see multiple leads, the wrong lead, or a clipped trace, adjust
the coordinates in `crop_lead_ii()`. AliveCor occasionally changes their
PDF layout, and a different version of the Kardia app may render
differently.

## Gotchas (and why the code is the way it is)

These bit me during development; the script handles them, but if you
modify the code, watch out:

### 1. The "II" lead label contaminates beat 1

The "II" text label printed in the rhythm strip overlaps the first
beat's area. The letterform's dark pixels get extracted as if they were
ECG signal, producing spurious spikes that fool the algorithm.

**Fix**: beat index 0 is always dropped before further filtering.

### 2. T waves in Lead II can be very low amplitude

T-wave amplitude is often a small fraction of QRS amplitude (e.g.,
~10 px vs. ~140 px for QRS). At that size, naive slope detection finds
noise spikes rather than the real T-wave descent.

**Fix**: Heavy Savitzky-Golay smoothing (window 11) before any derivative
work. Without it, the gradient is dominated by quantization noise.

### 3. Gridline artifacts produce false-short QT values

Vertical gridlines on the ECG paper occasionally cross the trace at
locations where the trace happens to dip. The waveform extractor picks
these up as dark pixels, producing false sharp spikes that look like
early T-wave descents. The algorithm then reports a falsely short QT.

This is **asymmetric** — contamination always *shortens* measured QT
(false early termination), never lengthens it.

**Fix**: Density-based cluster finding for beat filtering. Simple
median + tolerance filters fail when contamination is heavy because the
median falls between the clean cluster and the noise. The density
approach finds the densest cluster of mutually-close values, which is
robust to scattered contamination.

### 4. Q onset detection is tricky

The QRS upstroke is sharp but its very start (where ventricular
depolarization begins) is gradual. Walking *forward* from baseline
looking for "slope > X" fails because the threshold is ambiguous.

**Fix**: Walk *backward* from the R peak; find the first point where
slope drops below 0.8 px/px. This reliably catches the transition from
rapid upstroke back to baseline.

### 5. Calibration depends on render DPI

The pixel-to-time conversion (3.387 ms/px) is specific to 300 DPI
rendering of a standard 25 mm/s ECG. If you change `DPI` in the code,
the calibration constants auto-update. Don't hardcode time values in
pixels anywhere.

### 6. Crop coordinates are Kardia-specific

The Lead II crop is based on KardiaMobile 6L PDF layout as of mid-2026.
If AliveCor changes the PDF format, these will need adjustment. See
"Verifying the crop is correct" above.

## Limitations

- **Single readings are noisy.** Pixel-based measurement on a printed
  tracing has roughly ±10–15 ms uncertainty per beat. Average 3+ beats
  per reading (which the script does) and at least 3 readings over time
  for a more reliable estimate.

- **Not a substitute for clinical ECG.** If QTc is reported as elevated,
  or trending upward by >30 ms, get a clinical 12-lead ECG with
  cardiologist review before making drug decisions.

- **Lead II only.** Some QT prolongation patterns (notched T waves,
  prominent U waves) are better appreciated on precordial leads, which
  Kardia 6L recordings do not include.

- **Tangent method systematically underestimates** compared to threshold-
  based "T end at baseline crossing." This is by design — tangent is the
  international clinical standard. It reports ~10–20 ms shorter than
  threshold methods, which is fine for trend tracking as long as the same
  method is used consistently.

- **Validated on one subject's recordings only.** Detection parameters
  were tuned empirically. Other people's T-wave amplitudes, QRS
  morphology, or rhythm irregularities may require parameter tuning.
  Always inspect the diagnostic plot.

## Files

- `measure_qtc.py` — main script
- `README.md` — this file

## Contributing

PRs welcome, particularly for:

- Other Kardia models (single-lead KardiaMobile, future 12-lead devices)
- Other PDF formats (AliveCor layout changes, third-party ECG patch
  exports)
- Better T-end detection (alternative tangent variants, ML approaches)
- Automated crop verification

## License

[Add a license of your choice — MIT, Apache 2.0, GPL, etc.]
