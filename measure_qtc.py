#!/usr/bin/env python3
"""
KardiaMobile 6L QTc Measurement Tool
=====================================

Measures QTc from a KardiaMobile 6L PDF export using the tangent method on Lead II.

Usage:
    python measure_qtc.py <path_to_kardia_pdf>

Approach (refined through empirical trial-and-error on real Kardia exports):
  1. Rasterize one of the rhythm-strip pages (page 4 by default) at 300 DPI
  2. Crop Lead II using known coordinates for KardiaMobile PDF layout
  3. Extract waveform by finding the top edge of dark pixels per column
  4. Smooth with Savitzky-Golay filter to enable stable derivatives
  5. Find R peaks via amplitude + minimum-distance criteria
  6. For each beat: find Q onset (slope-based, walking back from R)
                   find T peak (max in 200-450ms window after R)
                   find T end (tangent method: steepest desc slope extrapolated to baseline)
  7. Exclude contaminated beats (gridline-spike artifacts, "II" lead label in first beat)
  8. Report mean/median QT and compute QTc via Bazett + Fridericia formulas

Calibration (validated empirically):
  - 300 DPI rendering of standard 25mm/s ECG gives ~11.81 px/mm horizontal
  - 1 small square (1mm = 40ms) = 11.81 px, so 1 px = 3.387 ms

Lead II crop coords for KardiaMobile 6L PDF rendered at 300 DPI:
  - Page is 2550 x 3300 pixels (8.5" x 11" at 300 DPI)
  - Lead II strip is the 2nd of 6 strips; on page 4 it's roughly y=850-1250
  - These coords work for KardiaMobile 6L exports — adjust if Kardia changes layout

Known gotchas (these all came up in development; fixes are in the code):
  - "II" lead label embedded in beat 1 area creates spurious dark pixels that look
    like extra QRS spikes. Exclude beat 1 by default (use indices 1+).
  - At column granularity, some "dark pixels" come from gridlines. The "first cluster
    from top" heuristic filters most but not all. Beats with visible vertical-line
    artifacts in the result plot should be excluded manually.
  - T waves in Lead II can be low amplitude (often only ~10 px ≈ 0.085 mV). Naive
    slope detection can miss them; need at least 11-sample Savitzky-Golay smoothing.
  - QRS onset detection: walk BACK from the R peak until slope drops below ~0.8 px/px.
    Simple threshold methods (slope > X) walking forward fail because the segment
    starts in flat baseline.

Dependencies: pillow, numpy, scipy, matplotlib (Python) + pdftoppm from poppler (system)
Install: pip install pillow numpy scipy matplotlib
         macOS: brew install poppler
         Debian/Ubuntu: sudo apt install poppler-utils
"""

import sys
import os
import subprocess
import argparse
import numpy as np
from PIL import Image
from scipy.signal import find_peaks, savgol_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# === Calibration constants (300 DPI, 25mm/s ECG paper) ===
DPI = 300
MM_PER_INCH = 25.4
PX_PER_MM = DPI / MM_PER_INCH  # 11.81
MS_PER_PX = 40 / PX_PER_MM     # 3.387 ms/px
PX_PER_MS = 1 / MS_PER_PX


def rasterize_page(pdf_path, page_num, output_path):
    """Render a single PDF page to JPEG at 300 DPI."""
    result = subprocess.run(
        ['pdftoppm', '-jpeg', '-r', str(DPI), '-f', str(page_num), '-l', str(page_num),
         pdf_path, output_path.replace('.jpg', '')],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"pdftoppm failed: {result.stderr}")
    # pdftoppm uses zero-padded names based on total page count
    # For a 5-page PDF that's pageX-N.jpg (no zero pad needed since total < 10)
    candidates = [
        f"{output_path.replace('.jpg', '')}-{page_num}.jpg",
        f"{output_path.replace('.jpg', '')}-0{page_num}.jpg",
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    raise FileNotFoundError(f"Output not found for page {page_num}")


def crop_lead_ii(page_image_path, output_crop_path):
    """
    Crop Lead II strip from a rasterized KardiaMobile rhythm-strip page.

    Lead II is the 2nd of 6 strips. On rhythm pages of a 300 DPI render, it
    occupies roughly y=850-1250. Width spans full strip area x=50-2500.

    If Kardia changes their PDF layout, adjust these coordinates.
    """
    img = Image.open(page_image_path)
    crop = img.crop((50, 850, 2500, 1250))
    crop.save(output_crop_path)
    return output_crop_path


def extract_waveform(lead_image_path):
    """
    Extract ECG waveform y-coordinate per column from a Lead II crop.

    For each column, find the topmost cluster of dark pixels (the ECG trace).
    The "first cluster" heuristic ignores stray dark pixels from gridlines.
    """
    img = Image.open(lead_image_path)
    arr = np.array(img.convert('L'))
    H, W = arr.shape
    trace_mask = arr < 100  # dark pixels = trace

    waveform = np.full(W, np.nan)
    for x in range(W):
        dark_y = np.where(trace_mask[:, x])[0]
        if len(dark_y) == 0:
            continue
        # Group consecutive dark pixels; take topmost cluster's top y
        gaps = np.where(np.diff(dark_y) > 3)[0]
        first_cluster = dark_y[:gaps[0] + 1] if len(gaps) > 0 else dark_y
        if len(first_cluster) >= 1:
            waveform[x] = first_cluster[0]

    # Interpolate any NaN columns and smooth
    valid = ~np.isnan(waveform)
    xx = np.arange(W)
    waveform = np.interp(xx, xx[valid], waveform[valid])
    waveform = savgol_filter(waveform, window_length=11, polyorder=3)

    # Invert so positive deflections go up; zero baseline
    ecg = -waveform
    ecg = ecg - np.median(ecg)
    return ecg


def detect_r_peaks(ecg):
    """Find R peaks. Min height 50px, min separation 150px (~500ms)."""
    r_peaks, _ = find_peaks(ecg, height=50, distance=150)
    return r_peaks


def measure_beat(ecg, r, r_peaks, beat_idx):
    """
    Measure QT for a single beat centered at R peak position `r`.

    Returns dict with q_onset, t_peak, t_end positions (in pixel coords)
    and qt_ms, or None if measurement fails.
    """
    W = len(ecg)
    start = max(0, r - 40)  # 40 px before R = ~135 ms
    end = min(W, r + int(550 * PX_PER_MS))
    seg = ecg[start:end]
    if len(seg) < 100:
        return None

    seg_slopes = np.gradient(seg)

    # === Q onset: walk backward from R peak (at idx 40) ===
    r_idx_in_seg = 40 if r >= 40 else r  # adjust if near segment start
    q_idx = r_idx_in_seg
    for j in range(r_idx_in_seg - 1, 0, -1):
        if seg_slopes[j] < 0.8:  # slope dropped below rapid-upstroke threshold
            q_idx = j
            break
    q_onset_ms = q_idx * MS_PER_PX

    # === T peak: max in 200-450ms window from segment start ===
    tp_start = int(200 * PX_PER_MS)
    tp_end = min(len(seg), int(450 * PX_PER_MS))
    if tp_end <= tp_start:
        return None
    t_peak_local = tp_start + np.argmax(seg[tp_start:tp_end])
    t_peak_ms = t_peak_local * MS_PER_PX

    # === T end: tangent method ===
    desc_start = t_peak_local + 2
    desc_end = min(len(seg), t_peak_local + int(150 * PX_PER_MS))
    desc_seg = seg[desc_start:desc_end]
    if len(desc_seg) < 5:
        return None

    desc_slopes = np.gradient(desc_seg)
    if len(desc_slopes) > 5:
        win = min(5, len(desc_slopes) // 2 * 2 + 1)
        if win >= 3:
            desc_slopes = savgol_filter(desc_slopes, window_length=win, polyorder=2)

    steepest = int(np.argmin(desc_slopes))
    steepest_idx = desc_start + steepest
    steepest_slope = desc_slopes[steepest]

    # Baseline from end of segment (TP segment)
    post_baseline = np.median(seg[-30:])

    if steepest_slope < -0.1:
        # Tangent line: y = seg[steepest_idx] + slope * (x - steepest_idx)
        # Solve for y = post_baseline:
        t_end_idx = steepest_idx + (post_baseline - seg[steepest_idx]) / steepest_slope
        t_end_ms = float(t_end_idx) * MS_PER_PX
    else:
        # T wave too flat for tangent — fall back to baseline crossing
        t_end_idx = t_peak_local
        for k in range(t_peak_local, len(seg)):
            if seg[k] <= post_baseline + 2:
                t_end_idx = k
                break
        t_end_ms = t_end_idx * MS_PER_PX

    qt_ms = t_end_ms - q_onset_ms

    return {
        'beat_idx': beat_idx,
        'q_onset_ms': q_onset_ms,
        't_peak_ms': t_peak_ms,
        't_end_ms': t_end_ms,
        'qt_ms': qt_ms,
        'seg': seg,
        'q_idx_seg': q_idx,
        't_peak_idx_seg': t_peak_local,
        't_end_idx_seg': t_end_idx if steepest_slope < -0.1 else t_end_idx,
        'post_baseline': post_baseline,
        'tangent_valid': steepest_slope < -0.1,
    }


def filter_contaminated_beats(measurements):
    """
    Drop contaminated beats using density-based cluster finding.

    Key insight from empirical work: contamination is ASYMMETRIC. Gridline
    artifacts (vertical grid lines crossing the trace) get detected as false
    early T-wave descents, which produces shorter-than-true QT values.
    Real beat-to-beat QT variation is small (~10-20 ms). So contaminated
    beats appear as scattered short values while clean beats cluster tightly.

    A simple median+tolerance filter fails when contamination is heavy because
    the median falls between the clean cluster and the noise. Instead:
      - Find the densest cluster: for each value, count peers within ±20 ms
      - Use the highest-density value as the anchor
      - Keep all values within ±20 ms of the anchor

    Also drop beat 0 (lead label "II" artifact contaminates this beat).
    """
    if not measurements:
        return []

    filtered = [m for m in measurements if m['beat_idx'] != 0]
    if len(filtered) < 3:
        return filtered

    qts = np.array([m['qt_ms'] for m in filtered])
    TOLERANCE_MS = 20

    # Density score for each value: how many peers are within tolerance
    densities = np.array([
        np.sum(np.abs(qts - q) <= TOLERANCE_MS) for q in qts
    ])

    # Anchor: value with highest density (ties broken by larger QT since
    # contamination only shortens)
    max_density = densities.max()
    candidates = np.where(densities == max_density)[0]
    # Among ties, pick the LONGEST QT (safer assumption given contamination direction)
    anchor_idx = candidates[np.argmax(qts[candidates])]
    anchor_qt = qts[anchor_idx]

    keep_mask = np.abs(qts - anchor_qt) <= TOLERANCE_MS
    return [m for m, k in zip(filtered, keep_mask) if k]


def compute_qtc(qt_ms, rr_ms):
    """Compute Bazett and Fridericia QTc."""
    rr_s = rr_ms / 1000.0
    qtc_b = qt_ms / np.sqrt(rr_s)
    qtc_f = qt_ms / (rr_s ** (1 / 3))
    return qtc_b, qtc_f


def plot_results(ecg, r_peaks, all_measurements, kept_measurements, output_path):
    """Diagnostic plot: full strip + 4 individual beats."""
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(3, 4)

    # Full strip on top row
    ax_full = fig.add_subplot(gs[0, :])
    ax_full.plot(ecg, color='black', linewidth=0.6)
    ax_full.axhline(0, color='gray', linewidth=0.3)
    ax_full.plot(r_peaks, ecg[r_peaks], 'rv', markersize=6)
    kept_idxs = {m['beat_idx'] for m in kept_measurements}
    for m in all_measurements:
        color = 'blue' if m['beat_idx'] in kept_idxs else 'gray'
        alpha = 1.0 if m['beat_idx'] in kept_idxs else 0.3
        r = r_peaks[m['beat_idx']]
        start = r - 40
        q = start + m['q_idx_seg']
        t = start + int(m['t_end_idx_seg'])
        ax_full.plot([q, t], [-35, -35], '-', color=color, linewidth=2, alpha=alpha)
        ax_full.text((q + t) / 2, -45, f"#{m['beat_idx']+1}",
                     ha='center', fontsize=8, color=color, alpha=alpha)
    ax_full.set_title(f"Full Lead II strip — kept beats in blue, excluded in gray")
    ax_full.grid(True, alpha=0.3)

    # 4 individual beats on bottom rows (prefer kept beats)
    to_plot = kept_measurements[:4] if len(kept_measurements) >= 4 else (
        kept_measurements + [m for m in all_measurements if m not in kept_measurements]
    )[:4]
    for i, m in enumerate(to_plot):
        ax = fig.add_subplot(gs[1 + i // 2, (i % 2) * 2:(i % 2) * 2 + 2])
        seg = m['seg']
        t_ms = np.arange(len(seg)) * MS_PER_PX
        ax.plot(t_ms, seg, 'k-', linewidth=1.2)
        ax.axhline(0, color='gray', linewidth=0.4)
        ax.axhline(m['post_baseline'], color='lightgray', linestyle=':', linewidth=0.5)
        ax.axvline(m['q_onset_ms'], color='green', linewidth=1.5,
                   label=f"Q ({m['q_onset_ms']:.0f}ms)")
        ax.axvline(m['t_peak_ms'], color='magenta', linewidth=1, linestyle='--',
                   label=f"T pk ({m['t_peak_ms']:.0f}ms)")
        ax.axvline(m['t_end_ms'], color='blue', linewidth=1.5,
                   label=f"T end ({m['t_end_ms']:.0f}ms)")
        kept_str = "KEPT" if m['beat_idx'] in kept_idxs else "EXCLUDED"
        ax.set_title(f"Beat {m['beat_idx']+1} [{kept_str}] — QT = {m['qt_ms']:.0f} ms")
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlabel('ms')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=100)
    plt.close()


def measure_qtc_from_pdf(pdf_path, page_num=4, workdir='/tmp/qtc_work'):
    """
    Full pipeline: PDF → QTc result dict.

    Returns:
        {
          'hr_bpm', 'rr_mean_ms',
          'qt_mean_ms', 'qt_median_ms', 'qt_values',
          'qtc_bazett', 'qtc_fridericia',
          'beats_measured', 'beats_kept',
          'diagnostic_image_path',
        }
    """
    os.makedirs(workdir, exist_ok=True)
    page_image = os.path.join(workdir, 'page.jpg')
    lead_image = os.path.join(workdir, 'lead_ii.jpg')
    plot_path = os.path.join(workdir, 'qtc_diagnostic.png')

    page_image = rasterize_page(pdf_path, page_num, page_image)
    crop_lead_ii(page_image, lead_image)
    ecg = extract_waveform(lead_image)

    r_peaks = detect_r_peaks(ecg)
    if len(r_peaks) < 3:
        raise RuntimeError(f"Only found {len(r_peaks)} R peaks; need at least 3")

    rr_ms = np.diff(r_peaks) * MS_PER_PX
    rr_mean_ms = float(rr_ms.mean())
    hr_bpm = 60000.0 / rr_mean_ms

    all_measurements = []
    for i, r in enumerate(r_peaks):
        m = measure_beat(ecg, r, r_peaks, i)
        if m is not None:
            all_measurements.append(m)

    kept = filter_contaminated_beats(all_measurements)
    if len(kept) < 2:
        raise RuntimeError(f"Only {len(kept)} clean beats; result unreliable")

    qts = np.array([m['qt_ms'] for m in kept])
    qt_mean = float(qts.mean())
    qt_median = float(np.median(qts))

    qtc_b, qtc_f = compute_qtc(qt_mean, rr_mean_ms)

    plot_results(ecg, r_peaks, all_measurements, kept, plot_path)

    return {
        'hr_bpm': hr_bpm,
        'rr_mean_ms': rr_mean_ms,
        'qt_mean_ms': qt_mean,
        'qt_median_ms': qt_median,
        'qt_values': qts.tolist(),
        'qtc_bazett': float(qtc_b),
        'qtc_fridericia': float(qtc_f),
        'beats_measured': len(all_measurements),
        'beats_kept': len(kept),
        'diagnostic_image_path': plot_path,
    }


def print_report(result):
    """Pretty-print the measurement result."""
    print()
    print("=" * 60)
    print("  KARDIA 6L QTc MEASUREMENT REPORT")
    print("=" * 60)
    print(f"  Heart rate:         {result['hr_bpm']:.1f} BPM")
    print(f"  Mean RR interval:   {result['rr_mean_ms']:.0f} ms")
    print()
    print(f"  Beats measured:     {result['beats_measured']}")
    print(f"  Beats kept (clean): {result['beats_kept']}")
    print(f"  QT values (ms):     {[round(q,1) for q in result['qt_values']]}")
    print(f"  Mean QT:            {result['qt_mean_ms']:.0f} ms")
    print(f"  Median QT:          {result['qt_median_ms']:.0f} ms")
    print()
    print(f"  QTc Bazett:         {result['qtc_bazett']:.0f} ms")
    print(f"  QTc Fridericia:     {result['qtc_fridericia']:.0f} ms")
    print()
    # Interpretation guide
    qtcb = result['qtc_bazett']
    if qtcb < 450:
        flag = "NORMAL"
    elif qtcb < 460:
        flag = "HIGH-NORMAL"
    elif qtcb < 470:
        flag = "BORDERLINE PROLONGED"
    else:
        flag = "PROLONGED"
    print(f"  Interpretation:     {flag}")
    print(f"  (Common reference thresholds: 450 ms caution, 470 ms prolonged.")
    print(f"   Not a medical device; consult your clinician.)")
    print()
    print(f"  Diagnostic plot:    {result['diagnostic_image_path']}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Measure QTc from KardiaMobile 6L PDF")
    parser.add_argument('pdf', help="Path to Kardia 6L PDF export")
    parser.add_argument('--page', type=int, default=4,
                        help="PDF page to analyze (default 4; pages 2-5 are rhythm strips)")
    parser.add_argument('--workdir', default='/tmp/qtc_work',
                        help="Working directory for intermediate files")
    args = parser.parse_args()

    if not os.path.exists(args.pdf):
        print(f"Error: PDF not found at {args.pdf}", file=sys.stderr)
        sys.exit(1)

    result = measure_qtc_from_pdf(args.pdf, page_num=args.page, workdir=args.workdir)
    print_report(result)


if __name__ == '__main__':
    main()
