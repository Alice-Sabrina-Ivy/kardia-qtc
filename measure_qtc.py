#!/usr/bin/env python3
"""
KardiaMobile 6L QTc Measurement Tool — Vector Method
=====================================================

Measures QTc from a KardiaMobile 6L PDF export by parsing the PDF's vector
path data directly, recovering the original 300 Hz sampled ECG without
image processing.

Usage:
    python measure_qtc.py <path_to_kardia_pdf>

How it works
============

The Kardia PDF stores the ECG as ~14,000 tiny line segments per rhythm-strip
page in pure black (0 0 0 RG) — these are the actual sample points. Other
colors (light purple gridlines, gray text) are ignored. Parsing the PDF's
vector operators directly recovers the original ADC samples instead of
estimating the trace from rendered pixels.

Gains over a rasterization-based approach:
  - Native 300 Hz sample rate (3.33 ms per sample), no pixel quantization
  - All 6 leads recovered, not just Lead II
  - All 4 rhythm pages combined → ~30 beats per reading instead of ~10
  - QT measured in Leads II, III, aVF independently then averaged
  - No image processing, no smoothing, no gridline contamination
  - ECG voltage values in real mV (10 mm/mV calibration)

How the PDF encodes ECG
=======================

The Kardia PDF is produced by Skia and stores the ECG as PDF vector path
operators:

    x y m          # moveto - start a path at (x, y)
    x y l          # lineto - draw line to (x, y)
    S              # stroke - render the path

Before the ECG paths the stream sets stroke color to "0 0 0 RG" (black).
Every m/l operator after that point until end-of-page-content is an ECG
sample. Each lineto endpoint is one sample of the waveform.

PDF coordinates use points (1/72 inch). At 25 mm/s ECG paper speed:
  1 mm = 72/25.4 = 2.8346 PDF points
  1 PDF point = 14.111 ms
At 10 mm/mV amplitude:
  1 mV = 28.346 PDF points
  1 PDF point = 0.0353 mV

The 6 leads are arranged vertically with Y centers near 172, 260, 352,
441, 529, 621 PDF points (top-to-bottom: I, II, III, aVR, aVL, aVF).
The Y values stored in the PDF are within ±40 points of each strip's
center, so partitioning by Y separates leads cleanly.

Within each lead the X spacing of samples is 0.2362 PDF points, which
is exactly 1/300 second = 3.33 ms. This is AliveCor's published sample
rate, confirming we're recovering the actual ADC samples.

Measurement method
==================

For each rhythm page (pages 2-5 of the PDF):

  1. Extract ECG paths in the black color, split into 6 leads by Y-center.
  2. For each of Lead II, III, aVF: build a (time, voltage) signal.
  3. Detect R peaks via scipy find_peaks.
  4. For each beat, measure Q onset, T peak, T end:
       Q onset: walk backward from R until slope drops below upstroke threshold
       T peak: maximum amplitude in 150-450 ms after R
       T end: tangent method — steepest descending slope after T peak,
              extrapolated to the post-T baseline
  5. Compute QT = T_end - Q_onset for each beat.

Across all pages and three leads, we get ~50-90 raw beat measurements per
reading. Filter contaminated beats (density-based clustering robust to
outliers) and average. Compute QTc with Bazett and Fridericia formulas
from mean QT and mean RR.

Multi-lead averaging cancels lead-specific T-wave noise: when Lead II's
T end is hard to pin down on a particular beat, Lead III or aVF often
shows it more clearly.

Known gotchas
=============

  1. Lead Y-center positions are page-layout-specific. If AliveCor changes
     the Kardia PDF format, the Y centers will shift. To verify, look at
     the per-page Lead II plot in the diagnostic image — clean signal with
     P-QRS-T morphology means crop coords are right.

  2. Edge-of-page beats can have truncated segments (algorithm needs ~550 ms
     post-R window). The script drops the first and last beats per page per
     lead before further filtering.

  3. T waves are sometimes very small in Lead III. The script checks
     T-peak amplitude > 0.05 mV and skips beats where the T wave is too
     small to measure reliably in that specific lead.

  4. The post-T baseline reference is the median of the last 30 samples of
     the beat segment. If RR is very short (HR > 120) this baseline can
     overlap with the next P wave. Not an issue at normal HRs.

  5. RR/HR are computed from all Lead-II R-peak intervals across all pages,
     NOT from the filtered beat set. This avoids biasing HR by which beats
     happened to have measurable QT.

Dependencies: pypdf, numpy, scipy, matplotlib
Install: pip install pypdf numpy scipy matplotlib
"""

import sys
import os
import re
import argparse
import numpy as np
from scipy.signal import find_peaks
import pypdf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# === Calibration constants ===
PDF_POINTS_PER_MM = 72 / 25.4            # 2.8346
MS_PER_POINT = 40 / PDF_POINTS_PER_MM     # 14.111 ms per PDF point at 25 mm/s
MV_PER_POINT = 1 / (10 * PDF_POINTS_PER_MM)  # 0.0353 mV at 10 mm/mV
SAMPLE_RATE_HZ = 300                      # AliveCor native rate
SAMPLE_INTERVAL_MS = 1000 / SAMPLE_RATE_HZ  # 3.333 ms

# Y-center of each lead strip on a Kardia 6L rhythm page (PDF points)
LEAD_Y_CENTERS = {
    'I':   171.8,
    'II':  260.3,
    'III': 352.1,
    'aVR': 440.6,
    'aVL': 529.0,
    'aVF': 620.9,
}
LEAD_HALF_HEIGHT = 40  # samples within ±this distance of center belong to that lead

# Leads we use for QT measurement (positive QRS deflections, usable T waves)
QT_LEADS = ['II', 'III', 'aVF']

# Pages 2-5 of the Kardia PDF are the 4 rhythm strips (0-indexed: 1,2,3,4)
RHYTHM_PAGES = [1, 2, 3, 4]


def extract_ecg_samples(pdf_path, page_idx):
    """Parse one rhythm page, return dict mapping lead → (times_ms, voltages_mv)."""
    reader = pypdf.PdfReader(pdf_path)
    page = reader.pages[page_idx]
    content_str = page.get_contents().get_data().decode('latin-1', errors='replace')

    # ECG paths come after "0 0 0 RG" (black stroke). Everything before that
    # is gridlines (light purple) and text (various grays).
    match = re.search(r'\b0\s+0\s+0\s+RG\b', content_str)
    if not match:
        raise RuntimeError(f"No ECG color block found on page {page_idx + 1}")
    ecg_section = content_str[match.end():]

    # Parse all 'l' lineto endpoints (these are the sample points).
    samples = []
    for m in re.finditer(r'([-\d.]+)\s+([-\d.]+)\s+l\b', ecg_section):
        samples.append((float(m.group(1)), float(m.group(2))))
    if not samples:
        raise RuntimeError(f"No ECG samples on page {page_idx + 1}")
    samples = np.array(samples)

    leads = {}
    for lead_name, y_center in LEAD_Y_CENTERS.items():
        mask = np.abs(samples[:, 1] - y_center) < LEAD_HALF_HEIGHT
        pts = samples[mask]
        if len(pts) < 100:
            continue
        order = np.argsort(pts[:, 0])
        pts = pts[order]
        # Convert to time (ms from page start) and voltage (mV, baseline-corrected)
        times_ms = (pts[:, 0] - pts[0, 0]) * MS_PER_POINT
        baseline_y = np.median(pts[:, 1])
        # PDF Y axis points DOWN; flip sign so positive deflections go up
        voltages_mv = -(pts[:, 1] - baseline_y) * MV_PER_POINT
        leads[lead_name] = (times_ms, voltages_mv)
    return leads


def detect_r_peaks(voltages_mv, min_height_mv=0.5):
    """Indices of R peaks. Min 400ms between peaks → allows HR up to 150."""
    min_dist = int(400 / SAMPLE_INTERVAL_MS)
    r_peaks, _ = find_peaks(voltages_mv, height=min_height_mv, distance=min_dist)
    return r_peaks


def measure_qt_for_beat(times_ms, voltages_mv, r_idx, beat_idx):
    """Measure QT for one beat. Returns dict or None on failure."""
    n = len(times_ms)
    samples_per_ms = SAMPLE_RATE_HZ / 1000.0

    pre_r = int(135 * samples_per_ms)
    post_r = int(550 * samples_per_ms)
    start = max(0, r_idx - pre_r)
    end = min(n, r_idx + post_r)
    if end - start < int(300 * samples_per_ms):
        return None

    seg_t = times_ms[start:end]
    seg_v = voltages_mv[start:end]
    r_local = r_idx - start

    slopes = np.gradient(seg_v)

    # Q onset: walk backward from R until slope drops below rapid-upstroke threshold
    SLOPE_THRESHOLD = 0.02  # mV/sample = 6 mV/s
    q_local = r_local
    for j in range(r_local - 1, 0, -1):
        if slopes[j] < SLOPE_THRESHOLD:
            q_local = j
            break
    q_onset_ms = seg_t[q_local] - seg_t[0]

    # T peak: max in 150-450 ms after R
    tp_start = r_local + int(150 * samples_per_ms)
    tp_end = min(len(seg_v), r_local + int(450 * samples_per_ms))
    if tp_end <= tp_start:
        return None
    t_peak_local = tp_start + int(np.argmax(seg_v[tp_start:tp_end]))
    t_peak_ms = seg_t[t_peak_local] - seg_t[0]
    t_peak_amp = float(seg_v[t_peak_local])
    if t_peak_amp < 0.05:  # T wave too small to measure reliably in this lead
        return None

    # T end: tangent method
    desc_start = t_peak_local + 1
    desc_end = min(len(seg_v), t_peak_local + int(150 * samples_per_ms))
    if desc_end - desc_start < 5:
        return None
    desc_slopes = slopes[desc_start:desc_end]
    steepest_local = int(np.argmin(desc_slopes))
    steepest_idx = desc_start + steepest_local
    steepest_slope = float(desc_slopes[steepest_local])

    post_baseline = float(np.median(seg_v[-30:]))

    if steepest_slope < -0.005:
        # Extrapolate tangent line to post-T baseline
        t_end_idx_frac = steepest_idx + (post_baseline - seg_v[steepest_idx]) / steepest_slope
        if t_end_idx_frac < 0 or t_end_idx_frac >= len(seg_v):
            return None
        t_end_ms = float(t_end_idx_frac * SAMPLE_INTERVAL_MS)
    else:
        # T too flat — baseline crossing fallback
        t_end_local = t_peak_local
        for k in range(t_peak_local, len(seg_v)):
            if seg_v[k] <= post_baseline + 0.01:
                t_end_local = k
                break
        t_end_ms = seg_t[t_end_local] - seg_t[0]

    qt_ms = t_end_ms - q_onset_ms

    return {
        'beat_idx': beat_idx,
        'q_onset_ms': q_onset_ms,
        't_peak_ms': t_peak_ms,
        't_end_ms': t_end_ms,
        'qt_ms': qt_ms,
        't_peak_amplitude_mv': t_peak_amp,
        'seg_t': seg_t,
        'seg_v': seg_v,
        'q_idx_seg': q_local,
        't_peak_idx_seg': t_peak_local,
        'post_baseline': post_baseline,
    }


def filter_contaminated_beats(measurements):
    """Density-based outlier rejection. Keep beats clustered within ±20ms of densest cluster."""
    if len(measurements) < 3:
        return measurements
    qts = np.array([m['qt_ms'] for m in measurements])
    TOLERANCE_MS = 20
    densities = np.array([np.sum(np.abs(qts - q) <= TOLERANCE_MS) for q in qts])
    max_density = densities.max()
    candidates = np.where(densities == max_density)[0]
    # Tie-break by longer QT (contamination shortens, so trust longer cluster)
    anchor_idx = candidates[np.argmax(qts[candidates])]
    anchor_qt = qts[anchor_idx]
    keep_mask = np.abs(qts - anchor_qt) <= TOLERANCE_MS
    return [m for m, k in zip(measurements, keep_mask) if k]


def compute_qtc(qt_ms, rr_ms):
    rr_s = rr_ms / 1000.0
    return qt_ms / np.sqrt(rr_s), qt_ms / (rr_s ** (1.0 / 3.0))


def measure_qtc_from_pdf(pdf_path, plot_path=None):
    all_measurements = []
    all_rr_ms = []
    per_page_diag = []

    for page_idx in RHYTHM_PAGES:
        try:
            leads = extract_ecg_samples(pdf_path, page_idx)
        except RuntimeError:
            continue
        if 'II' not in leads:
            continue

        t_ii, v_ii = leads['II']
        r_peaks_ii = detect_r_peaks(v_ii, min_height_mv=0.5)
        if len(r_peaks_ii) < 3:
            continue

        # RR intervals from Lead II R-peak times
        rr_intervals = np.diff(t_ii[r_peaks_ii])
        all_rr_ms.extend(rr_intervals.tolist())

        page_meas_by_lead = {}
        for lead in QT_LEADS:
            if lead not in leads:
                continue
            t_l, v_l = leads[lead]
            min_h = 0.3 if lead == 'III' else 0.5
            r_peaks_l = detect_r_peaks(v_l, min_height_mv=min_h)

            lead_meas = []
            for i, r in enumerate(r_peaks_l):
                m = measure_qt_for_beat(t_l, v_l, r, i)
                if m is not None:
                    m['lead'] = lead
                    m['page'] = page_idx + 1
                    lead_meas.append(m)
            # Drop first and last beat per page per lead (edge effects)
            if len(lead_meas) >= 3:
                lead_meas = lead_meas[1:-1]
            elif len(lead_meas) >= 2:
                lead_meas = lead_meas[1:]
            page_meas_by_lead[lead] = lead_meas
            all_measurements.extend(lead_meas)

        per_page_diag.append({
            'page': page_idx + 1,
            'lead_ii_signal': (t_ii, v_ii),
            'r_peaks_ii': r_peaks_ii,
            'measurements_by_lead': page_meas_by_lead,
        })

    if len(all_measurements) < 5:
        raise RuntimeError(f"Only {len(all_measurements)} beats measurable — too few")
    if len(all_rr_ms) < 3:
        raise RuntimeError("Too few RR intervals")

    kept = filter_contaminated_beats(all_measurements)
    if len(kept) < 5:
        raise RuntimeError(f"Only {len(kept)} clean beats after filtering")

    qts = np.array([m['qt_ms'] for m in kept])
    rr_mean = float(np.mean(all_rr_ms))
    hr_bpm = 60000.0 / rr_mean
    qt_mean = float(qts.mean())
    qt_median = float(np.median(qts))
    qt_std = float(qts.std(ddof=1))

    qtc_b, qtc_f = compute_qtc(qt_mean, rr_mean)
    # Approximate standard error of QTc Bazett based on QT scatter and N
    qtc_b_se = (qt_std / np.sqrt(len(kept))) / np.sqrt(rr_mean / 1000)

    qts_by_lead = {}
    for lead in QT_LEADS:
        lead_qts = [m['qt_ms'] for m in kept if m['lead'] == lead]
        if lead_qts:
            qts_by_lead[lead] = {
                'n': len(lead_qts),
                'mean': float(np.mean(lead_qts)),
                'std': float(np.std(lead_qts, ddof=1)) if len(lead_qts) > 1 else 0.0,
            }

    if plot_path is not None:
        _plot_diagnostic(per_page_diag, kept, qts_by_lead, plot_path)

    return {
        'hr_bpm': hr_bpm,
        'rr_mean_ms': rr_mean,
        'rr_count': len(all_rr_ms),
        'qt_mean_ms': qt_mean,
        'qt_median_ms': qt_median,
        'qt_std_ms': qt_std,
        'qtc_bazett': float(qtc_b),
        'qtc_bazett_se': float(qtc_b_se),
        'qtc_fridericia': float(qtc_f),
        'beats_measured': len(all_measurements),
        'beats_kept': len(kept),
        'qts_by_lead': qts_by_lead,
        'diagnostic_image_path': plot_path,
    }


def _plot_diagnostic(per_page_diag, kept_measurements, qts_by_lead, output_path):
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(5, 4)

    kept_set = set((m['page'], m['lead'], m['beat_idx']) for m in kept_measurements)

    for i, pd_ in enumerate(per_page_diag[:4]):
        ax = fig.add_subplot(gs[i, :])
        t, v = pd_['lead_ii_signal']
        ax.plot(t, v, 'k-', linewidth=0.6)
        ax.axhline(0, color='gray', linewidth=0.3)
        ax.plot(t[pd_['r_peaks_ii']], v[pd_['r_peaks_ii']], 'rv', markersize=5)
        ii_meas = pd_['measurements_by_lead'].get('II', [])
        for m in ii_meas:
            color = 'blue' if (m['page'], m['lead'], m['beat_idx']) in kept_set else 'lightgray'
            t0 = m['seg_t'][0]
            qx = t0 + m['q_onset_ms']
            ex = t0 + m['t_end_ms']
            ax.plot([qx, ex], [-0.3, -0.3], '-', color=color, linewidth=2)
        ax.set_ylabel(f"Page {pd_['page']}\nLead II (mV)")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.5, 1.5)

    for i, lead in enumerate(QT_LEADS):
        ax = fig.add_subplot(gs[4, i])
        lead_qts = [m['qt_ms'] for m in kept_measurements if m['lead'] == lead]
        if lead_qts:
            ax.hist(lead_qts, bins=10, color='steelblue', edgecolor='black')
            ax.axvline(np.mean(lead_qts), color='red', linewidth=2,
                       label=f'Mean {np.mean(lead_qts):.0f} ms')
            ax.set_xlabel('QT (ms)')
            ax.set_title(f'Lead {lead}: n={len(lead_qts)}')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[4, 3])
    all_kept_qts = [m['qt_ms'] for m in kept_measurements]
    ax.hist(all_kept_qts, bins=15, color='darkorange', edgecolor='black')
    ax.axvline(np.mean(all_kept_qts), color='red', linewidth=2,
               label=f'Mean {np.mean(all_kept_qts):.0f} ms')
    ax.set_xlabel('QT (ms)')
    ax.set_title(f'All kept beats (n={len(all_kept_qts)})')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=100)
    plt.close()


def print_report(result):
    print()
    print("=" * 64)
    print("  KARDIA 6L QTc MEASUREMENT (vector method)")
    print("=" * 64)
    print(f"  Heart rate:           {result['hr_bpm']:.1f} BPM")
    print(f"  Mean RR:              {result['rr_mean_ms']:.0f} ms ({result['rr_count']} intervals)")
    print()
    print(f"  Beats measured:       {result['beats_measured']}")
    print(f"  Beats kept (clean):   {result['beats_kept']}")
    print(f"  Mean QT:              {result['qt_mean_ms']:.1f} ms")
    print(f"  Median QT:            {result['qt_median_ms']:.1f} ms")
    print(f"  QT std dev:           {result['qt_std_ms']:.1f} ms")
    print()
    print(f"  Per-lead breakdown:")
    for lead, info in result['qts_by_lead'].items():
        print(f"    Lead {lead:>4}: n={info['n']:>2}, mean QT = {info['mean']:.1f} ± {info['std']:.1f} ms")
    print()
    print(f"  QTc Bazett:           {result['qtc_bazett']:.1f} ms  (±{result['qtc_bazett_se']:.1f} SE)")
    print(f"  QTc Fridericia:       {result['qtc_fridericia']:.1f} ms")
    print()

    qtcb = result['qtc_bazett']
    if qtcb < 450:
        flag = "NORMAL"
    elif qtcb < 460:
        flag = "HIGH-NORMAL"
    elif qtcb < 470:
        flag = "BORDERLINE PROLONGED"
    else:
        flag = "PROLONGED"
    print(f"  Interpretation:       {flag}")
    print(f"  (Common reference thresholds: 450 ms caution, 470 ms prolonged.")
    print(f"   Not a medical device; consult your clinician.)")
    print()
    print(f"  Diagnostic plot:      {result['diagnostic_image_path']}")
    print("=" * 64)


def main():
    parser = argparse.ArgumentParser(description="Measure QTc from KardiaMobile 6L PDF (vector method)")
    parser.add_argument('pdf', help="Path to Kardia 6L PDF export")
    parser.add_argument('--workdir', default='/tmp/qtc_work',
                        help="Working directory for diagnostic image")
    args = parser.parse_args()

    if not os.path.exists(args.pdf):
        print(f"Error: PDF not found at {args.pdf}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.workdir, exist_ok=True)
    plot_path = os.path.join(args.workdir, 'qtc_diagnostic.png')
    result = measure_qtc_from_pdf(args.pdf, plot_path=plot_path)
    print_report(result)


if __name__ == '__main__':
    main()
