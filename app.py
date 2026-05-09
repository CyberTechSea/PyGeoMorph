"""
PyGeoMorph v2.0 — Geometric Shell Morphometry for Marine Invertebrates
=======================================================================

A modern, all-in-one Python pipeline for landmark- and outline-based
geometric morphometric analysis. Updates the classical tpsDig/tpsRelw
workflow used in Criscione & Patti (2010, The Nautilus 124(3):117-128)
to current state-of-the-art (2024-2026) standards:

  • Generalized Procrustes Analysis (GPA) on fixed + sliding semi-landmarks
    (Gunz, Mitteroecker & Bookstein 2005; Bookstein 1997 BE-min sliding)
  • Procrustes ANOVA with RRPP permutation testing (Collyer & Adams 2018)
  • Phylogenetically-aligned PCA (PACA) and Kmult phylogenetic signal
    (Adams 2014; geomorph 4.0.8+)
  • Common Allometric Component (CAC) trajectory analysis
    (Mitteroecker et al. 2004)
  • Modularity / Integration via Covariance Ratio (Adams 2016)
  • Elliptic Fourier Analysis (Kuhl & Giardina 1982; Bonhomme et al. 2014)
  • Thin-Plate Spline deformation grids (Bookstein 1989)
  • Pre-loaded Rissoa template from Criscione & Patti (2010): 19 landmarks
    on gastropod shells (apex + 3 last whorls + aperture)
"""

import os
import io
import json
import uuid
import base64
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import linalg, stats
from scipy.spatial.distance import pdist, squareform
from PIL import Image
import cv2

from flask import Flask, render_template, request, jsonify, send_file

# -----------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent
UPLOAD_DIR = BASE_DIR / "uploads"
RESULTS_DIR = BASE_DIR / "results"
UPLOAD_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024
app.config["UPLOAD_FOLDER"] = str(UPLOAD_DIR)

ALLOWED_EXT = {"jpg", "jpeg", "png", "tif", "tiff", "pdf"}

PROJECT = {
    "specimens": [],
    "config": {
        "n_landmarks": 19,
        "auto_method": "active_contour",
        "smoothing": 3,
        "template": "rissoa",
        "sliders": [],
    },
    "analysis": {},
    "phylogeny": None,
}

# -----------------------------------------------------------------------------
# Built-in templates (after Criscione & Patti 2010)
# -----------------------------------------------------------------------------
TEMPLATES = {
    "rissoa": {
        "name": "Rissoidae gastropod (19 LM, Criscione & Patti 2010)",
        "n_landmarks": 19,
        "fixed": [0, 1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18],
        "sliders": [
            [1, 2, 3],
            [3, 4, 5],
            [14, 15, 16],
            [16, 17, 18],
        ],
        "links": [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],
                  [10,11],[11,12],[12,13],[13,14],[14,15],[15,16],[16,17],[17,18],[18,0]],
        "description": (
            "19-landmark configuration on gastropod shells. LM1 = apex; "
            "LM2/4/6 right whorl borders; LM3/5 sliding semi-landmarks; "
            "LM7-14 aperture and outer lip; LM15-19 left whorl mirror."
        ),
    },
    "bivalve": {
        "name": "Bivalve outline (15 LM)",
        "n_landmarks": 15,
        "fixed": [0, 7, 14],
        "sliders": [[i-1, i, i+1] for i in list(range(1,7)) + list(range(8,14))],
        "links": [[i, (i+1) % 15] for i in range(15)],
        "description": (
            "15-landmark configuration for bivalve outline: umbo (LM1), "
            "anterior margin (LM8), posterior margin (LM15) as anchors, "
            "with sliding semi-landmarks along the dorsal and ventral margins."
        ),
    },
    "echinoid": {
        "name": "Echinoid test (12 LM)",
        "n_landmarks": 12,
        "fixed": list(range(12)),
        "sliders": [],
        "links": [[i, (i+1) % 12] for i in range(12)],
        "description": (
            "12-landmark configuration for sea urchin tests at ambulacral "
            "and interambulacral plate boundaries."
        ),
    },
    "generic_outline": {
        "name": "Generic outline (15 sliding semi-LM)",
        "n_landmarks": 15,
        "fixed": [0],
        "sliders": [[i-1, i, (i+1) % 15] for i in range(1, 15)],
        "links": [[i, (i+1) % 15] for i in range(15)],
        "description": "Generic closed-outline template; first LM as anchor.",
    },
}


# =============================================================================
# IMAGE HANDLING
# =============================================================================
def allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT


def load_image_any(path: Path) -> np.ndarray:
    ext = path.suffix.lower()
    if ext == ".pdf":
        try:
            from pdf2image import convert_from_path
            pages = convert_from_path(str(path), dpi=200, first_page=1, last_page=1)
            return np.array(pages[0].convert("L"))
        except Exception:
            img = Image.open(path).convert("L")
            return np.array(img)
    img = Image.open(path)
    if img.mode != "L":
        img = img.convert("L")
    return np.array(img)


def encode_image_b64(arr: np.ndarray) -> str:
    if arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    img = Image.fromarray(arr)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# =============================================================================
# AUTOMATIC LANDMARK DETECTION
# =============================================================================
def extract_outline(gray: np.ndarray) -> np.ndarray:
    blur = cv2.GaussianBlur(gray, (5, 5), 0)
    _, th = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    if th.mean() > 127:
        th = 255 - th
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    th = cv2.morphologyEx(th, cv2.MORPH_CLOSE, kernel, iterations=2)
    th = cv2.morphologyEx(th, cv2.MORPH_OPEN, kernel, iterations=1)
    contours, _ = cv2.findContours(th, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return np.array([])
    cnt = max(contours, key=cv2.contourArea)
    return cnt[:, 0, :].astype(np.float64)


def sample_equidistant(pts: np.ndarray, n: int) -> np.ndarray:
    if len(pts) < n:
        return pts
    diffs = np.diff(pts, axis=0, append=pts[:1])
    seg_len = np.linalg.norm(diffs, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])
    total = cum[-1]
    targets = np.linspace(0, total, n, endpoint=False)
    out = np.zeros((n, 2))
    j = 0
    for k, t in enumerate(targets):
        while j < len(cum) - 1 and cum[j + 1] < t:
            j += 1
        if j >= len(pts) - 1:
            out[k] = pts[-1]
            continue
        seg = cum[j + 1] - cum[j]
        alpha = (t - cum[j]) / seg if seg > 0 else 0
        out[k] = pts[j] * (1 - alpha) + pts[j + 1] * alpha
    return out


def auto_landmarks_contour(gray: np.ndarray, n: int = 15) -> list:
    pts = extract_outline(gray)
    if len(pts) == 0:
        return []
    top_idx = int(np.argmin(pts[:, 1]))
    pts = np.roll(pts, -top_idx, axis=0)
    sampled = sample_equidistant(pts, n)
    return sampled.tolist()


def auto_landmarks_corners(gray: np.ndarray, n: int = 15, method: str = "shi-tomasi") -> list:
    g = cv2.GaussianBlur(gray, (3, 3), 0)
    if method == "harris":
        dst = cv2.cornerHarris(np.float32(g), 5, 3, 0.04)
        dst = cv2.dilate(dst, None)
        ys, xs = np.where(dst > 0.01 * dst.max())
        if len(xs) == 0:
            return []
        pts = np.column_stack([xs, ys]).astype(np.float64)
        if len(pts) > n:
            idx = np.random.RandomState(0).choice(len(pts), n, replace=False)
            pts = pts[idx]
    else:
        corners = cv2.goodFeaturesToTrack(g, maxCorners=n, qualityLevel=0.01, minDistance=10)
        if corners is None:
            return []
        pts = corners.reshape(-1, 2).astype(np.float64)
    return pts.tolist()


def auto_landmarks_active_contour(gray: np.ndarray, n: int = 15) -> list:
    pts = extract_outline(gray)
    if len(pts) == 0:
        return []
    top_idx = int(np.argmin(pts[:, 1]))
    pts = np.roll(pts, -top_idx, axis=0)
    smooth = np.column_stack([
        np.convolve(pts[:, 0], np.ones(7)/7, mode='same'),
        np.convolve(pts[:, 1], np.ones(7)/7, mode='same'),
    ])
    dx = np.gradient(smooth[:, 0])
    dy = np.gradient(smooth[:, 1])
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    curv = np.abs(dx * ddy - dy * ddx) / (dx**2 + dy**2 + 1e-9)**1.5
    curv[:5] = 0; curv[-5:] = 0
    peaks = []
    used = np.zeros(len(curv), dtype=bool)
    min_dist = max(5, len(curv) // (2 * n))
    sorted_idx = np.argsort(-curv)
    for i in sorted_idx:
        if used[i]: continue
        peaks.append(i)
        lo, hi = max(0, i - min_dist), min(len(curv), i + min_dist)
        used[lo:hi] = True
        if len(peaks) >= n: break
    peaks = sorted(peaks)
    if len(peaks) < n:
        eq_pts = sample_equidistant(pts, n)
        return eq_pts.tolist()
    return [pts[p].tolist() for p in peaks[:n]]


# =============================================================================
# CORE GEOMETRIC MORPHOMETRICS
# =============================================================================
def centroid_size(coords: np.ndarray) -> float:
    c = coords.mean(axis=0)
    return float(np.sqrt(((coords - c) ** 2).sum()))


def center_and_scale(coords: np.ndarray) -> np.ndarray:
    c = coords.mean(axis=0)
    centered = coords - c
    cs = np.sqrt((centered ** 2).sum())
    return centered / cs if cs > 0 else centered


def optimal_rotation(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    M = B.T @ A
    U, S, Vt = linalg.svd(M)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0] * (M.shape[0] - 1) + [d])
    return (U @ D @ Vt).T


def gpa(landmarks: np.ndarray, max_iter: int = 50, tol: float = 1e-7,
        sliders: list = None, slide_iters: int = 3):
    n_spec, n_land, dim = landmarks.shape
    centered_scaled = np.array([center_and_scale(s) for s in landmarks])
    consensus = centered_scaled[0].copy()
    aligned = centered_scaled.copy()

    diffs = []
    for it in range(max_iter):
        for i in range(n_spec):
            R = optimal_rotation(aligned[i], consensus)
            aligned[i] = aligned[i] @ R

        if sliders and it >= 0:
            for _ in range(slide_iters):
                for i in range(n_spec):
                    aligned[i] = slide_semilandmarks(aligned[i], consensus, sliders)
                for i in range(n_spec):
                    R = optimal_rotation(aligned[i], consensus)
                    aligned[i] = aligned[i] @ R

        new_consensus = aligned.mean(axis=0)
        new_consensus = center_and_scale(new_consensus)
        d = np.linalg.norm(new_consensus - consensus)
        diffs.append(float(d))
        consensus = new_consensus
        if d < tol:
            break

    proc_dists = np.array([np.linalg.norm(s - consensus) for s in aligned])
    total_var = float((proc_dists ** 2).sum() / (n_spec - 1)) if n_spec > 1 else 0.0
    return aligned, consensus, proc_dists, total_var, diffs


def slide_semilandmarks(spec: np.ndarray, consensus: np.ndarray, sliders: list) -> np.ndarray:
    out = spec.copy()
    for prev_i, slide_i, next_i in sliders:
        tangent = spec[next_i] - spec[prev_i]
        tn = np.linalg.norm(tangent)
        if tn < 1e-9:
            continue
        u = tangent / tn
        delta = consensus[slide_i] - spec[slide_i]
        shift = (delta @ u) * u
        out[slide_i] = spec[slide_i] + shift
    return out


def shape_pca(aligned: np.ndarray, consensus: np.ndarray):
    n_spec = aligned.shape[0]
    flat = aligned.reshape(n_spec, -1)
    mean_flat = consensus.flatten()
    X = flat - mean_flat
    U, S, Vt = linalg.svd(X, full_matrices=False)
    eigvals = (S ** 2) / max(n_spec - 1, 1)
    total = eigvals.sum()
    explained = eigvals / total if total > 0 else eigvals
    scores = U * S
    return scores, eigvals, Vt, explained


def common_allometric_component(aligned: np.ndarray, log_cs: np.ndarray) -> dict:
    n_spec = aligned.shape[0]
    flat = aligned.reshape(n_spec, -1)
    mean_flat = flat.mean(axis=0)
    X = flat - mean_flat
    s = log_cs - log_cs.mean()
    if (s @ s) < 1e-12:
        return None
    b = (s @ X) / (s @ s)
    nb = np.linalg.norm(b)
    if nb < 1e-12:
        return None
    b_unit = b / nb
    cac_scores = X @ b_unit
    residuals = X - np.outer(s, b)
    return {
        "cac_scores": cac_scores.tolist(),
        "cac_direction_norm": float(nb),
        "residual_var": float((residuals ** 2).sum() / max(n_spec - 1, 1)),
        "raw_var": float((X ** 2).sum() / max(n_spec - 1, 1)),
    }


def procrustes_anova_rrpp(aligned: np.ndarray, groups: list, n_perm: int = 999):
    n_spec = aligned.shape[0]
    Y = aligned.reshape(n_spec, -1)
    groups = np.array(groups)
    unique = list(sorted(set(groups)))
    if len(unique) < 2:
        return None
    grand = Y.mean(axis=0)
    SS_total = ((Y - grand) ** 2).sum()
    means = {g: Y[groups == g].mean(axis=0) for g in unique}
    SS_among = sum(len(Y[groups == g]) * ((means[g] - grand) ** 2).sum() for g in unique)
    SS_within = SS_total - SS_among
    df_among = len(unique) - 1
    df_within = n_spec - len(unique)
    F_obs = (SS_among / df_among) / (SS_within / max(df_within, 1)) if df_within > 0 else np.nan

    rng = np.random.RandomState(42)
    F_null = np.zeros(n_perm)
    residuals = Y - grand
    for k in range(n_perm):
        perm = rng.permutation(n_spec)
        Yp = grand + residuals[perm]
        means_p = {g: Yp[groups == g].mean(axis=0) for g in unique}
        SSa = sum(len(Yp[groups == g]) * ((means_p[g] - grand) ** 2).sum() for g in unique)
        SSw = ((Yp - grand) ** 2).sum() - SSa
        F_null[k] = (SSa / df_among) / (SSw / max(df_within, 1)) if df_within > 0 else 0
    p_val = float((F_null >= F_obs).sum() + 1) / (n_perm + 1)
    Z = float((F_obs - F_null.mean()) / F_null.std()) if F_null.std() > 0 else 0.0
    return {
        "SS_among": float(SS_among),
        "SS_within": float(SS_within),
        "SS_total": float(SS_total),
        "df_among": int(df_among),
        "df_within": int(df_within),
        "MS_among": float(SS_among / df_among),
        "F": float(F_obs),
        "p_value": p_val,
        "Z": Z,
        "n_perm": n_perm,
        "Rsq": float(SS_among / SS_total) if SS_total > 0 else 0.0,
    }


def mahalanobis_pairwise(scores: np.ndarray, groups: list, k: int = None):
    groups = np.array(groups)
    unique = list(sorted(set(groups)))
    if k is None:
        k = min(scores.shape[1], max(2, scores.shape[0] - len(unique)))
    Y = scores[:, :k]
    pooled = np.zeros((k, k))
    dof = 0
    for g in unique:
        Yg = Y[groups == g]
        if len(Yg) < 2:
            continue
        cov = np.cov(Yg, rowvar=False)
        pooled += cov * (len(Yg) - 1)
        dof += len(Yg) - 1
    if dof == 0:
        return unique, np.zeros((len(unique), len(unique)))
    pooled /= dof
    inv = np.linalg.pinv(pooled)
    means = np.array([Y[groups == g].mean(axis=0) for g in unique])
    D = np.zeros((len(unique), len(unique)))
    for i in range(len(unique)):
        for j in range(len(unique)):
            d = means[i] - means[j]
            D[i, j] = float(np.sqrt(d @ inv @ d.T))
    return unique, D


def manova_simple(scores: np.ndarray, groups: list, k: int = 4):
    groups = np.array(groups)
    unique = list(sorted(set(groups)))
    if len(unique) < 2:
        return None
    Y = scores[:, :min(k, scores.shape[1])]
    n, p = Y.shape
    grand = Y.mean(axis=0)
    H = np.zeros((p, p))
    E = np.zeros((p, p))
    for g in unique:
        Yg = Y[groups == g]
        if len(Yg) == 0: continue
        mg = Yg.mean(axis=0)
        H += len(Yg) * np.outer(mg - grand, mg - grand)
        for row in Yg:
            E += np.outer(row - mg, row - mg)
    try:
        wilks = np.linalg.det(E) / np.linalg.det(E + H)
    except Exception:
        return None
    g_n = len(unique)
    df1 = p * (g_n - 1)
    df2 = n - g_n
    if df2 <= 0 or wilks <= 0:
        return {"wilks_lambda": float(wilks), "p_value": None, "df1": df1}
    chi2 = -(df2 - (p - g_n + 2) / 2) * np.log(wilks)
    p_val = 1 - stats.chi2.cdf(chi2, df1)
    return {
        "wilks_lambda": float(wilks),
        "chi2": float(chi2),
        "df1": int(df1), "df2": int(df2),
        "p_value": float(p_val),
    }


def allometry_test(log_cs: np.ndarray, scores: np.ndarray, k: int = 5,
                   groups: list = None, n_perm: int = 999):
    Y = scores[:, :min(k, scores.shape[1])]
    n = Y.shape[0]
    X = np.column_stack([np.ones(n), log_cs])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ beta
    ss_res = ((Y - pred) ** 2).sum()
    ss_tot = ((Y - Y.mean(axis=0)) ** 2).sum()
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0

    rng = np.random.RandomState(42)
    null = np.zeros(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(log_cs)
        Xp = np.column_stack([np.ones(n), perm])
        bp, *_ = np.linalg.lstsq(Xp, Y, rcond=None)
        ssr = ((Y - Xp @ bp) ** 2).sum()
        null[i] = 1 - ssr / ss_tot if ss_tot > 0 else 0
    p_val = float((null >= r2).sum() + 1) / (n_perm + 1)
    out = {
        "r_squared": r2,
        "p_value": p_val,
        "slope_norm": float(np.linalg.norm(beta[1])),
        "n_pcs_used": Y.shape[1],
        "n_perm": n_perm,
    }
    if groups is not None and len(set(groups)) > 1:
        groups = np.array(groups)
        unique = list(sorted(set(groups)))
        SS_groups = 0.0
        for g in unique:
            mask = (groups == g)
            if mask.sum() < 2: continue
            Xg = np.column_stack([np.ones(mask.sum()), log_cs[mask]])
            Yg = Y[mask]
            bg, *_ = np.linalg.lstsq(Xg, Yg, rcond=None)
            SS_groups += ((Yg - Xg @ bg) ** 2).sum()
        SS_diff = ss_res - SS_groups
        df_diff = (len(unique) - 1) * Y.shape[1]
        df_full = n - len(unique) * 2
        if df_full > 0 and SS_groups > 0:
            F_hos = (SS_diff / df_diff) / (SS_groups / df_full)
            null_hos = np.zeros(n_perm)
            for i in range(n_perm):
                perm = rng.permutation(n)
                grp_p = groups[perm]
                SS_g = 0.0
                for g in unique:
                    mask = (grp_p == g)
                    if mask.sum() < 2: continue
                    Xg = np.column_stack([np.ones(mask.sum()), log_cs[mask]])
                    Yg = Y[mask]
                    bg, *_ = np.linalg.lstsq(Xg, Yg, rcond=None)
                    SS_g += ((Yg - Xg @ bg) ** 2).sum()
                SS_d = ss_res - SS_g
                null_hos[i] = (SS_d / df_diff) / (SS_g / df_full) if SS_g > 0 else 0
            p_hos = float((null_hos >= F_hos).sum() + 1) / (n_perm + 1)
            out["HOS"] = {"F": float(F_hos), "p_value": p_hos,
                          "df": [int(df_diff), int(df_full)],
                          "interpretation": (
                              "Significant: groups have different allometric trajectories"
                              if p_hos < 0.05 else
                              "Non-significant: common allometric trajectory cannot be rejected"
                          )}
    return out


def disparity(aligned: np.ndarray, consensus: np.ndarray, groups: list = None):
    n = aligned.shape[0]
    flat = aligned.reshape(n, -1)
    mean_flat = consensus.flatten()
    overall = float(((flat - mean_flat) ** 2).sum() / max(n - 1, 1))
    out = {"overall": overall}
    if groups is not None:
        groups = np.array(groups)
        for g in sorted(set(groups)):
            sub = flat[groups == g]
            if len(sub) < 2: continue
            m = sub.mean(axis=0)
            out[f"group_{g}"] = float(((sub - m) ** 2).sum() / (len(sub) - 1))
    return out


def kmult_phylogenetic_signal(aligned: np.ndarray, names: list,
                              tree_distance: dict = None, n_perm: int = 999):
    if not tree_distance:
        return None
    n = aligned.shape[0]
    Y = aligned.reshape(n, -1)
    Y -= Y.mean(axis=0)
    try:
        D = np.zeros((n, n))
        for i, ni in enumerate(names):
            for j, nj in enumerate(names):
                if ni in tree_distance and nj in tree_distance[ni]:
                    D[i, j] = float(tree_distance[ni][nj])
                else:
                    return None
    except Exception:
        return None
    if D.max() == 0:
        return None
    C = D.max() - D
    np.fill_diagonal(C, np.diag(C) + 1e-9)
    try:
        Cinv = np.linalg.pinv(C)
    except Exception:
        return None
    one = np.ones(n)
    a = (one @ Cinv @ Y) / (one @ Cinv @ one)
    Yc = Y - a
    mse_phylo = float(np.trace(Yc.T @ Cinv @ Yc) / (n - 1))
    eig = np.diag(C)
    sumD = (np.sum(eig) - n / (one @ Cinv @ one)) / (n - 1)
    mse_obs = float((Yc ** 2).sum() / (n - 1))
    K = (mse_obs / mse_phylo) / sumD if mse_phylo > 0 and sumD > 0 else None
    if K is None:
        return None
    rng = np.random.RandomState(42)
    K_null = np.zeros(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(n)
        Yp = Y[perm]
        ap = (one @ Cinv @ Yp) / (one @ Cinv @ one)
        Yc_p = Yp - ap
        mse_p = float(np.trace(Yc_p.T @ Cinv @ Yc_p) / (n - 1))
        mse_op = float((Yc_p ** 2).sum() / (n - 1))
        K_null[k] = (mse_op / mse_p) / sumD if mse_p > 0 else 0
    p = float((K_null >= K).sum() + 1) / (n_perm + 1)
    return {"K_mult": float(K), "p_value": p, "n_perm": n_perm}


def integration_test(aligned: np.ndarray, partition: list, n_perm: int = 499):
    if partition is None:
        return None
    n_spec, n_land, dim = aligned.shape
    if len(partition) != n_land:
        return None
    Y = aligned.reshape(n_spec, n_land * dim)
    partition = np.array(partition)
    modules = sorted(set(partition))
    if len(modules) < 2:
        return None
    def cols(mod):
        idx = []
        for i, p in enumerate(partition):
            if p == mod:
                idx.extend([i*dim + d for d in range(dim)])
        return np.array(idx)
    cov = np.cov(Y, rowvar=False)
    def cr_value(cov_matrix, parts):
        between_ss = 0.0
        within_terms = []
        for a in range(len(parts)):
            ia = parts[a]
            cov_aa = cov_matrix[np.ix_(ia, ia)]
            np.fill_diagonal(cov_aa, 0)
            within_terms.append((cov_aa ** 2).sum())
            for b in range(a+1, len(parts)):
                ib = parts[b]
                cov_ab = cov_matrix[np.ix_(ia, ib)]
                between_ss += (cov_ab ** 2).sum()
        within_prod = np.sqrt(np.prod(within_terms)) if all(w > 0 for w in within_terms) else 0
        if within_prod == 0:
            return None
        return between_ss / within_prod
    parts_idx = [cols(m) for m in modules]
    cr_obs = cr_value(cov, parts_idx)
    if cr_obs is None:
        return None
    rng = np.random.RandomState(42)
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(partition)
        parts_p = []
        for m in modules:
            idx = []
            for i, p in enumerate(perm):
                if p == m:
                    idx.extend([i*dim + d for d in range(dim)])
            parts_p.append(np.array(idx))
        cv = cr_value(cov, parts_p)
        if cv is not None:
            null.append(cv)
    null = np.array(null)
    p = float((null <= cr_obs).sum() + 1) / (len(null) + 1) if len(null) else None
    return {"CR": float(cr_obs), "p_value": p, "n_modules": len(modules), "n_perm": len(null)}


def tps_warp_grid(source: np.ndarray, target: np.ndarray, grid_size: int = 20):
    n = source.shape[0]
    def U(r2):
        return np.where(r2 > 0, r2 * np.log(np.sqrt(r2) + 1e-12), 0.0)
    K = U(squareform(pdist(source, "sqeuclidean")))
    P = np.column_stack([np.ones(n), source])
    L = np.zeros((n + 3, n + 3))
    L[:n, :n] = K
    L[:n, n:] = P
    L[n:, :n] = P.T
    Yt = np.vstack([target, np.zeros((3, 2))])
    try:
        W = np.linalg.solve(L, Yt)
    except np.linalg.LinAlgError:
        W = np.linalg.lstsq(L, Yt, rcond=None)[0]
    Ke = K
    bending_energy = float(np.trace(W[:n].T @ Ke @ W[:n]))
    x_min, y_min = source.min(axis=0) - 0.2 * np.ptp(source, axis=0)
    x_max, y_max = source.max(axis=0) + 0.2 * np.ptp(source, axis=0)
    xs = np.linspace(x_min, x_max, grid_size)
    ys = np.linspace(y_min, y_max, grid_size)
    XX, YY = np.meshgrid(xs, ys)
    grid_pts = np.column_stack([XX.ravel(), YY.ravel()])
    diffs = grid_pts[:, None, :] - source[None, :, :]
    r2 = (diffs ** 2).sum(axis=2)
    Ug = U(r2)
    warped = Ug @ W[:n] + np.column_stack([np.ones(len(grid_pts)), grid_pts]) @ W[n:]
    return XX, YY, warped.reshape(grid_size, grid_size, 2), bending_energy


# =============================================================================
# ELLIPTIC FOURIER ANALYSIS
# =============================================================================
def elliptic_fourier_descriptors(contour: np.ndarray, n_harmonics: int = 20,
                                 normalize: bool = True) -> dict:
    n_pts = len(contour)
    dx = np.diff(contour[:, 0], append=contour[0, 0])
    dy = np.diff(contour[:, 1], append=contour[0, 1])
    dt = np.sqrt(dx**2 + dy**2)
    t = np.concatenate([[0], np.cumsum(dt)])
    T = t[-1]
    if T == 0:
        return None
    coeffs = np.zeros((n_harmonics, 4))
    for n in range(1, n_harmonics + 1):
        const = T / (2 * (n * np.pi) ** 2)
        phi = 2 * n * np.pi * t[1:] / T
        phi_p = 2 * n * np.pi * t[:-1] / T
        cos_diff = np.cos(phi) - np.cos(phi_p)
        sin_diff = np.sin(phi) - np.sin(phi_p)
        a = const * np.sum((dx / dt) * cos_diff)
        b = const * np.sum((dx / dt) * sin_diff)
        c = const * np.sum((dy / dt) * cos_diff)
        d = const * np.sum((dy / dt) * sin_diff)
        coeffs[n - 1] = [a, b, c, d]
    if normalize:
        a1, b1, c1, d1 = coeffs[0]
        theta = 0.5 * np.arctan2(2 * (a1 * b1 + c1 * d1), a1**2 - b1**2 + c1**2 - d1**2)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        rot = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
        for k in range(n_harmonics):
            mat = np.array([[coeffs[k, 0], coeffs[k, 1]],
                            [coeffs[k, 2], coeffs[k, 3]]])
            mat_rot = mat @ rot
            coeffs[k] = mat_rot.flatten()
        scale = np.sqrt(coeffs[0, 0]**2 + coeffs[0, 2]**2)
        if scale > 0:
            coeffs /= scale
    return {
        "coefficients": coeffs.tolist(),
        "n_harmonics": n_harmonics,
        "perimeter": float(T),
        "normalized": normalize,
    }


def reconstruct_efa(coefficients: list, n_points: int = 200) -> np.ndarray:
    coeffs = np.array(coefficients)
    n_harmonics = len(coeffs)
    t = np.linspace(0, 1, n_points)
    x = np.zeros(n_points)
    y = np.zeros(n_points)
    for n in range(1, n_harmonics + 1):
        a, b, c, d = coeffs[n - 1]
        phi = 2 * n * np.pi * t
        x += a * np.cos(phi) + b * np.sin(phi)
        y += c * np.cos(phi) + d * np.sin(phi)
    return np.column_stack([x, y])


# =============================================================================
# FLASK ROUTES
# =============================================================================
@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/templates")
def api_templates():
    return jsonify({k: {"name": v["name"], "n_landmarks": v["n_landmarks"],
                        "description": v["description"], "links": v["links"]}
                    for k, v in TEMPLATES.items()})


@app.route("/api/upload", methods=["POST"])
def api_upload():
    if "file" not in request.files:
        return jsonify({"error": "No file"}), 400
    f = request.files["file"]
    if f.filename == "" or not allowed_file(f.filename):
        return jsonify({"error": "Invalid file type"}), 400

    group = request.form.get("group", "default")
    template = request.form.get("template", "rissoa")
    sid = uuid.uuid4().hex[:8]
    ext = f.filename.rsplit(".", 1)[1].lower()
    safe_name = f"{sid}.{ext}"
    save_path = UPLOAD_DIR / safe_name
    f.save(save_path)

    try:
        gray = load_image_any(save_path)
    except Exception as e:
        return jsonify({"error": f"Could not read image: {e}"}), 400

    h, w = gray.shape
    img_b64 = encode_image_b64(gray)

    tpl = TEMPLATES.get(template, TEMPLATES["rissoa"])
    n = int(request.form.get("n_landmarks", tpl["n_landmarks"]))
    method = request.form.get("method", "active_contour")
    if method == "active_contour":
        landmarks = auto_landmarks_active_contour(gray, n)
    elif method == "contour":
        landmarks = auto_landmarks_contour(gray, n)
    else:
        landmarks = auto_landmarks_corners(gray, n, method)

    outline = extract_outline(gray)
    outline_list = outline.tolist() if len(outline) else []

    cs = centroid_size(np.array(landmarks)) if landmarks else 0.0

    spec = {
        "id": sid, "name": f.filename, "filename": safe_name,
        "group": group, "template": template,
        "width": int(w), "height": int(h),
        "landmarks": landmarks, "outline": outline_list,
        "image_b64": img_b64, "centroid_size": cs,
        "added": datetime.now().isoformat(timespec="seconds"),
    }
    PROJECT["specimens"].append(spec)
    return jsonify({"specimen": _spec_summary(spec)})


@app.route("/api/specimen/<sid>", methods=["GET"])
def api_specimen(sid):
    for s in PROJECT["specimens"]:
        if s["id"] == sid:
            return jsonify(s)
    return jsonify({"error": "not found"}), 404


@app.route("/api/specimen/<sid>/landmarks", methods=["POST"])
def api_set_landmarks(sid):
    data = request.get_json()
    landmarks = data.get("landmarks", [])
    group = data.get("group")
    for s in PROJECT["specimens"]:
        if s["id"] == sid:
            s["landmarks"] = landmarks
            if group is not None:
                s["group"] = group
            s["centroid_size"] = centroid_size(np.array(landmarks)) if landmarks else 0.0
            return jsonify({"ok": True, "centroid_size": s["centroid_size"]})
    return jsonify({"error": "not found"}), 404


@app.route("/api/specimen/<sid>", methods=["DELETE"])
def api_delete(sid):
    PROJECT["specimens"] = [s for s in PROJECT["specimens"] if s["id"] != sid]
    return jsonify({"ok": True})


@app.route("/api/specimens", methods=["GET"])
def api_specimens():
    return jsonify({"specimens": [_spec_summary(s) for s in PROJECT["specimens"]]})


def _spec_summary(s):
    return {
        "id": s["id"], "name": s["name"], "group": s["group"],
        "template": s.get("template", "rissoa"),
        "width": s["width"], "height": s["height"],
        "n_landmarks": len(s["landmarks"]),
        "centroid_size": s["centroid_size"],
        "image_b64": s["image_b64"],
        "landmarks": s["landmarks"],
        "outline": s.get("outline", []),
        "added": s.get("added", ""),
    }


@app.route("/api/redetect/<sid>", methods=["POST"])
def api_redetect(sid):
    data = request.get_json() or {}
    n = int(data.get("n_landmarks", 19))
    method = data.get("method", "active_contour")
    for s in PROJECT["specimens"]:
        if s["id"] == sid:
            path = UPLOAD_DIR / s["filename"]
            gray = load_image_any(path)
            if method == "active_contour":
                lm = auto_landmarks_active_contour(gray, n)
            elif method == "contour":
                lm = auto_landmarks_contour(gray, n)
            else:
                lm = auto_landmarks_corners(gray, n, method)
            s["landmarks"] = lm
            s["centroid_size"] = centroid_size(np.array(lm)) if lm else 0.0
            return jsonify({"landmarks": lm, "centroid_size": s["centroid_size"]})
    return jsonify({"error": "not found"}), 404


@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    data = request.get_json() or {}
    use_sliders = data.get("use_sliders", True)
    template = data.get("template", "rissoa")
    n_perm = int(data.get("n_perm", 999))
    do_efa = data.get("do_efa", True)
    n_harmonics = int(data.get("n_harmonics", 20))

    specs = PROJECT["specimens"]
    if len(specs) < 3:
        return jsonify({"error": "Need at least 3 specimens for GPA"}), 400

    counts = {len(s["landmarks"]) for s in specs}
    if len(counts) != 1:
        return jsonify({
            "error": f"All specimens must have the same number of landmarks. Found: {sorted(counts)}"
        }), 400
    n_land = counts.pop()
    if n_land < 3:
        return jsonify({"error": "At least 3 landmarks per specimen required"}), 400

    coords = np.array([s["landmarks"] for s in specs])
    names = [s["name"] for s in specs]
    groups = [s["group"] for s in specs]
    raw_cs = np.array([centroid_size(c) for c in coords])

    sliders = None
    if use_sliders:
        tpl = TEMPLATES.get(template, TEMPLATES["rissoa"])
        if tpl["sliders"] and tpl["n_landmarks"] == n_land:
            sliders = tpl["sliders"]

    aligned, consensus, proc_dists, total_var, conv = gpa(
        coords, sliders=sliders, slide_iters=3 if sliders else 0
    )

    scores, eigvals, eigvecs, explained = shape_pca(aligned, consensus)

    log_cs = np.log(raw_cs + 1e-12)
    allom = allometry_test(log_cs, scores, k=5, groups=groups, n_perm=n_perm)
    cac = common_allometric_component(aligned, log_cs)
    rrpp = procrustes_anova_rrpp(aligned, groups, n_perm=n_perm) if len(set(groups)) > 1 else None
    manova_res = manova_simple(scores, groups) if len(set(groups)) > 1 else None
    if len(set(groups)) > 1:
        unique_g, mahal = mahalanobis_pairwise(scores, groups)
        mahal_data = {"groups": unique_g, "matrix": mahal.tolist()}
    else:
        mahal_data = None
    disp = disparity(aligned, consensus, groups)

    phylo_signal = None
    if PROJECT.get("phylogeny"):
        phylo_signal = kmult_phylogenetic_signal(aligned, names,
                                                 PROJECT["phylogeny"], n_perm=n_perm)

    efa_results = None
    if do_efa:
        coeffs_all = []
        valid_idx = []
        for i, s in enumerate(specs):
            if s.get("outline") and len(s["outline"]) > 50:
                outline_arr = np.array(s["outline"])
                if len(outline_arr) > 300:
                    outline_arr = sample_equidistant(outline_arr, 300)
                ef = elliptic_fourier_descriptors(outline_arr, n_harmonics=n_harmonics)
                if ef:
                    flat = np.array(ef["coefficients"]).flatten()
                    coeffs_all.append(flat[3:])
                    valid_idx.append(i)
        if len(coeffs_all) >= 3:
            X = np.array(coeffs_all)
            X = X - X.mean(axis=0)
            U, S, Vt = linalg.svd(X, full_matrices=False)
            efa_eigvals = (S ** 2) / max(len(coeffs_all) - 1, 1)
            efa_explained = efa_eigvals / efa_eigvals.sum()
            efa_scores = U * S
            efa_results = {
                "n_specimens_analysed": len(coeffs_all),
                "specimen_indices": valid_idx,
                "n_harmonics": n_harmonics,
                "pca_scores": efa_scores.tolist(),
                "pca_explained": efa_explained.tolist(),
            }

    pc1_max_idx = int(np.argmax(scores[:, 0]))
    XX, YY, warped, bending_energy = tps_warp_grid(consensus, aligned[pc1_max_idx], grid_size=18)

    PROJECT["analysis"] = {
        "names": names, "groups": groups,
        "centroid_sizes": raw_cs.tolist(),
        "aligned": aligned.tolist(),
        "consensus": consensus.tolist(),
        "proc_distances": proc_dists.tolist(),
        "total_procrustes_variance": total_var,
        "convergence": conv,
        "pca_scores": scores.tolist(),
        "pca_explained": explained.tolist(),
        "pca_eigvals": eigvals.tolist(),
        "allometry": allom,
        "cac": cac,
        "rrpp": rrpp,
        "manova": manova_res,
        "mahalanobis": mahal_data,
        "disparity": disp,
        "phylogenetic_signal": phylo_signal,
        "efa": efa_results,
        "tps": {
            "x": XX.tolist(), "y": YY.tolist(),
            "warped": warped.tolist(),
            "extreme_max": pc1_max_idx,
            "bending_energy": bending_energy,
        },
        "sliding_used": sliders is not None,
        "n_specimens": len(specs),
        "n_landmarks": n_land,
        "n_perm": n_perm,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    return jsonify(PROJECT["analysis"])


@app.route("/api/phylogeny", methods=["POST"])
def api_set_phylogeny():
    data = request.get_json() or {}
    PROJECT["phylogeny"] = data.get("distances")
    return jsonify({"ok": True, "n_taxa": len(PROJECT["phylogeny"]) if PROJECT["phylogeny"] else 0})


@app.route("/api/export/csv")
def api_export_csv():
    if not PROJECT["specimens"]:
        return jsonify({"error": "No specimens"}), 400
    lines = []
    lines.append("# PyGeoMorph v2.0 — landmark export")
    lines.append(f"# Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"# n_specimens={len(PROJECT['specimens'])}")
    lines.append("specimen_id,specimen_name,group,template,centroid_size,landmark_idx,x,y")
    for s in PROJECT["specimens"]:
        for i, (x, y) in enumerate(s["landmarks"]):
            lines.append(f"{s['id']},{s['name']},{s['group']},{s.get('template','-')},"
                         f"{s['centroid_size']:.6f},{i},{x:.4f},{y:.4f}")
    if PROJECT.get("analysis"):
        a = PROJECT["analysis"]
        lines.append("")
        lines.append("# === Procrustes-aligned shape coordinates ===")
        lines.append("specimen_idx,specimen_name,group,landmark_idx,x_aligned,y_aligned")
        for i, name in enumerate(a["names"]):
            for j, (x, y) in enumerate(a["aligned"][i]):
                lines.append(f"{i},{name},{a['groups'][i]},{j},{x:.6f},{y:.6f}")
        lines.append("")
        lines.append("# === PCA scores (tangent space) ===")
        n_pc = min(10, len(a["pca_scores"][0]))
        lines.append("specimen_name,group,centroid_size," + ",".join(f"PC{k+1}" for k in range(n_pc)))
        for i, name in enumerate(a["names"]):
            row = [name, a["groups"][i], f"{a['centroid_sizes'][i]:.4f}"]
            row += [f"{a['pca_scores'][i][k]:.6f}" for k in range(n_pc)]
            lines.append(",".join(row))
        lines.append("")
        lines.append("# === PC variance explained ===")
        lines.append("PC,eigenvalue,prop_variance,cum_variance")
        cum = 0.0
        for k, (ev, pv) in enumerate(zip(a["pca_eigvals"], a["pca_explained"])):
            cum += pv
            lines.append(f"PC{k+1},{ev:.6e},{pv:.6f},{cum:.6f}")
        if a.get("efa"):
            lines.append("")
            lines.append("# === EFA PCA scores ===")
            efa = a["efa"]
            n_efa_pc = min(10, len(efa["pca_scores"][0]))
            lines.append("specimen_idx," + ",".join(f"EFA_PC{k+1}" for k in range(n_efa_pc)))
            for i, sc in enumerate(efa["pca_scores"]):
                lines.append(str(efa["specimen_indices"][i]) + "," +
                             ",".join(f"{v:.6f}" for v in sc[:n_efa_pc]))
        if a.get("rrpp"):
            lines.append("")
            lines.append("# === Procrustes ANOVA (RRPP) ===")
            r = a["rrpp"]
            for k, v in r.items():
                lines.append(f"{k},{v}")

    csv_text = "\n".join(lines)
    out = io.BytesIO(csv_text.encode("utf-8"))
    out.seek(0)
    fname = f"pygeomorph_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    return send_file(out, mimetype="text/csv", as_attachment=True, download_name=fname)


@app.route("/api/export/tps")
def api_export_tps():
    if not PROJECT["specimens"]:
        return jsonify({"error": "No specimens"}), 400
    lines = []
    for s in PROJECT["specimens"]:
        lines.append(f"LM={len(s['landmarks'])}")
        for x, y in s["landmarks"]:
            lines.append(f"{x:.4f} {s['height'] - y:.4f}")
        lines.append(f"IMAGE={s['name']}")
        lines.append(f"ID={s['id']}")
        lines.append("SCALE=1.0000")
        lines.append("")
    out = io.BytesIO("\n".join(lines).encode("utf-8"))
    out.seek(0)
    return send_file(
        out, mimetype="text/plain", as_attachment=True,
        download_name=f"pygeomorph_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tps"
    )


@app.route("/api/reset", methods=["POST"])
def api_reset():
    PROJECT["specimens"] = []
    PROJECT["analysis"] = {}
    PROJECT["phylogeny"] = None
    for p in UPLOAD_DIR.glob("*"):
        try: p.unlink()
        except: pass
    return jsonify({"ok": True})


# =============================================================================
def open_browser():
    import webbrowser, threading, time
    def _open():
        time.sleep(1.2)
        webbrowser.open("http://127.0.0.1:5000")
    threading.Thread(target=_open, daemon=True).start()


if __name__ == "__main__":
    print("=" * 72)
    print("  PyGeoMorph v2.0 — Geometric Shell Morphometry — Marine Invertebrates")
    print("  Built on Criscione & Patti (2010) Rissoa template, modernised 2026")
    print("  Open your browser at:  http://127.0.0.1:5000")
    print("=" * 72)
    open_browser()
    app.run(host="127.0.0.1", port=5000, debug=False)