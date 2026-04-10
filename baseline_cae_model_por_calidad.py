#!/usr/bin/env python3
"""Modelo baseline CAE/IPMVP segregado por calidades de betún (F).

Construye líneas base separadas de LNG (H) para cada segmento de calidad,
definido por terciles de viscosidad (F).
"""

import math
import random
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from statistics import NormalDist

NS = {'a': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
NORM = NormalDist()


@dataclass
class Obs:
    b: float
    c: float
    d: float
    e: float
    f: float
    g: float
    h: float


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


def load_shared_strings(zf):
    if 'xl/sharedStrings.xml' not in zf.namelist():
        return []
    root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
    out = []
    for si in root.findall('a:si', NS):
        out.append(''.join((t.text or '') for t in si.findall('.//a:t', NS)))
    return out


def read_sheet(zf, sheet_path, shared_strings):
    root = ET.fromstring(zf.read(sheet_path))
    rows = []
    for row in root.findall('.//a:sheetData/a:row', NS):
        rec = {}
        for cell in row.findall('a:c', NS):
            ref = cell.get('r', '')
            col = ''.join(ch for ch in ref if ch.isalpha())
            v = cell.find('a:v', NS)
            if v is None:
                val = ''
            elif cell.get('t') == 's':
                val = shared_strings[int(v.text)]
            else:
                val = v.text
            rec[col] = val
        rows.append(rec)
    return rows


def build_dataset(path):
    with zipfile.ZipFile(path) as zf:
        ss = load_shared_strings(zf)
        corr = read_sheet(zf, 'xl/worksheets/sheet1.xml', ss)

    data = []
    for r in corr[2:]:
        vals = [fnum(r.get(k)) for k in ['B', 'C', 'D', 'E', 'F', 'G', 'H']]
        if any(v is None for v in vals):
            continue
        b, c, d, e, f, g, h = vals
        if any(math.isnan(v) or math.isinf(v) for v in vals):
            continue
        if not (b > 140 and c > 300):
            continue
        if f <= 0 or g <= 0 or h <= 0:
            continue
        data.append(Obs(b, c, d, e, f, g, h))
    return data


def percentile(sorted_vals, p):
    if not sorted_vals:
        raise ValueError('Sin datos para percentil')
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = (len(sorted_vals) - 1) * p
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def define_quality_segments(data):
    fvals = sorted(o.f for o in data)
    q33 = percentile(fvals, 1.0 / 3.0)
    q66 = percentile(fvals, 2.0 / 3.0)

    out = {
        'Q1_baja': [o for o in data if o.f <= q33],
        'Q2_media': [o for o in data if q33 < o.f <= q66],
        'Q3_alta': [o for o in data if o.f > q66],
    }
    return q33, q66, out


def mat_inv(a):
    n = len(a)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(a)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1e-12:
            raise ValueError('Matriz singular en inversión')
        aug[col], aug[pivot] = aug[pivot], aug[col]

        div = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= div

        for i in range(n):
            if i == col:
                continue
            fac = aug[i][col]
            if fac == 0:
                continue
            for j in range(2 * n):
                aug[i][j] -= fac * aug[col][j]

    return [row[n:] for row in aug]


def mat_vec_mul(a, x):
    return [sum(ai * xi for ai, xi in zip(row, x)) for row in a]


def ols_fit(x, y):
    n = len(x)
    p = len(x[0])
    xtx = [[0.0] * p for _ in range(p)]
    xty = [0.0] * p
    for row, yi in zip(x, y):
        for i in range(p):
            xty[i] += row[i] * yi
            for j in range(i, p):
                xtx[i][j] += row[i] * row[j]
    for i in range(p):
        for j in range(i):
            xtx[i][j] = xtx[j][i]

    xtx_inv = mat_inv(xtx)
    beta = mat_vec_mul(xtx_inv, xty)
    yhat = [sum(v * b for v, b in zip(row, beta)) for row in x]

    sse = sum((yi - yh) ** 2 for yi, yh in zip(y, yhat))
    dof = max(1, n - p)
    sigma2 = sse / dof

    se = [math.sqrt(max(0.0, sigma2 * xtx_inv[i][i])) for i in range(p)]
    tvals = [beta[i] / se[i] if se[i] > 0 else float('inf') for i in range(p)]
    pvals = [2 * (1 - NORM.cdf(abs(tv))) if math.isfinite(tv) else 0.0 for tv in tvals]

    return {'beta': beta, 'yhat': yhat, 'pvals': pvals}


def calc_metrics(y, yhat, p_params):
    n = len(y)
    yavg = sum(y) / n
    sst = sum((yi - yavg) ** 2 for yi in y)
    sse = sum((yi - yh) ** 2 for yi, yh in zip(y, yhat))
    r2 = 1 - (sse / sst if sst else 0.0)
    adj_r2 = 1 - (1 - r2) * (n - 1) / max(1, n - p_params - 1)
    rmse = math.sqrt(sse / n)
    cvrmse = (rmse / yavg * 100) if yavg else float('inf')
    nmbe = (sum(yi - yh for yi, yh in zip(y, yhat)) / max(1, n - p_params) / yavg * 100) if yavg else float('inf')
    return {'r2': r2, 'adj_r2': adj_r2, 'cvrmse': cvrmse, 'nmbe': nmbe}


def build_design(rows, feature_names):
    x, y = [], []
    for o in rows:
        fmap = {
            'B': o.b,
            'C': o.c,
            'D': o.d,
            'E': o.e,
            'G': o.g,
            'BxE': o.b * o.e,
            'CxE': o.c * o.e,
        }
        x.append([1.0] + [fmap[n] for n in feature_names])
        y.append(o.h)
    return x, y


def evaluate_model(rows, train_rows, test_rows, name, features):
    xtr, ytr = build_design(train_rows, features)
    xte, yte = build_design(test_rows, features)
    xal, yal = build_design(rows, features)

    fit = ols_fit(xtr, ytr)
    beta = fit['beta']

    yhat_tr = fit['yhat']
    yhat_te = [sum(v * b for v, b in zip(row, beta)) for row in xte]
    yhat_al = [sum(v * b for v, b in zip(row, beta)) for row in xal]

    return {
        'name': name,
        'features': features,
        'beta': beta,
        'pvals': fit['pvals'],
        'train': calc_metrics(ytr, yhat_tr, len(features)),
        'test': calc_metrics(yte, yhat_te, len(features)),
        'all': calc_metrics(yal, yhat_al, len(features)),
    }


def pick_best(models):
    def score(m):
        t = m['test']
        return (round(t['cvrmse'], 6), round(abs(t['nmbe']), 6), -round(t['r2'], 6), len(m['features']))

    return sorted(models, key=score)[0]


def equation_text(model):
    terms = [f"{model['beta'][0]:.6f}"]
    for f, b in zip(model['features'], model['beta'][1:]):
        terms.append(f"({b:.6f})·{f}")
    return 'LNG_hat = ' + ' + '.join(terms)


def write_report(path, n_total, q33, q66, outputs):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('# Baseline CAE/IPMVP por calidad de betún (F)\n\n')
        f.write('## Metodología de segregación\n')
        f.write(f'- Dataset válido global: **{n_total}** observaciones.\n')
        f.write('- Segmentación por terciles de viscosidad de betún (F):\n')
        f.write(f'  - **Q1 (baja): F ≤ {q33:.3f}**\n')
        f.write(f'  - **Q2 (media): {q33:.3f} < F ≤ {q66:.3f}**\n')
        f.write(f'  - **Q3 (alta): F > {q66:.3f}**\n\n')
        f.write('- En cada segmento se ajusta y valida una línea base independiente (split 80/20, semilla 42).\n\n')

        for out in outputs:
            name = out['segment']
            n_obs = out['n']
            final_model = out['best']
            candidates = out['models']

            f.write(f'## {name}\n\n')
            f.write(f'- Observaciones válidas en segmento: **{n_obs}**\n\n')

            f.write('| Modelo | Variables | R² test | CVRMSE test (%) | NMBE test (%) |\n')
            f.write('|---|---|---:|---:|---:|\n')
            for m in candidates:
                t = m['test']
                f.write(f"| {m['name']} | {', '.join(m['features'])} | {t['r2']:.4f} | {t['cvrmse']:.2f} | {t['nmbe']:.2f} |\n")

            f.write('\n')
            f.write(f"**Modelo base recomendado ({name}):** {final_model['name']}\n\n")
            f.write('```\n' + equation_text(final_model) + '\n```\n\n')

            f.write('| Escenario | R² | adj R² | CVRMSE (%) | NMBE (%) |\n')
            f.write('|---|---:|---:|---:|---:|\n')
            for label, key in [('Train', 'train'), ('Test', 'test'), ('Global', 'all')]:
                mm = final_model[key]
                f.write(f"| {label} | {mm['r2']:.4f} | {mm['adj_r2']:.4f} | {mm['cvrmse']:.2f} | {mm['nmbe']:.2f} |\n")
            f.write('\n')


def main():
    data = build_dataset('data/LB_CAE_Eyectores.xlsx')
    if len(data) < 90:
        raise RuntimeError('Muestra global insuficiente para segmentar en 3 calidades.')

    q33, q66, segments = define_quality_segments(data)
    rnd = random.Random(42)

    outputs = []
    for seg_name, seg_rows in segments.items():
        if len(seg_rows) < 60:
            raise RuntimeError(f'Segmento {seg_name} insuficiente (n={len(seg_rows)}).')

        rows = seg_rows[:]
        rnd.shuffle(rows)
        cut = int(0.8 * len(rows))
        train_rows = rows[:cut]
        test_rows = rows[cut:]

        candidates = [
            ('A_base', ['B', 'E', 'C']),
            ('B_interaccion_BE', ['B', 'E', 'C', 'BxE']),
            ('C_con_caudal_G', ['B', 'E', 'C', 'G']),
            ('D_reducido', ['B', 'E']),
        ]

        models = [evaluate_model(rows, train_rows, test_rows, name, feats) for name, feats in candidates]
        best = pick_best(models)

        outputs.append({'segment': seg_name, 'n': len(rows), 'models': models, 'best': best})

    write_report('baseline_model_results_por_calidad.md', len(data), q33, q66, outputs)
    print('Reporte generado: baseline_model_results_por_calidad.md')


if __name__ == '__main__':
    main()
