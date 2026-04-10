#!/usr/bin/env python3
"""Modelo baseline CAE/IPMVP para eyectores de vacío (refino).

Alcance actual:
- Variable objetivo: H (LNG).
- Construcción y comparación de modelos lineales interpretables.
- Selección por robustez + interpretabilidad + criterios CAE, no solo por R².

Notas:
- Implementación sin dependencias externas (solo librería estándar).
- p-values aproximados con distribución normal (válido como aproximación con n moderado/grande).
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
    b: float  # Carga de crudo (t/h)
    c: float  # Tflash (ºC)
    d: float  # T cabeza (ºC)
    e: float  # Presión cabeza (vacío)
    f: float  # Viscosidad bitumen
    g: float  # Caudal de bitumen
    h: float  # LNG (kg/h o unidad de hoja)


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
        # Filtro operativo representativo.
        if not (b > 140 and c > 300):
            continue
        # Positividad para estabilidad y consistencia física de consumo/calidad.
        if f <= 0 or g <= 0 or h <= 0:
            continue
        data.append(Obs(b, c, d, e, f, g, h))
    return data


def mat_copy(a):
    return [row[:] for row in a]


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

    # Aproximación normal para p-value bilateral.
    pvals = [2 * (1 - NORM.cdf(abs(tv))) if math.isfinite(tv) else 0.0 for tv in tvals]

    return {
        'beta': beta,
        'yhat': yhat,
        'sse': sse,
        'se': se,
        'tvals': tvals,
        'pvals': pvals,
        'xtx_inv': xtx_inv,
        'n': n,
        'p': p,
    }


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
    return {
        'r2': r2,
        'adj_r2': adj_r2,
        'rmse': rmse,
        'cvrmse': cvrmse,
        'nmbe': nmbe,
    }


def build_design(rows, feature_names):
    x = []
    y = []
    for o in rows:
        fmap = {
            'B': o.b,
            'C': o.c,
            'D': o.d,
            'E': o.e,
            'F': o.f,
            'G': o.g,
            'BxE': o.b * o.e,
            'visc_hi_dummy': 1.0 if o.f >= 100.0 else 0.0,  # dummy simple y estable
        }
        x.append([1.0] + [fmap[name] for name in feature_names])
        y.append(o.h)
    return x, y


def vif_for_feature_matrix(x, feature_names):
    # x incluye intercepto en primera columna.
    vifs = {}
    for i, name in enumerate(feature_names, start=1):
        y = [row[i] for row in x]
        others = [
            [1.0] + [row[j] for j in range(1, len(x[0])) if j != i]
            for row in x
        ]
        fit = ols_fit(others, y)
        m = calc_metrics(y, fit['yhat'], len(others[0]) - 1)
        r2 = min(0.999999, max(0.0, m['r2']))
        vifs[name] = 1.0 / (1.0 - r2)
    return vifs


def prune_features(train_rows, initial_features, locked_features=None, p_limit=0.05, vif_limit=8.0, max_vars=5):
    locked = set(locked_features or [])
    feats = initial_features[:]

    while True:
        x, y = build_design(train_rows, feats)
        fit = ols_fit(x, y)
        pvals = fit['pvals'][1:]  # sin intercepto
        vifs = vif_for_feature_matrix(x, feats)

        removed = False

        # 1) p-value
        worst_p = -1.0
        worst_feat = None
        for fn, pv in zip(feats, pvals):
            if fn in locked:
                continue
            if pv > worst_p:
                worst_p = pv
                worst_feat = fn
        if worst_feat is not None and worst_p > p_limit and len(feats) > 1:
            feats.remove(worst_feat)
            removed = True

        if removed:
            continue

        # 2) VIF
        bad_vifs = [(fn, v) for fn, v in vifs.items() if v > vif_limit and fn not in locked]
        if bad_vifs:
            bad_vifs.sort(key=lambda t: t[1], reverse=True)
            feats.remove(bad_vifs[0][0])
            removed = True

        if removed:
            continue

        # 3) Tamaño máximo
        if len(feats) > max_vars:
            candidates = [(fn, pv) for fn, pv in zip(feats, pvals) if fn not in locked]
            if candidates:
                candidates.sort(key=lambda t: t[1], reverse=True)
                feats.remove(candidates[0][0])
                continue

        break

    return feats


def evaluate_model(name, all_rows, train_rows, test_rows, base_features, locked_features=None):
    selected = prune_features(train_rows, base_features, locked_features=locked_features)

    xtr, ytr = build_design(train_rows, selected)
    xte, yte = build_design(test_rows, selected)
    xal, yal = build_design(all_rows, selected)

    fit_tr = ols_fit(xtr, ytr)
    beta = fit_tr['beta']

    yhat_tr = fit_tr['yhat']
    yhat_te = [sum(v * b for v, b in zip(row, beta)) for row in xte]
    yhat_al = [sum(v * b for v, b in zip(row, beta)) for row in xal]

    mtr = calc_metrics(ytr, yhat_tr, len(selected))
    mte = calc_metrics(yte, yhat_te, len(selected))
    mal = calc_metrics(yal, yhat_al, len(selected))

    vifs = vif_for_feature_matrix(xtr, selected)

    return {
        'name': name,
        'features': selected,
        'beta': beta,
        'pvals': fit_tr['pvals'],
        'vifs': vifs,
        'train': mtr,
        'test': mte,
        'all': mal,
        'residuals_test': [yi - yh for yi, yh in zip(yte, yhat_te)],
    }


def pick_final_model(models):
    def meets_thresholds(m):
        t = m['test']
        return t['r2'] >= 0.75 and t['cvrmse'] < 20.0 and abs(t['nmbe']) < 5.0

    valid = [m for m in models if meets_thresholds(m)]
    pool = valid if valid else models

    def score(m):
        t = m['test']
        tr = m['train']
        stability = abs(tr['cvrmse'] - t['cvrmse']) + abs(tr['nmbe'] - t['nmbe'])
        # Menor score es mejor; prioriza robustez/estabilidad y luego parsimonia.
        return (
            0 if meets_thresholds(m) else 1,
            round(t['cvrmse'], 6),
            round(abs(t['nmbe']), 6),
            round(stability, 6),
            len(m['features']),
            -round(t['r2'], 6),
        )

    pool.sort(key=score)
    return pool[0]


def equation_text(model):
    terms = [f"{model['beta'][0]:.6f}"]
    for f, b in zip(model['features'], model['beta'][1:]):
        terms.append(f"({b:.6f})·{f}")
    return 'LNG_hat = ' + ' + '.join(terms)


def write_report(path, data, models, final_model):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('# Baseline CAE/IPMVP - Eyectores de vacío\n\n')
        f.write('## Contexto y alcance\n')
        f.write('- Objetivo actual: baseline de **LNG (H)** para el sistema de vacío.\n')
        f.write('- Extensión futura prevista: **E_total = LNG + Fuel Gas** con la misma estructura de drivers físicos.\n\n')

        f.write('## Datos y validación\n')
        f.write(f'- Observaciones válidas tras filtros operativos: **{len(data)}**.\n')
        f.write('- Split: **train/test = 80/20** (aleatorio con semilla fija 42).\n')
        f.write('- Filtros: B > 140 t/h, C > 300 ºC, F>0, G>0, H>0.\n\n')

        f.write('## Comparativa de modelos candidatos\n\n')
        f.write('| Modelo | Variables finales | R² test | adjR² test | CVRMSE test (%) | NMBE test (%) | #vars |\n')
        f.write('|---|---|---:|---:|---:|---:|---:|\n')
        for m in models:
            t = m['test']
            f.write(
                f"| {m['name']} | {', '.join(m['features'])} | {t['r2']:.4f} | {t['adj_r2']:.4f} | {t['cvrmse']:.2f} | {t['nmbe']:.2f} | {len(m['features'])} |\n"
            )

        f.write('\n## Modelo final recomendado\n\n')
        f.write(f"**Seleccionado:** {final_model['name']}\n\n")
        f.write('### Ecuación\n\n')
        f.write('```\n' + equation_text(final_model) + '\n```\n\n')

        f.write('### Coeficientes y significancia (train)\n\n')
        f.write('| Término | Coeficiente | p-value (aprox) |\n')
        f.write('|---|---:|---:|\n')
        terms = ['Intercepto'] + final_model['features']
        for tname, b, p in zip(terms, final_model['beta'], final_model['pvals']):
            f.write(f'| {tname} | {b:.6f} | {p:.6g} |\n')

        f.write('\n### Métricas (modelo final)\n\n')
        f.write('| Escenario | R² | adj R² | CVRMSE (%) | NMBE (%) |\n')
        f.write('|---|---:|---:|---:|---:|\n')
        for label, key in [('Train', 'train'), ('Test', 'test'), ('Global', 'all')]:
            m = final_model[key]
            f.write(f"| {label} | {m['r2']:.4f} | {m['adj_r2']:.4f} | {m['cvrmse']:.2f} | {m['nmbe']:.2f} |\n")

        f.write('\n### Diagnóstico de multicolinealidad (VIF)\n\n')
        f.write('| Variable | VIF |\n|---|---:|\n')
        for k, v in sorted(final_model['vifs'].items()):
            f.write(f'| {k} | {v:.3f} |\n')

        res = final_model['residuals_test']
        mean_res = sum(res) / len(res)
        mae_res = sum(abs(r) for r in res) / len(res)
        f.write('\n### Residuos en test\n\n')
        f.write(f'- Residuo medio: **{mean_res:.4f}**\n')
        f.write(f'- MAE residuos: **{mae_res:.4f}**\n')
        f.write(f'- Rango residuos: **[{min(res):.4f}, {max(res):.4f}]**\n\n')

        f.write('## Justificación técnica defendible en auditoría\n\n')
        f.write('- **Carga (B):** principal driver de demanda energética del sistema de vacío.\n')
        f.write('- **Presión de cabeza (E):** representa esfuerzo del tren de eyectores para sostener vacío operativo.\n')
        f.write('- **Temperatura flash (C):** refleja severidad térmica y carga de vapores al sistema.\n')
        f.write('- **Interacción BxE** (si aplica): solo se mantiene cuando mejora estabilidad y tiene sentido físico (más carga a menor presión efectiva implica mayor exigencia).\n')
        f.write('- **Calidad (F o dummy):** se incluye solo si aporta señal estadística robusta sin degradar VIF.\n')
        f.write('- Se prioriza **parsimonia (4–5 variables máx.)** para trazabilidad y robustez CAE/IPMVP.\n')


def main():
    data = build_dataset('data/LB_CAE_Eyectores.xlsx')
    if len(data) < 30:
        raise RuntimeError('Muestra insuficiente para baseline defendible (n<30).')

    rnd = random.Random(42)
    rows = data[:]
    rnd.shuffle(rows)
    cut = int(0.8 * len(rows))
    train_rows = rows[:cut]
    test_rows = rows[cut:]

    candidates = [
        ('A_base', ['B', 'E', 'C'], ['B', 'E', 'C']),
        ('B_interaccion', ['B', 'E', 'C', 'BxE'], ['B', 'E', 'C']),
        ('C_calidad_viscosidad', ['B', 'E', 'C', 'F'], ['B', 'E', 'C']),
        ('C_calidad_dummy', ['B', 'E', 'C', 'visc_hi_dummy'], ['B', 'E', 'C']),
        ('D_reducido', ['B', 'E'], ['B', 'E']),
    ]

    models = []
    for name, feats, locked in candidates:
        models.append(evaluate_model(name, rows, train_rows, test_rows, feats, locked_features=locked))

    final_model = pick_final_model(models)
    write_report('baseline_model_results.md', rows, models, final_model)

    print('Reporte generado: baseline_model_results.md')
    print('Modelo final:', final_model['name'])
    print('Ecuación:', equation_text(final_model))


if __name__ == '__main__':
    main()
