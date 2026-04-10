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
from statistics import NormalDist, correlation

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


def betacf(a, b, x):
    """Continued fraction for incomplete beta function (Numerical Recipes)."""
    max_iter = 200
    eps = 3.0e-14
    fpmin = 1.0e-300

    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d

    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break

    return h


def betainc_reg(a, b, x):
    """Regularized incomplete beta I_x(a,b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_bt = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    bt = math.exp(ln_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b


def student_t_cdf(t, dof):
    if dof <= 0:
        return 0.5
    if t == 0:
        return 0.5
    x = dof / (dof + t * t)
    ib = betainc_reg(dof / 2.0, 0.5, x)
    if t > 0:
        return 1.0 - 0.5 * ib
    return 0.5 * ib


def f_cdf(x, dfn, dfd):
    if x <= 0.0:
        return 0.0
    if dfn <= 0 or dfd <= 0:
        return 0.0
    z = (dfn * x) / (dfn * x + dfd)
    return betainc_reg(dfn / 2.0, dfd / 2.0, z)


def student_t_ppf(p, dof):
    # Inversión numérica simple (búsqueda binaria), suficiente para IC al 95%.
    lo, hi = -20.0, 20.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        cmid = student_t_cdf(mid, dof)
        if cmid < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


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
    pvals = [2 * (1 - student_t_cdf(abs(tv), dof)) if math.isfinite(tv) else 0.0 for tv in tvals]

    yavg = sum(y) / n
    sst = sum((yi - yavg) ** 2 for yi in y)
    ssr = max(0.0, sst - sse)
    dfn = max(1, p - 1)
    dfd = dof
    msr = ssr / dfn
    mse = sse / dfd
    f_stat = (msr / mse) if mse > 0 else float('inf')
    f_pval = 1.0 - f_cdf(f_stat, dfn, dfd) if math.isfinite(f_stat) else 0.0

    t_crit = student_t_ppf(0.975, dof)
    ci95 = [(b - t_crit * s, b + t_crit * s) for b, s in zip(beta, se)]

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
        'dof': dof,
        'f_stat': f_stat,
        'f_pval': f_pval,
        'sse': sse,
        'sst': sst,
        'ssr': ssr,
        'ci95': ci95,
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
    corr_y_yhat_train = correlation(ytr, yhat_tr) if len(ytr) > 1 else 0.0

    return {
        'name': name,
        'features': selected,
        'beta': beta,
        'pvals': fit_tr['pvals'],
        'vifs': vifs,
        'f_stat': fit_tr['f_stat'],
        'f_pval': fit_tr['f_pval'],
        'tvals': fit_tr['tvals'],
        'se': fit_tr['se'],
        'ci95': fit_tr['ci95'],
        'corr_y_yhat_train': corr_y_yhat_train,
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

        f.write('\n### Pruebas estadísticas de correlación y significancia global (train)\n\n')
        f.write(f"- **Correlación Pearson (Y vs Ŷ):** r = **{final_model['corr_y_yhat_train']:.4f}**\n")
        f.write(f"- **ANOVA F global:** F = **{final_model['f_stat']:.4f}**, p-value = **{final_model['f_pval']:.6g}**\n")
        f.write('- Interpretación: p-value global < 0.05 respalda que el modelo explica varianza de LNG mejor que un modelo sin predictores.\n\n')

        f.write('### Coeficientes con error estándar, t-stat e IC95% (train)\n\n')
        f.write('| Término | Coeficiente | Error estándar | t-stat | p-value | IC95% inferior | IC95% superior |\n')
        f.write('|---|---:|---:|---:|---:|---:|---:|\n')
        for tname, b, se, tv, pv, ci in zip(
            terms,
            final_model['beta'],
            final_model['se'],
            final_model['tvals'],
            final_model['pvals'],
            final_model['ci95'],
        ):
            f.write(f'| {tname} | {b:.6f} | {se:.6f} | {tv:.4f} | {pv:.6g} | {ci[0]:.6f} | {ci[1]:.6f} |\n')

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
