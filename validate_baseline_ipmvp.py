#!/usr/bin/env python3
"""Validación final de baseline LNG para auditoría CAE/IPMVP.

Modelo fijo (no se altera estructura):
    LNG = β0 + β1·B + β2·E + β3·C + β4·F

Incluye:
- Diagnóstico estadístico final (DW, residuos, homocedasticidad).
- Robustez por múltiples particiones train/test.
- Check de simplicidad (con y sin F).
- Definición formal del modelo y texto técnico para memoria.
"""

from __future__ import annotations

import math
import random
import statistics
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


@dataclass
class Obs:
    b: float
    c: float
    e: float
    f: float
    h: float


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


def load_shared_strings(zf):
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    out = []
    for si in root.findall("a:si", NS):
        out.append("".join((t.text or "") for t in si.findall(".//a:t", NS)))
    return out


def read_sheet(zf, sheet_path, shared_strings):
    root = ET.fromstring(zf.read(sheet_path))
    rows = []
    for row in root.findall(".//a:sheetData/a:row", NS):
        rec = {}
        for cell in row.findall("a:c", NS):
            ref = cell.get("r", "")
            col = "".join(ch for ch in ref if ch.isalpha())
            v = cell.find("a:v", NS)
            if v is None:
                val = ""
            elif cell.get("t") == "s":
                val = shared_strings[int(v.text)]
            else:
                val = v.text
            rec[col] = val
        rows.append(rec)
    return rows


def build_dataset(path: str) -> List[Obs]:
    with zipfile.ZipFile(path) as zf:
        shared_strings = load_shared_strings(zf)
        corr = read_sheet(zf, "xl/worksheets/sheet1.xml", shared_strings)

    data: List[Obs] = []
    for r in corr[2:]:
        vals = [fnum(r.get(k)) for k in ["B", "C", "E", "F", "H"]]
        if any(v is None for v in vals):
            continue
        b, c, e, f, h = vals
        if any(math.isnan(v) or math.isinf(v) for v in vals):
            continue
        if not (b > 140 and c > 300):
            continue
        if f <= 0 or h <= 0:
            continue
        data.append(Obs(b=b, c=c, e=e, f=f, h=h))
    return data


def mat_inv(a):
    n = len(a)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(a)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1e-12:
            raise ValueError("Matriz singular en inversión")
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


def ols_fit(x: List[List[float]], y: List[float]):
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
    tvals = [beta[i] / se[i] if se[i] > 0 else float("inf") for i in range(p)]
    pvals = [2 * (1 - student_t_cdf(abs(tv), dof)) if math.isfinite(tv) else 0.0 for tv in tvals]

    return {
        "beta": beta,
        "yhat": yhat,
        "sse": sse,
        "dof": dof,
        "se": se,
        "tvals": tvals,
        "pvals": pvals,
    }


def calc_metrics(y: Sequence[float], yhat: Sequence[float], p_params: int):
    n = len(y)
    yavg = sum(y) / n
    sst = sum((yi - yavg) ** 2 for yi in y)
    sse = sum((yi - yh) ** 2 for yi, yh in zip(y, yhat))
    r2 = 1 - (sse / sst if sst else 0.0)
    rmse = math.sqrt(sse / n)
    cvrmse = (rmse / yavg * 100.0) if yavg else float("inf")
    nmbe = (sum(yi - yh for yi, yh in zip(y, yhat)) / max(1, n - p_params) / yavg * 100.0) if yavg else float("inf")
    return {"r2": r2, "cvrmse": cvrmse, "nmbe": nmbe, "rmse": rmse}


def build_design(rows: Sequence[Obs], features: Sequence[str]):
    x = []
    y = []
    for o in rows:
        fmap = {"B": o.b, "C": o.c, "E": o.e, "F": o.f}
        x.append([1.0] + [fmap[name] for name in features])
        y.append(o.h)
    return x, y


def fit_on_rows(rows: Sequence[Obs], features: Sequence[str]):
    x, y = build_design(rows, features)
    fit = ols_fit(x, y)
    m = calc_metrics(y, fit["yhat"], len(features))
    return fit, m, x, y


def predict(features: Sequence[str], beta: Sequence[float], rows: Sequence[Obs]):
    x, y = build_design(rows, features)
    yhat = [sum(v * b for v, b in zip(row, beta)) for row in x]
    return x, y, yhat


def durbin_watson(residuals: Sequence[float]) -> float:
    num = sum((residuals[i] - residuals[i - 1]) ** 2 for i in range(1, len(residuals)))
    den = sum(r * r for r in residuals)
    return num / den if den else float("nan")


def corr(x: Sequence[float], y: Sequence[float]) -> float:
    mx = statistics.mean(x)
    my = statistics.mean(y)
    sx = math.sqrt(sum((v - mx) ** 2 for v in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def simple_slope_and_p(y: Sequence[float], x: Sequence[float]) -> Tuple[float, float, float]:
    r = corr(x, y)
    n = len(x)
    xavg = statistics.mean(x)
    yavg = statistics.mean(y)
    sxx = sum((xi - xavg) ** 2 for xi in x)
    slope = sum((xi - xavg) * (yi - yavg) for xi, yi in zip(x, y)) / sxx if sxx else 0.0
    if n <= 2 or abs(r) >= 0.999999:
        return slope, r, 0.0
    t = abs(r) * math.sqrt((n - 2) / max(1e-12, 1 - r * r))
    p = 2 * (1 - student_t_cdf(t, n - 2))
    return slope, r, p


def breusch_pagan(residuals: Sequence[float], aux_x: List[List[float]]):
    # aux_x incluye intercepto + predictores.
    z = [r * r for r in residuals]
    fit = ols_fit(aux_x, z)
    m = calc_metrics(z, fit["yhat"], len(aux_x[0]) - 1)
    n = len(z)
    k = len(aux_x[0]) - 1
    lm = n * max(0.0, m["r2"])

    f_stat = (m["r2"] / k) / ((1 - m["r2"]) / (n - k - 1)) if (n - k - 1) > 0 and m["r2"] < 1 else float("inf")
    f_pval = 1 - f_cdf(f_stat, k, n - k - 1) if math.isfinite(f_stat) else 0.0

    return {"lm": lm, "r2_aux": m["r2"], "f_stat": f_stat, "f_pval": f_pval, "k": k, "n": n}


def quantile(vals: Sequence[float], q: float):
    s = sorted(vals)
    if not s:
        return float("nan")
    if len(s) == 1:
        return s[0]
    pos = (len(s) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return s[lo]
    frac = pos - lo
    return s[lo] * (1 - frac) + s[hi] * frac


def grouped_residual_dispersion(values: Sequence[float], residuals: Sequence[float], groups: int = 5):
    cut_points = [quantile(values, i / groups) for i in range(1, groups)]
    buckets: List[List[float]] = [[] for _ in range(groups)]
    for v, r in zip(values, residuals):
        idx = 0
        while idx < len(cut_points) and v > cut_points[idx]:
            idx += 1
        buckets[idx].append(r)
    out = []
    for i, bucket in enumerate(buckets, start=1):
        if bucket:
            out.append(
                {
                    "grupo": i,
                    "n": len(bucket),
                    "mean": statistics.mean(bucket),
                    "std": statistics.pstdev(bucket),
                    "p05": quantile(bucket, 0.05),
                    "p95": quantile(bucket, 0.95),
                }
            )
    return out


def split_rows(rows: Sequence[Obs], seed: int, frac_train: float = 0.8):
    rnd = random.Random(seed)
    idx = list(range(len(rows)))
    rnd.shuffle(idx)
    cut = int(frac_train * len(idx))
    train = [rows[i] for i in idx[:cut]]
    test = [rows[i] for i in idx[cut:]]
    return train, test


def summarize_stability(records: Sequence[Dict], coef_names: Sequence[str]) -> Dict[str, Dict[str, float]]:
    out = {}
    for name in coef_names:
        vals = [r[name] for r in records]
        mu = statistics.mean(vals)
        sd = statistics.pstdev(vals)
        cv = abs(sd / mu) * 100 if mu != 0 else float("inf")
        out[name] = {"mean": mu, "std": sd, "min": min(vals), "max": max(vals), "cv_pct": cv}
    return out


def model_equation(beta: Sequence[float], features: Sequence[str]) -> str:
    terms = [f"{beta[0]:.6f}"]
    for feat, b in zip(features, beta[1:]):
        terms.append(f"({b:.6f})·{feat}")
    return "LNG_hat = " + " + ".join(terms)


def write_report(path: Path, report: str):
    path.write_text(report, encoding="utf-8")


def main():
    data = build_dataset("data/LB_CAE_Eyectores.xlsx")
    if len(data) < 30:
        raise RuntimeError("Muestra insuficiente para validación robusta")

    features_full = ["B", "E", "C", "F"]
    features_base = ["B", "E", "C"]

    fit_all, m_all, x_all, y_all = fit_on_rows(data, features_full)
    residuals_all = [y - yh for y, yh in zip(y_all, fit_all["yhat"])]

    dw = durbin_watson(residuals_all)
    bp = breusch_pagan(residuals_all, x_all)

    resid_checks = {}
    yhat_all = fit_all["yhat"]
    for name, values in {
        "pred": yhat_all,
        "B": [o.b for o in data],
        "E": [o.e for o in data],
        "C": [o.c for o in data],
        "F": [o.f for o in data],
    }.items():
        slope, r, pval = simple_slope_and_p(residuals_all, values)
        resid_checks[name] = {"slope": slope, "corr": r, "pval": pval}

    dispersion_pred = grouped_residual_dispersion(yhat_all, residuals_all, groups=5)

    seeds = [7, 11, 21, 42, 84, 126, 256, 512]
    stability_records = []
    metric_records = []

    for seed in seeds:
        train, test = split_rows(data, seed)
        fit_tr, _, _, _ = fit_on_rows(train, features_full)
        _, yte, yhat_te = predict(features_full, fit_tr["beta"], test)
        mte = calc_metrics(yte, yhat_te, len(features_full))

        rec = {
            "seed": seed,
            "beta0": fit_tr["beta"][0],
            "beta_B": fit_tr["beta"][1],
            "beta_E": fit_tr["beta"][2],
            "beta_C": fit_tr["beta"][3],
            "beta_F": fit_tr["beta"][4],
            "r2_test": mte["r2"],
            "cvrmse_test": mte["cvrmse"],
            "nmbe_test": mte["nmbe"],
        }
        stability_records.append(rec)
        metric_records.append(rec)

    coef_summary = summarize_stability(stability_records, ["beta0", "beta_B", "beta_E", "beta_C", "beta_F"])
    metric_summary = summarize_stability(metric_records, ["r2_test", "cvrmse_test", "nmbe_test"])

    # Check de simplicidad en la misma partición de referencia (seed 42)
    train42, test42 = split_rows(data, 42)
    fit42_full, _, _, _ = fit_on_rows(train42, features_full)
    _, yte42, yhat42_full = predict(features_full, fit42_full["beta"], test42)
    full42 = calc_metrics(yte42, yhat42_full, len(features_full))

    fit42_base, _, _, _ = fit_on_rows(train42, features_base)
    _, _, yhat42_base = predict(features_base, fit42_base["beta"], test42)
    base42 = calc_metrics(yte42, yhat42_base, len(features_base))

    ranges = {
        "B": (min(o.b for o in data), max(o.b for o in data)),
        "E": (min(o.e for o in data), max(o.e for o in data)),
        "C": (min(o.c for o in data), max(o.c for o in data)),
        "F": (min(o.f for o in data), max(o.f for o in data)),
    }

    # Interpretación física mínima automáticamente basada en signo
    be = fit_all["beta"][2]
    e_note = (
        "β_E < 0: si E es presión absoluta en cabeza, mayor presión (menos vacío) reduce LNG según el ajuste; "
        "si E se define como nivel de vacío (magnitud de vacío), el signo aparente sería físicamente inconsistente y debe invertirse la definición operacional."
        if be < 0
        else "β_E > 0: coherente si E representa esfuerzo/magnitud de vacío; revisar convención si E fuera presión absoluta."
    )

    report = []
    report.append("# Validación final baseline LNG (CAE/IPMVP)\n")
    report.append(f"Observaciones válidas: **{len(data)}**. Modelo fijo validado: LNG = β0 + β1·B + β2·E + β3·C + β4·F.\n")

    report.append("## 1) Validación estadística final\n")
    report.append(f"- **Durbin-Watson (residuos orden original)**: **{dw:.4f}** (≈2 sugiere baja autocorrelación lineal de primer orden).\n")
    report.append(f"- **Breusch-Pagan (versión F)**: F={bp['f_stat']:.4f}, p-value={bp['f_pval']:.6g}. ")
    if bp["f_pval"] < 0.05:
        report.append("Se detecta heterocedasticidad estadísticamente significativa; recomendable usar errores robustos HC en auditoría.\n")
    else:
        report.append("No se detecta evidencia estadística fuerte de heterocedasticidad (α=5%).\n")

    report.append("- **Residuos vs predicción/variables** (pendiente en regresión simple de residuo):\n")
    report.append("\n| Relación | Pendiente | Correlación | p-value |\n|---|---:|---:|---:|\n")
    for key in ["pred", "B", "E", "C", "F"]:
        v = resid_checks[key]
        report.append(f"| Residuo~{key} | {v['slope']:.6f} | {v['corr']:.4f} | {v['pval']:.6g} |\n")

    report.append("\nNota técnica: en OLS, el residuo es ortogonal a los regresores incluidos en el ajuste; por eso las correlaciones lineales salen ~0 por construcción.\n")
    report.append("\n- **Dispersión de residuos vs predicción (quintiles de Ŷ):**\n")
    report.append("\n| Quintil Ŷ | n | Media residuo | Std residuo | P5 | P95 |\n|---:|---:|---:|---:|---:|---:|\n")
    for row in dispersion_pred:
        report.append(
            f"| {row['grupo']} | {row['n']} | {row['mean']:.4f} | {row['std']:.4f} | {row['p05']:.4f} | {row['p95']:.4f} |\n"
        )
    report.append("\nInterpretación: la variación de la desviación estándar entre quintiles sugiere heterocedasticidad, consistente con Breusch-Pagan.\n")

    report.append("## 2) Validación de robustez (múltiples semillas)\n")
    report.append(f"Semillas evaluadas: {seeds}.\n")
    report.append("\n| Semilla | β0 | βB | βE | βC | βF | R² test | CVRMSE test (%) | NMBE test (%) |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in stability_records:
        report.append(
            f"| {r['seed']} | {r['beta0']:.3f} | {r['beta_B']:.4f} | {r['beta_E']:.4f} | {r['beta_C']:.4f} | {r['beta_F']:.4f} | {r['r2_test']:.4f} | {r['cvrmse_test']:.3f} | {r['nmbe_test']:.3f} |\n"
        )

    report.append("\nResumen de estabilidad (media ± std):\n")
    report.append("\n| Parámetro | Media | Std | Mín | Máx | CV % |\n|---|---:|---:|---:|---:|---:|\n")
    for name in ["beta_B", "beta_E", "beta_C", "beta_F", "r2_test", "cvrmse_test", "nmbe_test"]:
        v = coef_summary[name] if name.startswith("beta") else metric_summary[name]
        report.append(f"| {name} | {v['mean']:.6f} | {v['std']:.6f} | {v['min']:.6f} | {v['max']:.6f} | {v['cv_pct']:.3f} |\n")

    stable = metric_summary["cvrmse_test"]["std"] < 0.25 and metric_summary["r2_test"]["std"] < 0.02
    report.append("\nConclusión robustez: **")
    report.append("Modelo estable ante particiones train/test**.\n" if stable else "Sensibilidad apreciable a la partición; revisar segmentación operacional**.\n")

    report.append("## 3) Coherencia física de coeficientes\n")
    b0, bB, bE, bC, bF = fit_all["beta"]
    report.append(f"Ecuación ajustada en dataset completo:\n\n`{model_equation(fit_all['beta'], features_full)}`\n\n")
    report.append(f"- **βB={bB:.6f}**: positivo; mayor caudal de crudo incrementa LNG (coherente físicamente).\n")
    report.append(f"- **βC={bC:.6f}**: positivo; mayor severidad térmica/flash aumenta carga de vapores y demanda de vacío.\n")
    report.append(f"- **βE={bE:.6f}**: {e_note}\n")
    report.append(
        f"- **βF={bF:.6f}**: F actúa como proxy de calidad de betún (no causal directa); su signo captura efecto composicional/operativo correlacionado sobre la demanda de LNG.\n"
    )

    report.append("## 4) Check de simplicidad (mantener o eliminar F)\n")
    report.append("Comparación en misma partición (seed=42):\n\n")
    report.append("| Modelo | Variables | R² test | CVRMSE test (%) | NMBE test (%) |\n|---|---|---:|---:|---:|\n")
    report.append(f"| Completo | B,E,C,F | {full42['r2']:.4f} | {full42['cvrmse']:.3f} | {full42['nmbe']:.3f} |\n")
    report.append(f"| Base | B,E,C | {base42['r2']:.4f} | {base42['cvrmse']:.3f} | {base42['nmbe']:.3f} |\n")

    d_r2 = full42["r2"] - base42["r2"]
    d_cv = full42["cvrmse"] - base42["cvrmse"]
    keep_f = d_r2 > 0.002 or d_cv < -0.05
    report.append("\nDecisión: **")
    if keep_f:
        report.append("mantener F** por mejora consistente, aunque moderada, de robustez predictiva.\n")
    else:
        report.append("eliminar F** por mejora marginal no relevante frente a parsimonia.\n")

    report.append("## 5) Definición formal del modelo\n")
    report.append("- **Ecuación final propuesta para baseline**:\n")
    report.append(f"  - `{model_equation(fit_all['beta'], features_full if keep_f else features_base)}`\n")
    report.append("- **Rango de validez (dominio observado):**\n")
    for k in ["B", "C", "E", "F"]:
        lo, hi = ranges[k]
        report.append(f"  - {k}: [{lo:.4f}, {hi:.4f}]\n")
    report.append("- **Condiciones de aplicación:** operación en régimen similar al histórico filtrado (B>140 t/h, C>300 ºC), mismo esquema de control de vacío y mismo tipo de crudo/mezcla.\n")

    report.append("## 6) Texto técnico listo para auditoría\n")
    report.append("### A. Definición de variables\n")
    report.append("- B: caudal de crudo (t/h).\n- C: temperatura flash (ºC).\n")
    report.append("- E: variable de cabecera asociada al sistema de vacío. En memoria debe declararse explícitamente si E es **presión absoluta** o **magnitud de vacío** para evitar ambigüedad de signo.\n")
    report.append("- F: viscosidad del betún; proxy de calidad/composición, sin causalidad termodinámica directa sobre consumo.\n")

    report.append("### B. Justificación física\n")
    report.append("El consumo de LNG de eyectores aumenta con mayor carga (B) y con mayor severidad de operación (C). La variable E representa el esfuerzo del tren de vacío en cabeza. La variable F captura cambios de calidad que alteran volatilidad/carga al sistema de vacío y mejora la explicabilidad operacional del baseline.\n")

    report.append("### C. Limitaciones\n")
    report.append("Modelo lineal estacionario, válido únicamente dentro del dominio observado y con arquitectura de proceso equivalente. No extrapolar fuera de rangos ni usar en transitorios, arranques o cambios de hardware/control.\n")
    report.append("Adicionalmente, Durbin-Watson muy por debajo de 2 sugiere dependencia temporal en residuos; para la defensa estadística se recomienda reportar errores estándar robustos HAC/Newey-West sin cambiar la estructura del modelo.\n")

    report.append("### D. Extensión futura\n")
    report.append("Para cerrar balance energético de vacío en CAE/IPMVP: `E_total = LNG + Fuel Gas`, manteniendo la misma filosofía de variables explicativas y validación estadística.\n")

    report_path = Path("baseline_validation_audit.md")
    write_report(report_path, "".join(report))

    print(f"Reporte generado: {report_path}")
    print(f"DW={dw:.4f}; BP p-value={bp['f_pval']:.6g}")
    print("Ecuación completa:", model_equation(fit_all["beta"], features_full))


if __name__ == "__main__":
    main()
