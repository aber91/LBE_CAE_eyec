# Línea Base CAE - Eyectores de vacío

## 1) Preparación de datos

- Archivo: `data/LB_CAE_Eyectores.xlsx`
- Filtros aplicados: **Crude feed (B) > 140 t/h** y **Tflash (C) > 300 ºC**.
- Exclusión adicional para estabilidad numérica: F, G y H > 0.
- Observaciones válidas finales: **9937**.

## 2) Calidad de betún (viscosidad -> PEN -> grado)

- La hoja `Visc_data` se usa para interpolar PEN desde viscosidad (F).
- Grados considerados: 15/25, 35/50, 50/70, 70/100 y 160/220.

| Grado | Nº muestras |
|---|---:|
| 15/25 | 1498 |
| 35/50 | 3859 |
| 50/70 | 2045 |
| 70/100 | 1062 |
| 160/220 | 1473 |

## 3) Modelo LB (OLS con no linealidades e interacciones)

**Variable objetivo:** H (LNG flow rate, kg/h).

Variables usadas: 1, BC, EF, DF, grade_160_220, G_logF, grade_35_50, pen2, B2, B, C2.

| Término | Coeficiente |
|---|---:|
| 1 | 1287.80722535 |
| BC | 0.03875603 |
| EF | 0.00753193 |
| DF | -0.00136287 |
| grade_160_220 | -14.86668323 |
| G_logF | -0.14568771 |
| grade_35_50 | 10.39789756 |
| pen2 | 0.00005543 |
| B2 | 0.02012243 |
| B | -13.40596180 |
| C2 | -0.00524184 |

## 4) Métricas estadísticas

| Escenario | R² | R² ajustado | CVRMSE (%) | NMBE (%) |
|---|---:|---:|---:|---:|
| Ajuste completo LB | 0.755338 | 0.755067 | 4.604 | -0.000 |
| Train (80%) | 0.712741 | 0.712343 | 4.747 | -0.000 |
| Test (20%) | 0.844173 | 0.843306 | 3.933 | 0.192 |

## 5) Criterio de aceptación CAE (estadístico)

- Referencia habitual M&V: **R² > 0.75**, **CVRMSE < 20%**, **|NMBE| < 5%**.
- El modelo **cumple holgadamente CVRMSE y NMBE**; y cumple R² en calibración completa LB.

## 6) Ecuación de predicción

LNG_hat = 1287.807225 + (0.038756)·BC + (0.007532)·EF + (-0.001363)·DF + (-14.866683)·grade_160_220 + (-0.145688)·G_logF + (10.397898)·grade_35_50 + (0.000055)·pen2 + (0.020122)·B2 + (-13.405962)·B + (-0.005242)·C2
