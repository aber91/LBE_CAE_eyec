# Baseline CAE/IPMVP por calidad de betún (F)

## Metodología de segregación
- Dataset válido global: **9937** observaciones.
- Segmentación por terciles de viscosidad de betún (F):
  - **Q1 (baja): F ≤ 471.654**
  - **Q2 (media): 471.654 < F ≤ 627.431**
  - **Q3 (alta): F > 627.431**

- En cada segmento se ajusta y valida una línea base independiente (split 80/20, semilla 42).

## Q1_baja

- Observaciones válidas en segmento: **3313**

| Modelo | Variables | R² test | CVRMSE test (%) | NMBE test (%) |
|---|---|---:|---:|---:|
| A_base | B, E, C | 0.6698 | 5.64 | 0.01 |
| B_interaccion_BE | B, E, C, BxE | 0.6718 | 5.62 | -0.01 |
| C_con_caudal_G | B, E, C, G | 0.6766 | 5.58 | -0.00 |
| D_reducido | B, E | 0.5651 | 6.47 | -0.00 |

**Modelo base recomendado (Q1_baja):** C_con_caudal_G

```
LNG_hat = -778.449925 + (6.177628)·B + (1.609011)·E + (2.978283)·C + (-1.640918)·G
```

| Escenario | R² | adj R² | CVRMSE (%) | NMBE (%) |
|---|---:|---:|---:|---:|
| Train | 0.7259 | 0.7255 | 5.16 | -0.00 |
| Test | 0.6766 | 0.6746 | 5.58 | -0.00 |
| Global | 0.7162 | 0.7159 | 5.24 | -0.00 |

## Q2_media

- Observaciones válidas en segmento: **3312**

| Modelo | Variables | R² test | CVRMSE test (%) | NMBE test (%) |
|---|---|---:|---:|---:|
| A_base | B, E, C | 0.7352 | 4.83 | -0.34 |
| B_interaccion_BE | B, E, C, BxE | 0.7357 | 4.82 | -0.33 |
| C_con_caudal_G | B, E, C, G | 0.7399 | 4.78 | -0.34 |
| D_reducido | B, E | 0.6915 | 5.21 | -0.16 |

**Modelo base recomendado (Q2_media):** C_con_caudal_G

```
LNG_hat = -861.574395 + (7.111213)·B + (-0.016012)·E + (2.721425)·C + (-1.227930)·G
```

| Escenario | R² | adj R² | CVRMSE (%) | NMBE (%) |
|---|---:|---:|---:|---:|
| Train | 0.7394 | 0.7390 | 4.75 | 0.00 |
| Test | 0.7399 | 0.7383 | 4.78 | -0.34 |
| Global | 0.7396 | 0.7393 | 4.75 | -0.07 |

## Q3_alta

- Observaciones válidas en segmento: **3312**

| Modelo | Variables | R² test | CVRMSE test (%) | NMBE test (%) |
|---|---|---:|---:|---:|
| A_base | B, E, C | 0.7774 | 3.79 | 0.13 |
| B_interaccion_BE | B, E, C, BxE | 0.7743 | 3.82 | 0.13 |
| C_con_caudal_G | B, E, C, G | 0.7770 | 3.80 | 0.12 |
| D_reducido | B, E | 0.6869 | 4.50 | 0.15 |

**Modelo base recomendado (Q3_alta):** A_base

```
LNG_hat = -1000.606668 + (5.965616)·B + (-0.925459)·E + (3.455321)·C
```

| Escenario | R² | adj R² | CVRMSE (%) | NMBE (%) |
|---|---:|---:|---:|---:|
| Train | 0.7366 | 0.7363 | 4.15 | 0.00 |
| Test | 0.7774 | 0.7764 | 3.79 | 0.13 |
| Global | 0.7449 | 0.7447 | 4.08 | 0.03 |

