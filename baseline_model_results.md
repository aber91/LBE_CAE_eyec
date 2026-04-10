# Baseline CAE/IPMVP - Eyectores de vacío

## Contexto y alcance
- Objetivo actual: baseline de **LNG (H)** para el sistema de vacío.
- Extensión futura prevista: **E_total = LNG + Fuel Gas** con la misma estructura de drivers físicos.

## Datos y validación
- Observaciones válidas tras filtros operativos: **9937**.
- Split: **train/test = 80/20** (aleatorio con semilla fija 42).
- Filtros: B > 140 t/h, C > 300 ºC, F>0, G>0, H>0.

## Comparativa de modelos candidatos

| Modelo | Variables finales | R² test | adjR² test | CVRMSE test (%) | NMBE test (%) | #vars |
|---|---|---:|---:|---:|---:|---:|
| C_calidad_viscosidad | B, E, C, F | 0.7336 | 0.7331 | 4.90 | 0.08 | 4 |
| A_base | B, E, C | 0.7331 | 0.7327 | 4.91 | 0.06 | 3 |
| B_interaccion | B, E, C | 0.7331 | 0.7327 | 4.91 | 0.06 | 3 |
| C_calidad_dummy | B, E, C, visc_hi_dummy | 0.7313 | 0.7308 | 4.92 | 0.06 | 4 |
| D_reducido | B, E | 0.5458 | 0.5453 | 6.40 | 0.07 | 2 |

## Modelo final recomendado

**Seleccionado:** C_calidad_viscosidad

### Ecuación

```
LNG_hat = -1118.184492 + (5.726244)·B + (0.608335)·E + (3.762346)·C + (0.023843)·F
```

### Coeficientes y significancia (train)

| Término | Coeficiente | p-value (aprox) |
|---|---:|---:|
| Intercepto | -1118.184492 | 0 |
| B | 5.726244 | 0 |
| E | 0.608335 | 2.49109e-05 |
| C | 3.762346 | 0 |
| F | 0.023843 | 9.53237e-12 |

### Métricas (modelo final)

| Escenario | R² | adj R² | CVRMSE (%) | NMBE (%) |
|---|---:|---:|---:|---:|
| Train | 0.7345 | 0.7343 | 4.77 | 0.00 |
| Test | 0.7336 | 0.7331 | 4.90 | 0.08 |
| Global | 0.7343 | 0.7342 | 4.80 | 0.02 |

### Diagnóstico de multicolinealidad (VIF)

| Variable | VIF |
|---|---:|
| B | 1.148 |
| C | 1.736 |
| E | 1.049 |
| F | 1.861 |

### Residuos en test

- Residuo medio: **0.9334**
- MAE residuos: **43.1848**
- Rango residuos: **[-214.7337, 486.0523]**

## Justificación técnica defendible en auditoría

- **Carga (B):** principal driver de demanda energética del sistema de vacío.
- **Presión de cabeza (E):** representa esfuerzo del tren de eyectores para sostener vacío operativo.
- **Temperatura flash (C):** refleja severidad térmica y carga de vapores al sistema.
- **Interacción BxE** (si aplica): solo se mantiene cuando mejora estabilidad y tiene sentido físico (más carga a menor presión efectiva implica mayor exigencia).
- **Calidad (F o dummy):** se incluye solo si aporta señal estadística robusta sin degradar VIF.
- Se prioriza **parsimonia (4–5 variables máx.)** para trazabilidad y robustez CAE/IPMVP.
