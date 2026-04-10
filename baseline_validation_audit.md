# Validación final baseline LNG (CAE/IPMVP)
Observaciones válidas: **9937**. Modelo fijo validado: LNG = β0 + β1·B + β2·E + β3·C + β4·F.
## 1) Validación estadística final
- **Durbin-Watson (residuos orden original)**: **0.2802** (≈2 sugiere baja autocorrelación lineal de primer orden).
- **Breusch-Pagan (versión F)**: F=51.2730, p-value=0. Se detecta heterocedasticidad estadísticamente significativa; recomendable usar errores robustos HC en auditoría.
- **Residuos vs predicción/variables** (pendiente en regresión simple de residuo):

| Relación | Pendiente | Correlación | p-value |
|---|---:|---:|---:|
| Residuo~pred | -0.000000 | -0.0000 | 1 |
| Residuo~B | -0.000000 | -0.0000 | 1 |
| Residuo~E | 0.000000 | 0.0000 | 1 |
| Residuo~C | -0.000000 | -0.0000 | 1 |
| Residuo~F | -0.000000 | -0.0000 | 1 |

Nota técnica: en OLS, el residuo es ortogonal a los regresores incluidos en el ajuste; por eso las correlaciones lineales salen ~0 por construcción.

- **Dispersión de residuos vs predicción (quintiles de Ŷ):**

| Quintil Ŷ | n | Media residuo | Std residuo | P5 | P95 |
|---:|---:|---:|---:|---:|---:|
| 1 | 1988 | 7.6589 | 70.2118 | -88.6361 | 121.3619 |
| 2 | 1987 | -14.6042 | 54.4367 | -104.5162 | 68.5793 |
| 3 | 1987 | -3.2721 | 55.8572 | -97.5906 | 81.9795 |
| 4 | 1987 | -5.4647 | 52.7933 | -90.5561 | 82.8505 |
| 5 | 1988 | 15.6703 | 47.4724 | -54.7624 | 93.3489 |

Interpretación: la variación de la desviación estándar entre quintiles sugiere heterocedasticidad, consistente con Breusch-Pagan.
## 2) Validación de robustez (múltiples semillas)
Semillas evaluadas: [7, 11, 21, 42, 84, 126, 256, 512].

| Semilla | β0 | βB | βE | βC | βF | R² test | CVRMSE test (%) | NMBE test (%) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | -1104.774 | 5.7129 | 0.5611 | 3.7311 | 0.0245 | 0.7397 | 4.802 | 0.102 |
| 11 | -1128.617 | 5.7266 | 0.6759 | 3.7886 | 0.0236 | 0.7191 | 4.929 | 0.185 |
| 21 | -1124.631 | 5.7363 | 0.6196 | 3.7765 | 0.0231 | 0.7279 | 4.849 | 0.120 |
| 42 | -1118.184 | 5.7262 | 0.6083 | 3.7623 | 0.0238 | 0.7336 | 4.903 | 0.078 |
| 84 | -1127.794 | 5.7881 | 0.3441 | 3.7747 | 0.0235 | 0.7212 | 5.003 | 0.195 |
| 126 | -1123.706 | 5.7399 | 0.5203 | 3.7830 | 0.0204 | 0.7414 | 4.760 | 0.119 |
| 256 | -1122.360 | 5.7426 | 0.6593 | 3.7645 | 0.0224 | 0.7446 | 4.696 | 0.162 |
| 512 | -1114.830 | 5.7668 | 0.3981 | 3.7437 | 0.0236 | 0.7335 | 4.840 | 0.180 |

Resumen de estabilidad (media ± std):

| Parámetro | Media | Std | Mín | Máx | CV % |
|---|---:|---:|---:|---:|---:|
| beta_B | 5.742429 | 0.022677 | 5.712874 | 5.788124 | 0.395 |
| beta_E | 0.548344 | 0.113176 | 0.344143 | 0.675931 | 20.640 |
| beta_C | 3.765556 | 0.018417 | 3.731137 | 3.788621 | 0.489 |
| beta_F | 0.023140 | 0.001157 | 0.020449 | 0.024498 | 5.000 |
| r2_test | 0.732621 | 0.008717 | 0.719081 | 0.744595 | 1.190 |
| cvrmse_test | 4.847602 | 0.091262 | 4.695739 | 5.002513 | 1.883 |
| nmbe_test | 0.142384 | 0.040590 | 0.077920 | 0.194709 | 28.507 |

Conclusión robustez: **Modelo estable ante particiones train/test**.
## 3) Coherencia física de coeficientes
Ecuación ajustada en dataset completo:

`LNG_hat = -1119.741341 + (5.733735)·B + (0.607600)·E + (3.766043)·C + (0.022418)·F`

- **βB=5.733735**: positivo; mayor caudal de crudo incrementa LNG (coherente físicamente).
- **βC=3.766043**: positivo; mayor severidad térmica/flash aumenta carga de vapores y demanda de vacío.
- **βE=0.607600**: β_E > 0: coherente si E representa esfuerzo/magnitud de vacío; revisar convención si E fuera presión absoluta.
- **βF=0.022418**: F actúa como proxy de calidad de betún (no causal directa); su signo captura efecto composicional/operativo correlacionado sobre la demanda de LNG.
## 4) Check de simplicidad (mantener o eliminar F)
Comparación en misma partición (seed=42):

| Modelo | Variables | R² test | CVRMSE test (%) | NMBE test (%) |
|---|---|---:|---:|---:|
| Completo | B,E,C,F | 0.7336 | 4.903 | 0.078 |
| Base | B,E,C | 0.7331 | 4.908 | 0.060 |

Decisión: **eliminar F** por mejora marginal no relevante frente a parsimonia.
## 5) Definición formal del modelo
- **Ecuación final propuesta para baseline**:
  - `LNG_hat = -1119.741341 + (5.733735)·B + (0.607600)·E + (3.766043)·C`
- **Rango de validez (dominio observado):**
  - B: [140.8821, 218.9184]
  - C: [300.1408, 355.1150]
  - E: [15.2801, 101.3323]
  - F: [52.5193, 1591.8065]
- **Condiciones de aplicación:** operación en régimen similar al histórico filtrado (B>140 t/h, C>300 ºC), mismo esquema de control de vacío y mismo tipo de crudo/mezcla.
## 6) Texto técnico listo para auditoría
### A. Definición de variables
- B: caudal de crudo (t/h).
- C: temperatura flash (ºC).
- E: variable de cabecera asociada al sistema de vacío. En memoria debe declararse explícitamente si E es **presión absoluta** o **magnitud de vacío** para evitar ambigüedad de signo.
- F: viscosidad del betún; proxy de calidad/composición, sin causalidad termodinámica directa sobre consumo.
### B. Justificación física
El consumo de LNG de eyectores aumenta con mayor carga (B) y con mayor severidad de operación (C). La variable E representa el esfuerzo del tren de vacío en cabeza. La variable F captura cambios de calidad que alteran volatilidad/carga al sistema de vacío y mejora la explicabilidad operacional del baseline.
### C. Limitaciones
Modelo lineal estacionario, válido únicamente dentro del dominio observado y con arquitectura de proceso equivalente. No extrapolar fuera de rangos ni usar en transitorios, arranques o cambios de hardware/control.
Adicionalmente, Durbin-Watson muy por debajo de 2 sugiere dependencia temporal en residuos; para la defensa estadística se recomienda reportar errores estándar robustos HAC/Newey-West sin cambiar la estructura del modelo.
### D. Extensión futura
Para cerrar balance energético de vacío en CAE/IPMVP: `E_total = LNG + Fuel Gas`, manteniendo la misma filosofía de variables explicativas y validación estadística.
