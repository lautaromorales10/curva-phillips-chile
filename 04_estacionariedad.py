"""Paso 4. Pruebas ADF y KPSS por serie y sexo, en niveles y en primeras diferencias.
Salida: resultados/estacionariedad.csv"""
import numpy as np, pandas as pd
from config import RES
from econometria import adf, kpss

b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
b["salario_real_hora"] = np.log(b.salario_real_hora)
# Componentes determinísticos en niveles: constante (c) o constante y tendencia (ct) para series con tendencia
SERIES = {"desempleo": "c", "inflacion": "c", "educ_sup": "ct", "salario_real_hora": "ct"}

filas = []
for s, etq in [(0, "H"), (1, "M")]:
    d = b[b.sexo == s]
    for serie, det in SERIES.items():
        if serie == "inflacion" and s == 1:           # la inflación es común a ambos sexos
            continue
        x = d[serie].values
        for forma, z, reg in [("nivel", x, det), ("diferencia", np.diff(x), "c")]:
            t, cr = adf(z, reg); k, kc = kpss(z, reg)
            filas.append(dict(serie=serie, sexo="común" if serie == "inflacion" else etq, forma=forma, det=reg,
                              adf=round(t, 2), crit5=round(cr[5], 2), adf_rechaza=t < cr[5],
                              kpss=round(k, 3), kpss_crit5=kc[5], kpss_rechaza=k > kc[5]))
e = pd.DataFrame(filas)
e.to_csv(RES / "estacionariedad.csv", index=False)
print(e.to_string(index=False))
