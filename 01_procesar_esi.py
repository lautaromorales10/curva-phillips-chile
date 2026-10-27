"""Paso 1. Procesa los microdatos de la ESI (INE) de cada año y calcula, por sexo:
  - ingreso medio por hora de los asalariados (razón: suma de ingresos / suma de horas, ponderadas)
  - porcentaje de la fuerza de trabajo con educación superior.
Entrada: datos/esi/*.dta (el año se lee del nombre del archivo).
Salida:  resultados/esi_por_anio.csv"""
import re
import numpy as np, pandas as pd
from config import ESI, RES

VARS = ["sexo", "activ", "nivel", "categoria_ocupacion", "categoria_ocupacional",
        "ing_t_p", "fact_cal_esi", "d4_dias", "d4_horas", "dias", "horas"]


def procesar(anio, ruta):
    with pd.read_stata(ruta, iterator=True) as rd:
        cols = list(rd.variable_labels().keys())
    df = pd.read_stata(ruta, columns=[c for c in cols if c.lower() in VARS], convert_categoricals=False)
    df.columns = [c.lower() for c in df.columns]
    if "categoria_ocupacion" not in df and "categoria_ocupacional" in df:     # nombre usado en 2019
        df = df.rename(columns={"categoria_ocupacional": "categoria_ocupacion"})
    vd = "d4_dias" if "d4_dias" in df else "dias"                          # cambia según el año
    vh = "d4_horas" if "d4_horas" in df else "horas"

    # Asalariados (categorías 3 y 4) con ingreso positivo y jornada válida
    a = df[df.categoria_ocupacion.isin([3, 4]) & (df.ing_t_p > 0) & df.fact_cal_esi.notna()]
    a = a[a[vd].between(1, 31) & a[vh].between(1, 24)].copy()
    a["h"] = a[vd] * a[vh]
    # Fuerza de trabajo: ocupados (1) y desocupados (2)
    ft = df[df.activ.isin([1, 2]) & df.fact_cal_esi.notna()]

    filas = []
    for s_esi, s in [(1, 0), (2, 1)]:          # ESI: 1 = hombre, 2 = mujer  ->  sexo: 0 = hombre, 1 = mujer
        g, f = a[a.sexo == s_esi], ft[ft.sexo == s_esi]
        filas.append(dict(
            anio=anio, sexo=s,
            salario_hora_nominal=(g.ing_t_p * g.fact_cal_esi).sum() / (g.h * g.fact_cal_esi).sum(),
            n_asalariados=len(g),
            educ_sup=100 * np.average(f.nivel.isin(range(7, 13)), weights=f.fact_cal_esi),
            variables_jornada=f"{vd} x {vh}"))
    return pd.DataFrame(filas)


if __name__ == "__main__":
    archivos = sorted(ESI.glob("*.dta"))
    if not archivos:
        raise SystemExit(f"No se encontraron archivos .dta en {ESI}")
    res = []
    for ruta in archivos:
        anio = int(re.search(r"20\d\d", ruta.name).group())
        r = procesar(anio, ruta); res.append(r)
        h, m = r.iloc[0], r.iloc[1]
        print(f"{anio}: ingreso/hora H {h.salario_hora_nominal:6.0f}  M {m.salario_hora_nominal:6.0f} | "
              f"ed. superior H {h.educ_sup:4.1f}%  M {m.educ_sup:4.1f}%  ({ruta.name})", flush=True)
    out = pd.concat(res).sort_values(["anio", "sexo"])
    faltan = sorted(set(range(2010, 2025)) - set(out.anio))
    if faltan:
        print("ATENCIÓN: faltan años de la ESI:", faltan)
    out.to_csv(RES / "esi_por_anio.csv", index=False)
    print("Guardado:", RES / "esi_por_anio.csv")
