"""Paso 2. Construye la base anual apilada por sexo, 2010-2024 (30 observaciones).
Entradas: datos/desocupacion_hombres.xlsx, datos/desocupacion_mujeres.xlsx, datos/ipc_general.xlsx,
          datos/poblacion_hombres_y_mujeres.xlsx (descargas de la BDE del Banco Central)
          resultados/esi_por_anio.csv (paso 1)
Salida:   resultados/base_apilada_2010_2024.csv"""
import numpy as np, pandas as pd
from config import DATOS, RES


def leer(archivo, nombres):
    d = pd.read_excel(DATOS / archivo, sheet_name="Cuadro", header=None, skiprows=3)
    d = d.iloc[:, : len(nombres) + 1]
    d.columns = ["fecha"] + nombres
    d["fecha"] = pd.to_datetime(d["fecha"], errors="coerce")
    d = d.dropna(subset=["fecha"])
    d["anio"], d["mes"] = d.fecha.dt.year, d.fecha.dt.month
    return d


um = leer("desocupacion_mujeres.xlsx", ["u"])
uh = leer("desocupacion_hombres.xlsx", ["u"])
ipc = leer("ipc_general.xlsx", ["ipc"])
pob = leer("poblacion_hombres_y_mujeres.xlsx", ["pob_h", "pob_m"])


def anual(d, s):
    """Promedio anual de los trimestres móviles que terminan en mar, jun, sep y dic (trimestres calendario)."""
    x = d[d.mes.isin([3, 6, 9, 12])].groupby("anio").u.agg(["mean", "size"]).reset_index()
    x.columns = ["anio", "desempleo", "n_trim"]; x["sexo"] = s
    return x


des = pd.concat([anual(um, 1), anual(uh, 0)])
ia = ipc.groupby("anio").ipc.mean()
infl = (100 * (ia / ia.shift(1) - 1)).rename("inflacion").reset_index()            # var. del IPC promedio anual
iond = ipc[ipc.mes.isin([10, 11, 12])].groupby("anio").ipc.mean().rename("ipc_ond").reset_index()
pl = pd.concat([pob[["anio", "pob_h"]].rename(columns={"pob_h": "poblacion"}).assign(sexo=0),
                pob[["anio", "pob_m"]].rename(columns={"pob_m": "poblacion"}).assign(sexo=1)])
esi = pd.read_csv(RES / "esi_por_anio.csv")

b = (des[des.anio.between(2010, 2024)]
     .merge(infl, on="anio").merge(iond, on="anio")
     .merge(esi[["anio", "sexo", "salario_hora_nominal", "educ_sup"]], on=["anio", "sexo"])
     .merge(pl, on=["anio", "sexo"], how="left"))
# Deflactado con el IPC promedio oct-dic (trimestre de la ESI); base 2023 = 100 -> pesos de 2023
b["salario_real_hora"] = b.salario_hora_nominal * 100 / b.ipc_ond
b["inflacion_x_sexo"] = b.inflacion * b.sexo
b["educ_sup_x_sexo"] = b.educ_sup * b.sexo
b["d2019"] = (b.anio == 2019).astype(int)
b["dpandemia"] = b.anio.isin([2020, 2021]).astype(int)
b = b.sort_values(["sexo", "anio"])
g = b.groupby("sexo")
b["d_desempleo"] = g.desempleo.diff()
b["d_inflacion"] = g.inflacion.diff()
b["d_educ_sup"] = g.educ_sup.diff()
b["d_salario"] = 100 * g.salario_real_hora.transform(lambda s: np.log(s).diff())
b["d_salario_nominal"] = 100 * g.salario_hora_nominal.transform(lambda s: np.log(s).diff())
b["d_inflacion_x_sexo"] = b.d_inflacion * b.sexo

cols = ["anio", "sexo", "desempleo", "inflacion", "inflacion_x_sexo", "educ_sup", "educ_sup_x_sexo",
        "salario_real_hora", "salario_hora_nominal", "d2019", "dpandemia", "poblacion",
        "d_desempleo", "d_inflacion", "d_inflacion_x_sexo", "d_educ_sup", "d_salario", "d_salario_nominal"]
assert len(b) == 30 and (b.n_trim == 4).all(), "La base debe tener 30 filas con 4 trimestres por año"
b[cols].round(4).to_csv(RES / "base_apilada_2010_2024.csv", index=False)
pd.set_option("display.width", 200)
print(b[["anio", "sexo", "desempleo", "inflacion", "educ_sup", "salario_real_hora"]].round(2).to_string(index=False))
print("Guardado:", RES / "base_apilada_2010_2024.csv")
