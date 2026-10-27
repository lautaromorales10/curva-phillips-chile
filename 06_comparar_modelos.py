"""Paso 6 (objetivo específico 5). ¿Qué variable explica mejor el desempleo: la inflación o la
composición educacional? Comparación de modelos anidados mediante R2 ajustado, AIC, BIC, R2 parcial
y pruebas F conjuntas con valor p por bootstrap salvaje restringido.
Salida: resultados/comparacion_modelos.json"""
import json
import numpy as np, pandas as pd
from config import RES, B, SEMILLA
from econometria import ols, wild_F

rng = np.random.default_rng(SEMILLA)
b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
for a in (2020, 2021, 2022):
    b[f"d{a}"] = (b.anio == a).astype(int)
b["d_educ_sup_x_sexo"] = b.d_educ_sup * b.sexo


def especificacion(forma, control_2022=False):
    if forma == "dif":
        y, base = "d_desempleo", ["sexo", "d2019", "d2020", "d2021"] + (["d2022"] if control_2022 else [])
        inf, edu = ["d_inflacion", "d_inflacion_x_sexo"], ["d_educ_sup", "d_educ_sup_x_sexo"]
    else:
        y, base = "desempleo", ["sexo", "d2019", "dpandemia"]
        inf, edu = ["inflacion", "inflacion_x_sexo"], ["educ_sup", "educ_sup_x_sexo"]
    return y, base, inf, edu


def restriccion(xs, vs):
    R = np.zeros((len(vs), len(xs) + 1))
    for i, v in enumerate(vs):
        R[i, xs.index(v) + 1] = 1
    return R


def comparar(forma, control_2022=False, bootstrap=True):
    y, base, inf, edu = especificacion(forma, control_2022)
    d = b.dropna(subset=[y] + inf + edu)
    Y, g = d[y].values, d.anio.values
    mods = {"Base (solo controles)": base, "Base + inflación": base + inf,
            "Base + educación": base + edu, "Base + inflación + educación": base + inf + edu}
    filas, rss0 = [], None
    for nom, xs in mods.items():
        X = np.column_stack([np.ones(len(d))] + [d[x].values for x in xs])
        r = ols(Y, X, cluster=g); rss = float(r["resid"] @ r["resid"])
        rss0 = rss if rss0 is None else rss0
        fila = dict(modelo=nom, k=r["k"], r2a=r["r2a"], aic=r["aic"], bic=r["bic"],
                    r2_parcial=1 - rss / rss0, F=None, p=None)
        if bootstrap and nom in ("Base + inflación", "Base + educación"):
            fila["F"], fila["p"] = wild_F(X, Y, g, restriccion(xs, inf if "inflación" in nom else edu), rng, B)
        filas.append(fila)
    extra = {}
    if bootstrap:      # aporte de cada bloque en el modelo completo, controlando por el otro
        xs = base + inf + edu
        X = np.column_stack([np.ones(len(d))] + [d[x].values for x in xs])
        for nom, vs in (("inflación | educación", inf), ("educación | inflación", edu)):
            extra[nom] = wild_F(X, Y, g, restriccion(xs, vs), rng, B)
    return pd.DataFrame(filas), extra, len(d)


out = {}
pd.set_option("display.width", 200)
for forma in ("dif", "niv"):
    t, extra, n = comparar(forma)
    out[forma] = dict(tabla=t.to_dict("records"), extra=extra, n=n)
    print(f"\n===== {'PRIMERAS DIFERENCIAS' if forma == 'dif' else 'NIVELES'} (n = {n}) =====")
    print(t.drop(columns="k").round(3).to_string(index=False))
    for k, (F, p) in extra.items():
        print(f"  Aporte de {k}: F = {F:.2f}, p bootstrap = {p:.3f}")

# Robustez: diferencias controlando por 2022 (solo criterios de ajuste)
t, _, n = comparar("dif", control_2022=True, bootstrap=False)
out["dif_control_2022"] = dict(tabla=t.to_dict("records"), n=n)
print(f"\n===== DIFERENCIAS CONTROLANDO POR 2022 (n = {n}) =====")
print(t.drop(columns=["k", "F", "p"]).round(3).to_string(index=False))

# Correlación entre las series de hombres y mujeres (sección 2.8.2): justifica usar la estructura
# apilada en lugar de incluir ambas series como regresores separados
print("\n===== Correlación entre las series de hombres y mujeres (niveles) =====")
corr = {}
for v in ("desempleo", "educ_sup", "salario_real_hora"):
    w = b.pivot(index="anio", columns="sexo", values=v)
    corr[v] = float(w[0].corr(w[1]))
    print(f"{v:20s} {corr[v]:.3f}")
out["correlaciones_sexo"] = corr

json.dump(out, open(RES / "comparacion_modelos.json", "w"), indent=1, default=float)
print("\nGuardado:", RES / "comparacion_modelos.json")
