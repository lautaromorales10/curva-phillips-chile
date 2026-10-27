"""Paso 5. Estimación de la Curva de Phillips y de los determinantes del desempleo.
MCO sobre la base apilada, errores agrupados por año e inferencia por bootstrap salvaje restringido
por año (pesos de Rademacher, 9.999 réplicas), adecuado con solo 15 grupos (Cameron et al., 2008).
Salida: resultados/resultados_modelos.json"""
import json
import numpy as np, pandas as pd
from config import RES, B, SEMILLA
from econometria import ols, wald, wild_t, durbin_watson

b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
for a in (2020, 2021, 2022):
    b[f"d{a}"] = (b.anio == a).astype(int)
b["desempleo_x_sexo"] = b.desempleo * b.sexo
b["d_desempleo_x_sexo"] = b.d_desempleo * b.sexo
b["d_educ_sup_x_sexo"] = b.d_educ_sup * b.sexo
for a in (2019, 2020, 2021):
    b[f"d{a}_x_sexo"] = b[f"d{a}"] * b.sexo
DD = ["d2019", "d2020", "d2021"]   # en diferencias, 2020 y 2021 por separado (efectos de signo opuesto)

# Columnas (1) a (6) de la Tabla de resultados: (variable dependiente, regresores). El primer regresor
# es la pendiente para hombres y el segundo su interacción con sexo (pendiente mujeres = suma de ambos).
# M1-M2: niveles; M3-M5: primeras diferencias; M6: curva salarial (ecuación 5.4).
MODELOS = {
    "M1": ("desempleo",   ["inflacion", "inflacion_x_sexo", "sexo", "d2019", "dpandemia"]),
    "M2": ("desempleo",   ["inflacion", "inflacion_x_sexo", "educ_sup", "educ_sup_x_sexo", "sexo", "d2019", "dpandemia"]),
    "M3": ("d_desempleo", ["d_inflacion", "d_inflacion_x_sexo", "sexo"] + DD),
    "M4": ("d_desempleo", ["d_inflacion", "d_inflacion_x_sexo", "d_educ_sup", "d_educ_sup_x_sexo", "sexo"] + DD),
    "M5": ("d_desempleo", ["d_inflacion", "d_inflacion_x_sexo", "sexo"] + DD + ["d2022"]),
    "M6": ("d_salario",   ["desempleo", "desempleo_x_sexo", "sexo"] + DD),
}
INDICADORAS = {"d2019", "dpandemia", "d2020", "d2021", "d2022"}
rng = np.random.default_rng(SEMILLA)


def preparar(datos, y, xs):
    d = datos.dropna(subset=[y] + xs)
    X = np.column_stack([np.ones(len(d))] + [d[x].values for x in xs])
    return d, X, d[y].values.astype(float), d.anio.values


out = {}
for nombre, (y, xs) in MODELOS.items():
    d, X, Y, g = preparar(b, y, xs)
    r = ols(Y, X, cluster=g)
    coefs = {}
    for i, v in enumerate(["const"] + xs):
        p = None
        if v not in INDICADORAS:     # las indicadoras se identifican con un solo año: solo controles
            R = np.zeros(X.shape[1]); R[i] = 1
            p = wild_t(X, Y, g, R, rng, B)[1]
        coefs[v] = dict(b=float(r["b"][i]), se=float(r["se"][i]), p_cl=float(r["p"][i]), p_wb=p)
    R = np.zeros(X.shape[1]); R[1] = R[2] = 1
    pm = dict(b=float(r["b"][1] + r["b"][2]), p_cl=float(wald(r, R)[1]), p_wb=wild_t(X, Y, g, R, rng, B)[1])
    dw = [durbin_watson(r["resid"][d.sexo.values == s]) for s in (0, 1)]
    out[nombre] = dict(y=y, xs=xs, coefs=coefs, mujeres=pm, n=int(r["n"]), r2a=float(r["r2a"]), dw=dw)
    print(f"\n({nombre}) {y}  n = {r['n']}  R2 ajustado = {r['r2a']:.3f}  DW (H/M) = {dw[0]:.2f} / {dw[1]:.2f}")
    for v, c in coefs.items():
        pw = "   --" if c["p_wb"] is None else f"{c['p_wb']:.3f}"
        print(f"   {v:22s} {c['b']:8.3f}  ee {c['se']:.3f}  p bootstrap {pw}")
    print(f"   {'pendiente mujeres':22s} {pm['b']:8.3f}  p bootstrap {pm['p_wb']:.3f}", flush=True)

# Sensibilidad del modelo M3: excluir episodios
print("\n===== Sensibilidad (modelo M3) =====")
sens = {}
for etq, excl, ds in [("sin_2022_2023", [2022, 2023], DD), ("sin_2020_2021", [2020, 2021], ["d2019"])]:
    xs = ["d_inflacion", "d_inflacion_x_sexo", "sexo"] + ds
    d, X, Y, g = preparar(b[~b.anio.isin(excl)], "d_desempleo", xs)
    r = ols(Y, X, cluster=g)
    Rh = np.zeros(X.shape[1]); Rh[1] = 1; Rm = Rh.copy(); Rm[2] = 1
    sens[etq] = dict(n=int(r["n"]), h=float(r["b"][1]), p_h=wild_t(X, Y, g, Rh, rng, B)[1],
                     m=float(r["b"][1] + r["b"][2]), p_m=wild_t(X, Y, g, Rm, rng, B)[1])
    s = sens[etq]
    print(f"{etq}: n = {s['n']} | hombres {s['h']:.3f} (p {s['p_h']:.3f}) | mujeres {s['m']:.3f} (p {s['p_m']:.3f})")
out["sens"] = sens

# Prueba unilateral (H1: pendiente < 0), como complemento de la bilateral. Semilla propia por prueba.
print("\n===== Valores p unilaterales (H1: pendiente negativa) =====")
unil = {}
for nombre in ("M3", "M4", "M5"):
    y, xs = MODELOS[nombre]
    d, X, Y, g = preparar(b, y, xs)
    Rh = np.zeros(X.shape[1]); Rh[1] = 1; Rm = Rh.copy(); Rm[2] = 1
    unil[nombre] = {s: wild_t(X, Y, g, R, np.random.default_rng(SEMILLA), B, unilateral=True)[1]
                    for s, R in (("hombres", Rh), ("mujeres", Rm))}
    print(f"({nombre}) hombres p = {unil[nombre]['hombres']:.3f} | mujeres p = {unil[nombre]['mujeres']:.3f}")
out["unilateral"] = unil

# Sensibilidad de los modelos en niveles al tratamiento de la pandemia
print("\n===== Sensibilidad de los modelos en niveles =====")
sens_niv = {}
for nombre in ("M1", "M2"):
    y, xs0 = MODELOS[nombre]
    xs0 = [x for x in xs0 if x != "dpandemia"]
    for etq, datos, extra in (("solo_2020", b, ["d2020"]), ("sin_2020_2021", b[~b.anio.isin([2020, 2021])], [])):
        d, X, Y, g = preparar(datos, y, xs0 + extra)
        r = ols(Y, X, cluster=g)
        Rh = np.zeros(X.shape[1]); Rh[1] = 1; Rm = Rh.copy(); Rm[2] = 1
        res = dict(h=float(r["b"][1]), p_h=wild_t(X, Y, g, Rh, np.random.default_rng(SEMILLA), B)[1],
                   m=float(r["b"][1] + r["b"][2]), p_m=wild_t(X, Y, g, Rm, np.random.default_rng(SEMILLA), B)[1])
        sens_niv[f"{nombre}_{etq}"] = res
        print(f"({nombre}) {etq}: hombres {res['h']:.3f} (p {res['p_h']:.3f}) | mujeres {res['m']:.3f} (p {res['p_m']:.3f})")
out["sens_niveles"] = sens_niv


def pendientes(datos, y, xs, etq):
    """Pendiente de hombres (primer regresor) y de mujeres (suma con su interacción), con valor p por
    bootstrap salvaje; cada prueba usa su propia semilla para que el resultado no dependa del orden."""
    d, X, Y, g = preparar(datos, y, xs)
    r = ols(Y, X, cluster=g)
    Rh = np.zeros(X.shape[1]); Rh[1] = 1; Rm = Rh.copy(); Rm[2] = 1
    res = dict(n=int(r["n"]), r2a=float(r["r2a"]),
               dw=[durbin_watson(r["resid"][d.sexo.values == s]) for s in (0, 1)],
               h=float(r["b"][1]), p_h=wild_t(X, Y, g, Rh, np.random.default_rng(SEMILLA), B)[1],
               m=float(r["b"][1] + r["b"][2]), p_m=wild_t(X, Y, g, Rm, np.random.default_rng(SEMILLA), B)[1])
    print(f"{etq:42s} hombres {res['h']:7.3f} (p {res['p_h']:.3f}) | mujeres {res['m']:7.3f} (p {res['p_m']:.3f})"
          f" | R2a {res['r2a']:.3f} | DW {res['dw'][0]:.2f} / {res['dw'][1]:.2f}")
    return res


# Sensibilidad del modelo M3: indicadoras de eventos diferenciadas por sexo. Con indicadoras comunes, los
# residuos de hombres y mujeres de 2019, 2020 y 2021 suman cero; al interactuarlas con sexo, ambos residuos
# de esos años se anulan, lo que equivale a excluirlos de la estimación de las pendientes.
print("\n===== Sensibilidad (modelo M3): indicadoras diferenciadas por sexo =====")
y, xs = MODELOS["M3"]
out["sens_indicadoras_sexo"] = pendientes(b, y, xs + [f"d{a}_x_sexo" for a in (2019, 2020, 2021)],
                                          "M3 con indicadoras x sexo")

# Curva de Phillips salarial: especificaciones alternativas (ingreso real o nominal; desempleo en
# nivel, como en Phillips (1958), o en primeras diferencias)
print("\n===== Curva salarial: especificaciones alternativas =====")
m6 = out["M6"]      # la especificación base es la columna 6: se reutiliza para que ambas tablas coincidan
sal = {"real_nivel": dict(n=m6["n"], r2a=m6["r2a"], dw=m6["dw"],
                          h=m6["coefs"]["desempleo"]["b"], p_h=m6["coefs"]["desempleo"]["p_wb"],
                          m=m6["mujeres"]["b"], p_m=m6["mujeres"]["p_wb"])}
for etq, yv, forma, xv in (("nominal", "d_salario_nominal", "nivel", "desempleo"),
                           ("real", "d_salario", "dif", "d_desempleo"),
                           ("nominal", "d_salario_nominal", "dif", "d_desempleo")):
    sal[f"{etq}_{forma}"] = pendientes(b, yv, [xv, f"{xv}_x_sexo", "sexo"] + DD, f"salario {etq}, desempleo en {forma}")
out["salarial"] = sal
json.dump(out, open(RES / "resultados_modelos.json", "w"), indent=1)
print("\nGuardado:", RES / "resultados_modelos.json")
