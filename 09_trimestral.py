"""Paso 9. Curva de Phillips con datos trimestrales (prueba de robustez con una muestra más grande).
La tasa de desocupación por sexo y la inflación están disponibles con frecuencia mensual, por lo que la
especificación en diferencias puede estimarse con trimestres calendario (enero-marzo, abril-junio,
julio-septiembre y octubre-diciembre). La composición educacional y el ingreso por hora provienen de la
ESI, que es anual, por lo que no forman parte de esta estimación.
  - Desocupación: trimestre móvil que termina en marzo, junio, septiembre y diciembre (sin traslape).
  - Inflación: variación porcentual del IPC promedio del trimestre respecto del mismo trimestre del año
    anterior (análoga a la definición anual).
  - Variaciones interanuales (respecto del mismo trimestre del año anterior), que eliminan la estacionalidad
    de las series no desestacionalizadas: 2011T1-2025T4, 60 trimestres por sexo.
Inferencia: errores agrupados por año y bootstrap salvaje restringido por año, como en la estimación anual;
como verificación, errores de Driscoll y Kraay (1998) con 4 rezagos, que consideran la autocorrelación que
genera el traslape de las variaciones interanuales.
Salidas: resultados/base_trimestral_2010_2025.csv, resultados/trimestral.json,
         resultados/tablas/tabla_trimestral.tex"""
import json
import numpy as np, pandas as pd
from scipy import stats
from config import DATOS, RES, TEX, B, SEMILLA, ANIO_FIN
from econometria import ols, wild_t, driscoll_kraay, durbin_watson

REZAGOS_DK = 4


def leer(archivo):
    d = pd.read_excel(DATOS / archivo, sheet_name="Cuadro", header=None, skiprows=3).iloc[:, :2]
    d.columns = ["fecha", "valor"]
    d["fecha"] = pd.to_datetime(d["fecha"], errors="coerce")
    return d.dropna(subset=["fecha"]).set_index("fecha").valor.astype(float)


ipc = leer("ipc_general.xlsx")
ipc_t = ipc.groupby(pd.PeriodIndex(ipc.index, freq="Q")).mean()
inflacion = 100 * (ipc_t / ipc_t.shift(4) - 1)

partes = []
for archivo, sexo in (("desocupacion_hombres.xlsx", 0), ("desocupacion_mujeres.xlsx", 1)):
    u = leer(archivo)
    u = u[u.index.month.isin([3, 6, 9, 12])]
    u.index = pd.PeriodIndex(u.index, freq="Q")
    d = pd.DataFrame({"desempleo": u, "inflacion": inflacion.reindex(u.index)})
    d["d_desempleo"] = d.desempleo - d.desempleo.shift(4)
    d["d_inflacion"] = d.inflacion - d.inflacion.shift(4)
    d["sexo"] = sexo
    partes.append(d)
b = pd.concat(partes)
b = b[(b.index >= pd.Period("2010Q1")) & (b.index <= pd.Period(f"{ANIO_FIN}Q4"))].copy()
b["anio"], b["trimestre"] = b.index.year, b.index.quarter
b["periodo"] = b.anio * 4 + b.trimestre                        # índice trimestral consecutivo (Driscoll-Kraay)
assert len(b) == 8 * (ANIO_FIN - 2009), "La base trimestral debe tener cuatro trimestres por año y sexo"
b.reset_index(drop=True)[["anio", "trimestre", "sexo", "desempleo", "inflacion", "d_desempleo", "d_inflacion"]] \
    .round(4).to_csv(RES / f"base_trimestral_2010_{ANIO_FIN}.csv", index=False)

d = b.dropna(subset=["d_desempleo", "d_inflacion"]).copy()
d["d_inflacion_x_sexo"] = d.d_inflacion * d.sexo
IND = {"2019T4": [(2019, 4)],
       "pandemia": [(a, t) for a in (2020, 2021) for t in (1, 2, 3, 4)],
       "2022": [(2022, t) for t in (1, 2, 3, 4)]}
for clave, trims in IND.items():
    for a, t in trims:
        d[f"d{a}t{t}"] = ((d.anio == a) & (d.trimestre == t)).astype(int)
nombres = lambda clave: [f"d{a}t{t}" for a, t in IND[clave]]
BASE = ["d_inflacion", "d_inflacion_x_sexo", "sexo"] + nombres("2019T4") + nombres("pandemia")
MODELOS = {"T1": (d, BASE),                                          # especificación base
           "T2": (d, BASE + nombres("2022")),                        # controla por los trimestres de 2022
           "T3": (d[~d.anio.isin([2022, 2023])], BASE)}              # excluye 2022 y 2023

rng = np.random.default_rng(SEMILLA)
out = {}
for m, (datos, xs) in MODELOS.items():
    X = np.column_stack([np.ones(len(datos))] + [datos[x].values for x in xs])
    Y, g = datos.d_desempleo.values, datos.anio.values
    r = ols(Y, X, cluster=g)
    dk = driscoll_kraay(Y, X, datos.periodo.values, REZAGOS_DK)
    k = X.shape[1]
    Rh = np.zeros(k); Rh[1] = 1
    Ri = np.zeros(k); Ri[2] = 1
    Rm = Rh + Ri
    se_dk = lambda R: float(np.sqrt(R @ dk["V"] @ R))
    p_dk = lambda R: float(2 * stats.t.sf(abs(R @ dk["b"]) / se_dk(R), dk["df"]))
    out[m] = dict(
        n=int(r["n"]), grupos=int(len(np.unique(g))), r2a=float(r["r2a"]),
        dw=[durbin_watson(r["resid"][datos.sexo.values == s]) for s in (0, 1)],
        coefs={v: dict(b=float(r["b"][i]), se=float(r["se"][i])) for i, v in enumerate(["const"] + xs[:3])},
        h=dict(b=float(r["b"][1]), p_wb=wild_t(X, Y, g, Rh, rng, B)[1], p_dk=p_dk(Rh)),
        m=dict(b=float(r["b"][1] + r["b"][2]), p_wb=wild_t(X, Y, g, Rm, rng, B)[1],
               p_unil=wild_t(X, Y, g, Rm, rng, B, unilateral=True)[1], p_dk=p_dk(Rm)),
        inter=dict(b=float(r["b"][2]), p_wb=wild_t(X, Y, g, Ri, rng, B)[1]),
        const_p=wild_t(X, Y, g, np.eye(k)[0], rng, B)[1], sexo_p=wild_t(X, Y, g, np.eye(k)[3], rng, B)[1])
    o = out[m]
    print(f"({m}) n = {o['n']}, {o['grupos']} años | hombres {o['h']['b']:.3f} (p bootstrap {o['h']['p_wb']:.3f}; "
          f"DK {o['h']['p_dk']:.3f}) | mujeres {o['m']['b']:.3f} (p bootstrap {o['m']['p_wb']:.3f}; unilateral "
          f"{o['m']['p_unil']:.3f}; DK {o['m']['p_dk']:.3f}) | interacción p {o['inter']['p_wb']:.3f} | "
          f"R2 ajustado {o['r2a']:.3f}", flush=True)
out["correlacion"] = {s: float(d[d.sexo == s][["d_desempleo", "d_inflacion"]].corr().iloc[0, 1]) for s in (0, 1)}
out["n_trimestres"] = int(d.periodo.nunique())
json.dump(out, open(RES / "trimestral.json", "w"), indent=1)

# ------------------------------------------------------------------ Tabla LaTeX
def num(x, dec=3):
    s = f"{x:.{dec}f}".replace(".", ",")
    return s.replace("-", "$-$") if s.startswith("-") else s
def est(p): return "$^{***}$" if p < .01 else "$^{**}$" if p < .05 else "$^{*}$" if p < .1 else ""
def pv(p): return f"[{p:.3f}]".replace(".", ",")
cols = ["T1", "T2", "T3"]
fila = lambda etq, f: " & ".join([etq] + [f(m) for m in cols]) + r" \\"
P = {"d_inflacion": lambda m: out[m]["h"]["p_wb"], "d_inflacion_x_sexo": lambda m: out[m]["inter"]["p_wb"],
     "sexo": lambda m: out[m]["sexo_p"], "const": lambda m: out[m]["const_p"]}
L = []
for etq, v in (("Inflación", "d_inflacion"), ("Inflación $\\times$ Mujer", "d_inflacion_x_sexo"),
               ("Mujer", "sexo"), ("Constante", "const")):
    L.append(fila(etq, lambda m: num(out[m]["coefs"][v]["b"]) + est(P[v](m))))
    L.append(fila("", lambda m: f"({out[m]['coefs'][v]['se']:.3f})".replace(".", ",")).replace(r" \\", r" \\[2pt]"))
pie = [fila("Pendiente para hombres: valor $p$ (bootstrap)", lambda m: pv(out[m]["h"]["p_wb"])),
       fila("Pendiente para mujeres", lambda m: num(out[m]["m"]["b"]) + est(out[m]["m"]["p_wb"])),
       fila("\\quad valor $p$ (bootstrap)", lambda m: pv(out[m]["m"]["p_wb"])),
       fila("Valor $p$ Driscoll--Kraay (H / M)",
            lambda m: f"{out[m]['h']['p_dk']:.3f} / {out[m]['m']['p_dk']:.3f}".replace(".", ",")),
       r"\midrule",
       fila("Indicadora estallido (4.\\textsuperscript{o} trim. 2019)", lambda m: "Sí"),
       fila("Indicadoras pandemia (trimestres 2020--2021)", lambda m: "Sí"),
       fila("Indicadoras trimestres 2022", lambda m: "Sí" if m == "T2" else "No"),
       fila("Años excluidos", lambda m: "2022--2023" if m == "T3" else "--"),
       r"\midrule",
       fila("Observaciones", lambda m: str(out[m]["n"])),
       fila("Grupos (años)", lambda m: str(out[m]["grupos"])),
       fila("$R^2$ ajustado", lambda m: num(out[m]["r2a"])),
       fila("Durbin--Watson (H / M)", lambda m: f"{out[m]['dw'][0]:.2f} / {out[m]['dw'][1]:.2f}".replace(".", ","))]
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{5pt}
\caption{Desocupación por sexo e inflación con datos trimestrales, Chile 2011--2025}
\label{tab:trimestral}
\begin{tabular}{lccc}
\toprule
 & (1) & (2) & (3) \\
Variable dependiente & $\Delta_4 U$ & $\Delta_4 U$ & $\Delta_4 U$ \\
\midrule
""" + "\n".join(L) + "\n\\midrule\n" + "\n".join(pie) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Mínimos cuadrados ordinarios sobre la base apilada por sexo con datos trimestrales
(trimestres calendario). $\Delta_4$ denota la variación respecto del mismo trimestre del año anterior; la inflación
corresponde a la variación interanual del IPC promedio de cada trimestre y entra en la regresión como $\Delta_4 \pi$. Entre paréntesis, errores estándar agrupados
por año. La significancia se evalúa mediante bootstrap salvaje restringido por año con pesos de Rademacher y 9.999
réplicas \citep{cameron2008}. Los valores $p$ de Driscoll--Kraay usan 4 rezagos \citep{driscoll1998}. Los coeficientes de
las variables indicadoras se omiten. $^{***}p<0{,}01$; $^{**}p<0{,}05$; $^{*}p<0{,}1$. Fuente: elaboración propia a partir
de INE y Banco Central de Chile.
\end{minipage}
\end{table}
"""
(TEX / "tabla_trimestral.tex").write_text(tex, encoding="utf-8")
print("Guardado:", RES / "trimestral.json", "y", TEX / "tabla_trimestral.tex")
