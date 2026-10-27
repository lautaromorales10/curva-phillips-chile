"""Paso 8. Pruebas de robustez adicionales de la Curva de Phillips (especificación en diferencias):
  A. Rezagos de la inflación (la respuesta del desempleo podría no ser contemporánea).
  B. Exclusión de un año a la vez (influencia de cada año en la pendiente estimada).
  C. Brecha de desempleo respecto de su tendencia (filtro de Hodrick-Prescott) en lugar de diferencias.
Salidas: resultados/robustez.json, resultados/tablas/tabla_robustez.tex, resultados/figuras/fig6_exclusion_anual.png"""
import json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import RES, FIG, TEX, B, SEMILLA
from econometria import ols, wild_t, durbin_watson

b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
for a in range(2019, 2024):
    b[f"d{a}"] = (b.anio == a).astype(int)
g = b.groupby("sexo")
b["d_inflacion_l1"] = g.d_inflacion.shift(1)                      # rezago de la variación de la inflación
b["d_inflacion_l1_x_sexo"] = b.d_inflacion_l1 * b.sexo
DD = ["d2019", "d2020", "d2021"]
C = ["d_inflacion", "d_inflacion_x_sexo"]
L = ["d_inflacion_l1", "d_inflacion_l1_x_sexo"]


def estimar(datos, y, xs, pruebas):
    """MCO con errores agrupados por año; valor p de cada combinación lineal por bootstrap salvaje.
    Las indicadoras que quedan en cero (por ejemplo, al excluir su año) se eliminan."""
    d = datos.dropna(subset=[y] + xs)
    xs = [x for x in xs if d[x].abs().sum() > 0]
    X = np.column_stack([np.ones(len(d))] + [d[x].values for x in xs])
    Y, G = d[y].values.astype(float), d.anio.values
    r = ols(Y, X, cluster=G)
    res = dict(n=int(r["n"]), r2a=float(r["r2a"]),
               dw=[durbin_watson(r["resid"][d.sexo.values == s]) for s in (0, 1)])
    for etq, comb in pruebas.items():
        R = np.zeros(X.shape[1])
        for v in comb:
            R[xs.index(v) + 1] = 1
        res[etq] = dict(b=float(R @ r["b"]), p=wild_t(X, Y, G, R, np.random.default_rng(SEMILLA), B)[1])
    return res


out = {}
# ------------------------------------------------------------------ A. Rezagos
pr_l = {"rez_H": ["d_inflacion_l1"], "rez_M": L}
pr_cl = {"con_H": ["d_inflacion"], "con_M": C, "rez_H": ["d_inflacion_l1"], "rez_M": L,
         "sum_H": ["d_inflacion", "d_inflacion_l1"], "sum_M": C + L}
out["A1_rezago"] = estimar(b, "d_desempleo", L + ["sexo"] + DD, pr_l)
out["A2_contemp_rezago"] = estimar(b, "d_desempleo", C + L + ["sexo"] + DD, pr_cl)
out["A3_contemp_rezago_2022"] = estimar(b, "d_desempleo", C + L + ["sexo"] + DD + ["d2022"], pr_cl)

# ------------------------------------------------------------------ B. Exclusión de un año a la vez
pr_c = {"con_H": ["d_inflacion"], "con_M": C}
base = estimar(b, "d_desempleo", C + ["sexo"] + DD, pr_c)
jk = {}
for a in range(2011, 2025):
    jk[a] = estimar(b[b.anio != a], "d_desempleo", C + ["sexo"] + DD, pr_c)
out["B_base"] = base
out["B_exclusion"] = {str(a): v for a, v in jk.items()}


# ------------------------------------------------------------------ C. Brecha de desempleo (filtro HP)
def filtro_hp(y, lam):
    """Tendencia de Hodrick-Prescott: minimiza la suma de desviaciones al cuadrado más lam veces
    la suma de las segundas diferencias de la tendencia al cuadrado."""
    n = len(y); K = np.zeros((n - 2, n))
    for i in range(n - 2):
        K[i, i:i + 3] = [1, -2, 1]
    return np.linalg.solve(np.eye(n) + lam * K.T @ K, y)


def brecha(datos, lam, interpolar=True):
    """Brecha = desempleo observado - tendencia HP. Como el filtro es bilateral, el salto de 2020 desplazaría
    hacia arriba la tendencia de los años adyacentes; por ello, por defecto, la tendencia se estima tras
    reemplazar 2020 y 2021 por una interpolación lineal entre 2019 y 2022 (esos años quedan cubiertos por
    sus indicadoras en la regresión)."""
    x = pd.Series(np.nan, index=datos.index)
    for s in (0, 1):
        m = datos.sexo == s
        u, anios = datos.loc[m, "desempleo"].values.astype(float), datos.loc[m, "anio"].values
        ui = u.copy()
        if interpolar:
            k = np.isin(anios, [2020, 2021])
            ui[k] = np.interp(anios[k], anios[~k], u[~k])
        x[m] = u - filtro_hp(ui, lam)
    return x


pr_n = {"con_H": ["inflacion"], "con_M": ["inflacion", "inflacion_x_sexo"]}
XB = ["inflacion", "inflacion_x_sexo", "sexo"] + DD
for lam in (6.25, 100):
    b["brecha"] = brecha(b, lam)
    out[f"C_hp{lam}"] = estimar(b, "brecha", XB, pr_n)
    out[f"C_hp{lam}_2022"] = estimar(b, "brecha", XB + ["d2022"], pr_n)
    out[f"C_hp{lam}_serie"] = b[["anio", "sexo", "brecha"]].round(4).values.tolist()
    b["brecha"] = brecha(b, lam, interpolar=False)          # tendencia con la serie observada (comparación)
    out[f"C_hp{lam}_obs"] = estimar(b, "brecha", XB, pr_n)
    out[f"C_hp{lam}_obs_serie"] = b[["anio", "sexo", "brecha"]].round(4).values.tolist()

json.dump(out, open(RES / "robustez.json", "w"), indent=1)
for k, v in out.items():
    if k.startswith("B_exclusion") or k.endswith("_serie"):
        continue
    print(k, {kk: (round(vv["b"], 3), round(vv["p"], 3)) for kk, vv in v.items() if isinstance(vv, dict)},
          f"n={v['n']} R2aj={v['r2a']:.3f}")
for a, v in jk.items():
    print(f"sin {a}: hombres {v['con_H']['b']:.3f} (p {v['con_H']['p']:.3f}) | mujeres {v['con_M']['b']:.3f} (p {v['con_M']['p']:.3f})")


# ------------------------------------------------------------------ Tabla LaTeX
def num(x, d=3):
    s = f"{x:.{d}f}".replace(".", ",")
    return s.replace("-", "$-$") if s.startswith("-") else s


def fila(etq, r, kH, kM, extra=True):
    c = [etq, num(r[kH]["b"]), "[" + num(r[kH]["p"]) + "]", num(r[kM]["b"]), "[" + num(r[kM]["p"]) + "]"]
    c += [num(r["r2a"]), str(r["n"])] if extra else ["", ""]
    return " & ".join(c) + r" \\"


A1, A2, A3 = out["A1_rezago"], out["A2_contemp_rezago"], out["A3_contemp_rezago_2022"]
C1, C2, C3 = out["C_hp6.25"], out["C_hp6.25_2022"], out["C_hp6.25_obs"]
filas = [r"\multicolumn{7}{l}{\textit{Panel A. Rezagos de la inflación (variables en primeras diferencias)}} \\",
         r"\quad (A1) Solo rezago: $\Delta\pi_{t-1}$ & " + " & ".join(fila("", A1, "rez_H", "rez_M").split(" & ")[1:]),
         r"\quad (A2) Contemporánea y rezago: $\Delta\pi_t$ & " + " & ".join(fila("", A2, "con_H", "con_M").split(" & ")[1:]),
         fila(r"\qquad $\Delta\pi_{t-1}$", A2, "rez_H", "rez_M", False),
         fila(r"\qquad Efecto acumulado", A2, "sum_H", "sum_M", False),
         r"\quad (A3) A2 con indicadora 2022: $\Delta\pi_t$ & " + " & ".join(fila("", A3, "con_H", "con_M").split(" & ")[1:]),
         fila(r"\qquad $\Delta\pi_{t-1}$", A3, "rez_H", "rez_M", False),
         fila(r"\qquad Efecto acumulado", A3, "sum_H", "sum_M", False),
         r"\midrule",
         r"\multicolumn{7}{l}{\textit{Panel B. Brecha de desempleo respecto de su tendencia (filtro de Hodrick--Prescott)}} \\",
         r"\quad (B1) Brecha sobre inflación: $\pi_t$ & " + " & ".join(fila("", C1, "con_H", "con_M").split(" & ")[1:]),
         r"\quad (B2) B1 con indicadora 2022: $\pi_t$ & " + " & ".join(fila("", C2, "con_H", "con_M").split(" & ")[1:]),
         r"\quad (B3) Tendencia con la serie observada: $\pi_t$ & " + " & ".join(fila("", C3, "con_H", "con_M").split(" & ")[1:])]
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{5pt}
\caption{Pruebas de robustez: rezagos de la inflación y brecha de desempleo}
\label{tab:robustez}
\begin{tabular}{lcccccc}
\toprule
 & \multicolumn{2}{c}{Hombres} & \multicolumn{2}{c}{Mujeres} & & \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Especificación y coeficiente & Coef. & Valor $p$ & Coef. & Valor $p$ & $R^2$ ajust. & $n$ \\
\midrule
""" + "\n".join(filas) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Mínimos cuadrados ordinarios sobre la base apilada por sexo, con errores agrupados por año; valores $p$
(entre corchetes) obtenidos mediante bootstrap salvaje restringido por año con 9.999 réplicas. Todas las especificaciones incluyen
la variable de sexo y las indicadoras de 2019, 2020 y 2021. Panel A: la variable dependiente es la variación de la tasa de
desocupación; el efecto acumulado corresponde a la suma de los coeficientes contemporáneo y rezagado. Panel B: la variable
dependiente es la desviación de la tasa de desocupación de cada sexo respecto de su tendencia de Hodrick--Prescott
($\lambda = 6{,}25$). En B1 y B2, la tendencia se estima tras reemplazar los valores de 2020 y 2021 por una interpolación
lineal entre 2019 y 2022, para que el salto de la pandemia no desplace la tendencia de los años adyacentes; B3 la estima
con la serie observada. Elaboración propia a partir de INE y Banco Central de Chile.
\end{minipage}
\end{table}
"""
open(TEX / "tabla_robustez.tex", "w", encoding="utf-8").write(tex)

# ------------------------------------------------------------------ Figura: exclusión de un año a la vez
C_H, C_M, INK, INK2, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.family": "DejaVu Serif", "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "savefig.dpi": 300, "savefig.bbox": "tight"})
anios = list(jk.keys())
from matplotlib.lines import Line2D
fig, ax = plt.subplots(figsize=(6.3, 3.6))
ax.axhline(0, color=INK2, lw=0.8)
for clave, col, etq, dx in [("con_H", C_H, "Hombres", -0.17), ("con_M", C_M, "Mujeres", 0.17)]:
    ax.axhline(base[clave]["b"], color=col, lw=1, ls=(0, (4, 3)), zorder=1)
    for a in anios:
        sig = jk[a][clave]["p"] < 0.10
        ax.plot(a + dx, jk[a][clave]["b"], marker="o", ms=7, mew=1.6, color=col,
                mfc=col if sig else "white", mec=col, ls="none", zorder=3)
    ax.annotate(etq, (2024.6, base[clave]["b"]), xytext=(4, 0), textcoords="offset points",
                va="center", fontsize=9, color=INK, annotation_clip=False)
h = [Line2D([], [], marker="o", ms=7, color=C_H, ls="none", label="Hombres"),
     Line2D([], [], marker="o", ms=7, color=C_M, ls="none", label="Mujeres"),
     Line2D([], [], marker="o", ms=7, color=INK2, ls="none", label="Significativa al 10%"),
     Line2D([], [], marker="o", ms=7, mfc="white", mec=INK2, mew=1.6, ls="none", label="No significativa")]
ax.set_xticks(anios); ax.set_xticklabels([str(a) for a in anios], rotation=45, fontsize=8.5)
ax.set_xlim(2010.4, 2024.6)
ax.set_xlabel("Año excluido de la estimación")
ax.set_ylabel("Pendiente estimada\n(pp. de desempleo por pp. de inflación)")
ax.legend(handles=h, frameon=False, fontsize=8, ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0))
fig.savefig(FIG / "fig6_exclusion_anual.png"); plt.close(fig)
print("Guardado:", RES / "robustez.json", TEX / "tabla_robustez.tex", FIG / "fig6_exclusion_anual.png")
