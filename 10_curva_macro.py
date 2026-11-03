"""Paso 10. Parte I: Curva de Phillips con datos agregados nacionales (la inflación es la variable dependiente).
Datos: series mensuales del Banco Central de Chile (datos/bcch_mensual_1995_2026.csv; fuentes y códigos en
datos/FUENTES_bcch_mensual.md) y el IPC empalmado de datos/ipc_general.xlsx.

Muestra principal: trimestres calendario 2010T1-2025T4 (64 observaciones).
  pi      inflación: variación del IPC promedio del trimestre respecto del mismo trimestre del año anterior (%)
  pie_a   expectativa formada un año antes: mediana de la Encuesta de Expectativas Económicas (EEE) para la
          inflación en 12 meses, a 11 meses plazo, registrada 11 meses antes (promedio de los meses del trimestre)
  pie_l1  expectativa a 11 meses plazo del trimestre anterior (formulación híbrida)
  u       tasa de desocupación nacional desestacionalizada (trimestre móvil que termina en el último mes)
  brecha  u menos su tendencia de Hodrick-Prescott (lambda = 1.600), estimada hasta 2026T2 con los trimestres
          2020T2-2021T1 interpolados (como en la Parte II)
  dep     depreciación: variación logarítmica en 12 meses del dólar observado (promedio del trimestre)
  pet     variación logarítmica en 12 meses del precio del petróleo Brent (promedio del trimestre)
  indicadoras: estallido social (2019T4) y pandemia (2020T1-2021T4), como controles de nivel.
Inferencia: errores de Newey-West con 4 rezagos (12 con datos mensuales).
Robustez: (i) datos mensuales 2010-2025; (ii) muestra 1997-2025, con la desocupación de la antigua ENE hasta 2009T4
          desestacionalizada con efectos trimestrales y una indicadora de nivel desde 2010 (cambio de encuesta).
Salidas: resultados/base_macro_trimestral.csv, resultados/macro.json, resultados/tablas/tabla_macro*.tex,
         resultados/figuras/fig7_*.png, fig8_*.png"""
import json
import numpy as np, pandas as pd, matplotlib
from scipy import stats
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import DATOS, RES, TEX, FIG
from econometria import newey_west, filtro_hp, adf, kpss, durbin_watson

INI, FIN = "2010Q1", "2025Q4"
POST = "2022Q1"                    # inicio del periodo posterior a la pandemia
L_T, L_M = 4, 12                   # rezagos de Newey-West (trimestral, mensual)

# ------------------------------------------------------------------ Datos mensuales
m = pd.read_csv(DATOS / "bcch_mensual_1995_2026.csv", index_col=0)
m.index = pd.PeriodIndex(m.index, freq="M")
d = pd.read_excel(DATOS / "ipc_general.xlsx", sheet_name="Cuadro", header=None, skiprows=3).iloc[:, :2]
d.columns = ["f", "v"]; d["f"] = pd.to_datetime(d.f, errors="coerce")
ipc = d.dropna().set_index("f").v.astype(float); ipc.index = ipc.index.to_period("M")
m["ipc"] = ipc
m["pi"] = 100 * (m.ipc / m.ipc.shift(12) - 1)
m["pie_a"] = m.eee_11m.shift(11)
m["dep"] = 100 * np.log(m.dolar_obs / m.dolar_obs.shift(12))
m["pet"] = 100 * np.log(m.brent / m.brent.shift(12))
m["u_nsa"] = m.u_ine.combine_first(m.u_ene_antigua)          # antigua ENE hasta feb. 2010, ENE desde mar. 2010


def promedio_q(s):
    s = s.dropna(); return s.groupby(s.index.asfreq("Q")).mean()


def fin_q(s):
    s = s.dropna(); s = s[s.index.month.isin([3, 6, 9, 12])]; return s.set_axis(s.index.asfreq("Q"))


def brecha_hp(u, ini, fin, interp, lam):
    """u menos su tendencia HP; los periodos 'interp' se interpolan antes de filtrar."""
    x = u.loc[ini:fin].copy(); xi = x.copy(); xi.loc[interp[0]:interp[1]] = np.nan; xi = xi.interpolate()
    return x - pd.Series(filtro_hp(xi.values, lam), index=x.index)


# ------------------------------------------------------------------ Datos trimestrales
q = pd.DataFrame(index=pd.period_range("1996Q1", "2026Q3", freq="Q"))
ipc_q = promedio_q(m.ipc); q["pi"] = 100 * (ipc_q / ipc_q.shift(4) - 1)
q["pi_sv"] = promedio_q(m.inflacion_sv); q["pi_ine"] = promedio_q(m.inflacion_ine)
for v in ("pie_a", "dep", "pet"):
    q[v] = promedio_q(m[v])
q["pie"] = promedio_q(m.eee_11m)
q["u"] = fin_q(m.u_ine_ae); q["u_nsa"] = fin_q(m.u_nsa)
q["brecha"] = brecha_hp(q.u, INI, "2026Q2", ("2020Q2", "2021Q1"), 1600)
q["pi_l1"], q["pi_l4"], q["pie_l1"] = q.pi.shift(1), q.pi.shift(4), q.pie.shift(1)
q["pi_sv_l1"], q["pi_sv_l4"] = q.pi_sv.shift(1), q.pi_sv.shift(4)
q["sorpresa"] = q.pi - q.pie_a
q["sorpresa_sv"] = q.pi_sv - q.pie_a
q["acel"] = q.pi - q.pi_l4
q["acel_sv"] = q.pi_sv - q.pi_sv_l4
q["brecha_l1"] = q.brecha.shift(1)
q["inv_u"] = 1 / q.u
q["estrecho"] = q.brecha * (q.brecha < 0)               # tramo con el mercado laboral más estrecho que su tendencia
per = q.index
q["d_estallido"] = (per == pd.Period("2019Q4")).astype(int)
q["d_pandemia"] = per.year.isin([2020, 2021]).astype(int)

DUM = ["d_estallido", "d_pandemia"]
CTR = ["dep", "pet"]
FORM = {  # nombre: (variable dependiente, regresores distintos de la holgura, holgura)
    "tradicional": ("pi", [], "brecha"),
    "expectativas": ("pi", ["pie_a"], "brecha"),
    "expectativas_r": ("sorpresa", [], "brecha"),
    "aceleracionista": ("acel", [], "brecha"),
    "hibrida": ("pi", ["pie_l1", "pi_l1"], "brecha"),
}


def estimar(datos, y, xs, L, ini=None, fin=None):
    dd = datos.loc[ini or INI:fin or FIN, [y] + xs].dropna()
    paso = np.diff(dd.index.asi8)
    assert (paso == 1).all(), f"Periodos no consecutivos en {y}: Newey-West requiere una serie continua"
    X = np.column_stack([np.ones(len(dd))] + [dd[x].values for x in xs])
    r = newey_west(dd[y].values, X, L); r["xs"] = ["const"] + xs; r["ini"], r["fin"] = str(dd.index[0]), str(dd.index[-1])
    r["dw"] = durbin_watson(r["resid"])
    return r


def comb(r, pesos, valor=0.0):
    """Estimación, error estándar y valor p de la combinación lineal sum(pesos[v] * b[v]) - valor."""
    R = np.zeros(len(r["b"]))
    for v, w in pesos.items():
        R[r["xs"].index(v)] = w
    est = float(R @ r["b"]); se = float(np.sqrt(R @ r["V"] @ R))
    return dict(b=est, se=se, p=float(2 * stats.t.sf(abs(est - valor) / se, r["df"])))


def resumen(r, foco):
    out = dict(n=int(r["n"]), ini=r["ini"], fin=r["fin"], r2a=float(r["r2a"]), dw=float(r["dw"]),
               coefs={v: dict(b=float(b), se=float(s), p=float(p)) for v, b, s, p in zip(r["xs"], r["b"], r["se"], r["p"])})
    out["pendiente"] = out["coefs"][foco]
    return out


res = {}
# ------------------------------------------------------------------ A. Formulaciones (Tabla 1 de la memoria)
for nombre, (y, xs, h) in FORM.items():
    r = estimar(q, y, xs + [h] + CTR + DUM, L_T); res[nombre] = resumen(r, h)
    if nombre == "expectativas":
        res[nombre]["delta_igual_1"] = comb(r, {"pie_a": 1}, 1.0)
    if nombre == "hibrida":
        res[nombre]["suma_igual_1"] = comb(r, {"pie_l1": 1, "pi_l1": 1}, 1.0)
r = estimar(q, "brecha", ["sorpresa", "brecha_l1"] + DUM, L_T)        # curva de oferta de Lucas
res["lucas"] = resumen(r, "sorpresa")

# ------------------------------------------------------------------ B. Estabilidad de la pendiente y no linealidad
def con_post(datos, y, xs, h, ini_post, L, ini=None, fin=None, extra=CTR + DUM):
    dd = datos.copy(); dd["P"] = (dd.index >= pd.Period(ini_post, freq=dd.index.freq)).astype(int)
    dd["HP"] = dd[h] * dd.P
    r = estimar(dd, y, xs + [h, "HP", "P"] + extra, L, ini, fin)
    return dict(antes=comb(r, {h: 1}), despues=comb(r, {h: 1, "HP": 1}), diferencia=comb(r, {"HP": 1}),
                n=int(r["n"]), r2a=float(r["r2a"]), dw=float(r["dw"]), coefs={v: float(b) for v, b in zip(r["xs"], r["b"])})


res["post"] = {}
for nombre in ("expectativas_r", "aceleracionista", "hibrida"):
    y, xs, h = FORM[nombre]
    res["post"][nombre] = con_post(q, y, xs, h, POST, L_T)
res["post_fechas"] = {f: con_post(q, "pi", ["pie_l1", "pi_l1"], "brecha", f, L_T) for f in ("2021Q1", "2021Q3", "2022Q1", "2022Q3")}
res["post_sv"] = con_post(q, "pi_sv", ["pie_l1", "pi_sv_l1"], "brecha", POST, L_T)
# exclusión de un trimestre a la vez del periodo posterior
loo = []
for p in pd.period_range(POST, FIN, freq="Q"):
    dd = q.copy(); dd["P"] = (dd.index >= pd.Period(POST)).astype(int); dd["HP"] = dd.brecha * dd.P
    dd = dd.loc[INI:FIN].drop(p)
    X = np.column_stack([np.ones(len(dd))] + [dd[x].values for x in ["pie_l1", "pi_l1", "brecha", "HP", "P"] + CTR + DUM])
    rr = newey_west(dd.pi.values, X, L_T); rr["xs"] = ["const", "pie_l1", "pi_l1", "brecha", "HP", "P"] + CTR + DUM
    c = comb(rr, {"brecha": 1, "HP": 1}); loo.append(dict(trimestre=str(p), b=c["b"], p=c["p"]))
res["post_excluye_uno"] = loo
r = estimar(q, "pi", ["pie_l1", "pi_l1", "brecha", "estrecho"] + CTR + DUM, L_T)
res["no_lineal_tramo"] = dict(holgado=comb(r, {"brecha": 1}), estrecho=comb(r, {"brecha": 1, "estrecho": 1}),
                              diferencia=comb(r, {"estrecho": 1}), n=int(r["n"]), r2a=float(r["r2a"]), dw=float(r["dw"]))
r = estimar(q, "pi", ["pie_l1", "pi_l1", "inv_u"] + CTR + DUM, L_T)
res["no_lineal_inversa"] = dict(inv_u=comb(r, {"inv_u": 1}), n=int(r["n"]), r2a=float(r["r2a"]), dw=float(r["dw"]))
# inflación sin volátiles
res["sv"] = {}
for nombre, (y, xs) in {"tradicional": ("pi_sv", []), "expectativas_r": ("sorpresa_sv", []),
                        "aceleracionista": ("acel_sv", []), "hibrida": ("pi_sv", ["pie_l1", "pi_sv_l1"])}.items():
    res["sv"][nombre] = resumen(estimar(q, y, xs + ["brecha"] + CTR + DUM, L_T), "brecha")
# sensibilidad a los rezagos de Newey-West (pendiente posterior, formulación híbrida)
res["post_rezagos_nw"] = {}
for L in (2, 4, 6, 8):
    y, xs, h = FORM["hibrida"]; res["post_rezagos_nw"][L] = con_post(q, y, xs, h, POST, L)["despues"]

# ------------------------------------------------------------------ C. Raíz unitaria (muestra principal)
res["raiz"] = {}
for v in ("pi", "pie_a", "sorpresa", "acel", "u", "brecha"):
    x = q.loc[INI:FIN, v].dropna().values
    t0, c0 = adf(x, "c", 4); k0, ck = kpss(x, "c")
    res["raiz"][v] = dict(adf=float(t0), adf_c5=float(c0[5]), adf_c10=float(c0[10]), kpss=float(k0), kpss_c5=ck[5], kpss_c10=ck[10])

# ------------------------------------------------------------------ D. Robustez con datos mensuales, 2010-2025
mm = m.copy()
mm["brecha"] = brecha_hp(mm.u_ine_ae, "2010-03", "2026-09", ("2020-04", "2021-03"), 129600)
mm["pie_l3"], mm["pi_l3"], mm["pi_l12"] = mm.eee_11m.shift(3), mm.pi.shift(3), mm.pi.shift(12)
mm["sorpresa"], mm["acel"] = mm.pi - mm.pie_a, mm.pi - mm.pi_l12
pm = mm.index
mm["d_estallido"] = ((pm >= pd.Period("2019-10")) & (pm <= pd.Period("2019-12"))).astype(int)
mm["d_pandemia"] = pm.year.isin([2020, 2021]).astype(int)
FORM_M = {"tradicional": ("pi", []), "expectativas_r": ("sorpresa", []), "aceleracionista": ("acel", []),
          "hibrida": ("pi", ["pie_l3", "pi_l3"])}
res["mensual"] = {}
for nombre, (y, xs) in FORM_M.items():
    res["mensual"][nombre] = resumen(estimar(mm, y, xs + ["brecha"] + CTR + DUM, L_M, "2010-03", "2025-12"), "brecha")
res["mensual"]["hibrida_post"] = con_post(mm, "pi", ["pie_l3", "pi_l3"], "brecha", "2022-01", L_M, "2010-03", "2025-12")

# ------------------------------------------------------------------ E. Robustez con la muestra 1997-2025
lq = q.copy()
lq["pi_l"] = lq.pi.where(lq.index >= pd.Period("2000Q1"), lq.pi_ine)        # antes de 2000: IPC del INE
lq["d2010"] = (lq.index >= pd.Period("2010Q1")).astype(int)
us = lq.u_nsa.copy(); trim = pd.Series(lq.index.quarter, index=lq.index)
for reg in (0, 1):                                                           # efectos trimestrales por encuesta
    sel = (lq.d2010 == reg) & us.notna()
    ef = us[sel].groupby(trim[sel]).mean() - us[sel].mean(); us[sel] = us[sel] - trim[sel].map(ef)
lq["u_ae"] = us
lq["brecha"] = brecha_hp(lq.u_ae, "1996Q1", "2026Q2", ("2020Q2", "2021Q1"), 1600)
lq["pi_l_l1"], lq["pi_l_l4"] = lq.pi_l.shift(1), lq.pi_l.shift(4)
lq["sorpresa"], lq["acel"] = lq.pi_l - lq.pie_a, lq.pi_l - lq.pi_l_l4
FORM_L = {"tradicional": ("pi_l", [], "1997Q1"), "expectativas_r": ("sorpresa", [], "2002Q3"),
          "aceleracionista": ("acel", [], "1997Q1"), "adaptativa": ("pi_l", ["pi_l_l1"], "1997Q1"),
          "hibrida": ("pi_l", ["pie_l1", "pi_l_l1"], "2002Q1")}
res["larga"], res["larga_pre2010"] = {}, {}
for nombre, (y, xs, ini) in FORM_L.items():
    res["larga"][nombre] = resumen(estimar(lq, y, xs + ["brecha"] + CTR + DUM + ["d2010"], L_T, ini, FIN), "brecha")
    res["larga_pre2010"][nombre] = resumen(estimar(lq, y, xs + ["brecha"] + CTR, L_T, ini, "2009Q4"), "brecha")
res["larga"]["hibrida_post"] = con_post(lq, "pi_l", ["pie_l1", "pi_l_l1"], "brecha", POST, L_T, "2002Q1", FIN,
                                        extra=CTR + DUM + ["d2010"])

json.dump(res, open(RES / "macro.json", "w"), indent=1)
cols = ["pi", "pi_sv", "pie_a", "pie_l1", "u", "brecha", "dep", "pet", "d_estallido", "d_pandemia"]
q.loc[INI:FIN, cols].round(4).rename_axis("trimestre").to_csv(RES / "base_macro_trimestral.csv")

for nombre in list(FORM) + ["lucas"]:
    p = res[nombre]["pendiente"]
    print(f"{nombre:16s} pendiente {p['b']:7.3f} (EE {p['se']:.3f}; p = {p['p']:.3f}) | n = {res[nombre]['n']}, R2 aj. {res[nombre]['r2a']:.3f}")
for nombre, o in res["post"].items():
    print(f"post {nombre:16s} antes {o['antes']['b']:7.3f} (p = {o['antes']['p']:.3f}) después {o['despues']['b']:7.3f} (p = {o['despues']['p']:.3f})")
print("mensual:", {k: round(v["pendiente"]["b"], 3) for k, v in res["mensual"].items() if "pendiente" in v})
print("larga:", {k: round(v["pendiente"]["b"], 3) for k, v in res["larga"].items() if "pendiente" in v})


# ------------------------------------------------------------------ Tablas LaTeX
def num(x, dec=3):
    s = f"{x:.{dec}f}".replace(".", ",")
    return s.replace("-", "$-$") if s.startswith("-") else s


def est(p): return "$^{***}$" if p < .01 else "$^{**}$" if p < .05 else "$^{*}$" if p < .1 else ""
def ee(x): return f"({x:.3f})".replace(".", ",")
def pv(p): return "[$<$0,001]" if p < 0.0005 else f"[{p:.3f}]".replace(".", ",")


COLS = ["tradicional", "expectativas", "expectativas_r", "aceleracionista", "hibrida", "lucas"]
FILAS = [("Brecha de desempleo, $\\tilde U_t$", "brecha"), ("Inflación no anticipada", "sorpresa"),
         ("Expectativa, $\\pi^e_{t\\mid t-4}$", "pie_a"),
         ("Expectativa, $\\pi^e_{t-1}$", "pie_l1"), ("Inflación rezagada, $\\pi_{t-1}$", "pi_l1"),
         ("Brecha rezagada, $\\tilde U_{t-1}$", "brecha_l1"),
         ("Depreciación del peso", "dep"), ("Precio del petróleo", "pet"), ("Constante", "const")]
L = []
for etq, v in FILAS:
    celdas_b, celdas_s = [], []
    for c in COLS:
        cf = res[c]["coefs"].get(v)
        celdas_b.append(num(cf["b"]) + est(cf["p"]) if cf else ""); celdas_s.append(ee(cf["se"]) if cf else "")
    L.append(" & ".join([etq] + celdas_b) + r" \\")
    L.append(" & ".join([""] + celdas_s) + r" \\[2pt]")
f1 = lambda etq, fn: " & ".join([etq] + [fn(c) for c in COLS]) + r" \\"
pie = [f1("$H_0$: $\\delta = 1$ (2); $\\gamma_f + \\gamma_b = 1$ (5)", lambda c: pv(res[c]["delta_igual_1"]["p"]) if c == "expectativas"
          else (pv(res[c]["suma_igual_1"]["p"]) if c == "hibrida" else "")),
       r"\midrule",
       f1("Indicadoras de eventos", lambda c: "Sí"),
       f1("Observaciones", lambda c: str(res[c]["n"])),
       f1("$R^2$ ajustado", lambda c: num(res[c]["r2a"])),
       f1("Durbin--Watson", lambda c: num(res[c]["dw"], 2))]
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{3.5pt}
\caption{Formulaciones de la Curva de Phillips con datos trimestrales, Chile 2010--2025}
\label{tab:macro}
\begin{tabular}{>{\raggedright\arraybackslash}p{4.4cm}cccccc}
\toprule
 & (1) & (2) & (3) & (4) & (5) & (6) \\
 & Tradicional & \multicolumn{2}{c}{Aumentada por expectativas} & Aceleracionista & Híbrida & Lucas \\
Variable dependiente & $\pi_t$ & $\pi_t$ & $\pi_t - \pi^e_{t\mid t-4}$ & $\pi_t - \pi_{t-4}$ & $\pi_t$ & $\tilde U_t$ \\
\midrule
""" + "\n".join(L) + "\n\\midrule\n" + "\n".join(pie) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Mínimos cuadrados ordinarios con datos trimestrales, 2010T1--2025T4. $\pi_t$: variación del IPC promedio
del trimestre respecto del mismo trimestre del año anterior; $\pi^e$: mediana de la inflación esperada a 12 meses en la Encuesta
de Expectativas Económicas del Banco Central de Chile; $\tilde U_t$: tasa de desocupación desestacionalizada menos su tendencia de
Hodrick--Prescott. La depreciación y el precio del petróleo se miden como variaciones logarítmicas en 12 meses (\%). Entre paréntesis,
errores estándar de Newey--West con 4 rezagos; entre corchetes, valores $p$. La inflación no anticipada corresponde a $\pi_t - \pi^e_{t\mid t-4}$. Las indicadoras de eventos corresponden al estallido social (2019T4) y a la pandemia (2020T1--2021T4); sus coeficientes se omiten.
$^{***}p<0{,}01$; $^{**}p<0{,}05$; $^{*}p<0{,}1$. Fuente: elaboración propia a partir del Banco Central de Chile e INE.
\end{minipage}
\end{table}
"""
(TEX / "tabla_macro.tex").write_text(tex, encoding="utf-8")


def celda(c): return f"{num(c['b'])}{est(c['p'])} {pv(c['p'])}"


filas_e = [("Aumentada por expectativas, ecuación \\eqref{eq:m_expr}", res["post"]["expectativas_r"]),
           ("Aceleracionista, ecuación \\eqref{eq:m_acel}", res["post"]["aceleracionista"]),
           ("Híbrida, ecuación \\eqref{eq:m_hib}", res["post"]["hibrida"]),
           ("Híbrida, IPC sin volátiles", res["post_sv"])]
L = [f"{etq} & {celda(o['antes'])} & {celda(o['despues'])} & {pv(o['diferencia']['p'])} \\\\" for etq, o in filas_e]
nl = res["no_lineal_tramo"]
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{5pt}
\caption{Estabilidad y forma de la pendiente de la Curva de Phillips, Chile 2010--2025}
\label{tab:macro_estabilidad}
\begin{tabular}{>{\raggedright\arraybackslash}p{6.0cm}ccc}
\toprule
\multicolumn{4}{l}{\textit{A. Pendiente antes y después de la pandemia}} \\
 & 2010T1--2021T4 & 2022T1--2025T4 & Diferencia \\
\midrule
""" + "\n".join(L) + r"""
\midrule
\multicolumn{4}{l}{\textit{B. No linealidad (formulación híbrida)}} \\
 & Holgura ($\tilde U_t \geq 0$) & Estrechez ($\tilde U_t < 0$) & Diferencia \\
\midrule
""" + f"Pendiente por tramo & {celda(nl['holgado'])} & {celda(nl['estrecho'])} & {pv(nl['diferencia']['p'])} \\\\\n" + \
    f"Coeficiente de $1/U_t$ (forma convexa) & \\multicolumn{{2}}{{c}}{{{celda(res['no_lineal_inversa']['inv_u'])}}} & \\\\" + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Panel A: cada fila corresponde a la ecuación indicada de la Tabla \ref{tab:macro}, a la que se agregan una
indicadora para el periodo 2022T1--2025T4 y su interacción con la brecha de desempleo; la pendiente posterior es la suma de ambos
coeficientes. Panel B: la pendiente se permite distinta cuando la desocupación se ubica bajo su tendencia; en la última fila,
la brecha se reemplaza por el inverso de la tasa de desocupación, como en la forma convexa de \cite{phillips1958} (un coeficiente
positivo indica convexidad). Entre corchetes, valores $p$ con errores de Newey--West (4 rezagos). $^{***}p<0{,}01$; $^{**}p<0{,}05$;
$^{*}p<0{,}1$. Fuente: elaboración propia a partir del Banco Central de Chile e INE.
\end{minipage}
\end{table}
"""
(TEX / "tabla_macro_estabilidad.tex").write_text(tex, encoding="utf-8")

ETQ = [("Tradicional", "tradicional"), ("Aumentada por expectativas (restringida)", "expectativas_r"),
       ("Aceleracionista", "aceleracionista"), ("Adaptativa ($\\pi_{t-1}$)", "adaptativa"), ("Híbrida", "hibrida")]
base_q = {"tradicional": res["tradicional"], "expectativas_r": res["expectativas_r"], "aceleracionista": res["aceleracionista"],
          "hibrida": res["hibrida"]}


def c2(dic, k):
    if k not in dic:
        return "--"
    o = dic[k]; p = o["pendiente"]; return f"{num(p['b'])}{est(p['p'])} {pv(p['p'])}"


L = [f"{etq} & {c2(base_q, k)} & {c2(res['mensual'], k)} & {c2(res['larga'], k)} & {c2(res['larga_pre2010'], k)} \\\\" for etq, k in ETQ]
L.append(r"\midrule")
L.append("Híbrida: pendiente desde 2022 & " + celda(res["post"]["hibrida"]["despues"]) + " & " +
         celda(res["mensual"]["hibrida_post"]["despues"]) + " & " + celda(res["larga"]["hibrida_post"]["despues"]) + r" & -- \\")
ini_txt = lambda dic, k: dic[k]["ini"].replace("Q", "T") if k in dic else ""
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{Pendiente de la Curva de Phillips con distintas frecuencias y muestras}
\label{tab:macro_robustez}
\begin{tabular}{>{\raggedright\arraybackslash}p{4.4cm}cccc}
\toprule
 & Trimestral & Mensual & Trimestral & Trimestral \\
 & 2010--2025 & 2010--2025 & 1997--2025 & 1997--2009 \\
\midrule
""" + "\n".join(L) + r"""
\midrule
""" + f"Observaciones (tradicional) & {res['tradicional']['n']} & {res['mensual']['tradicional']['n']} & {res['larga']['tradicional']['n']} & {res['larga_pre2010']['tradicional']['n']} \\\\" + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Coeficiente de la brecha de desempleo en cada formulación, con los mismos controles de la Tabla
\ref{tab:macro}; entre corchetes, valores $p$ con errores de Newey--West (4 rezagos con datos trimestrales y 12 con datos mensuales).
Con datos mensuales, la inflación y las expectativas rezagadas corresponden a tres meses antes y la tendencia de la desocupación se
estima con $\lambda = 129{.}600$. En la muestra 1997--2025, la desocupación proviene de la antigua ENE hasta el cuarto trimestre de 2009
y de la ENE desde entonces, se desestacionaliza con efectos trimestrales propios de cada encuesta y la regresión incluye una
indicadora de nivel desde 2010; antes de 2000, la inflación corresponde al IPC del INE. Las formulaciones que usan la EEE comienzan
en """ + ini_txt(res["larga"], "hibrida") + " (híbrida) y " + ini_txt(res["larga"], "expectativas_r") + r""" (aumentada por expectativas), por la disponibilidad
de la encuesta. La formulación adaptativa corresponde a la ecuación \eqref{eq:m_hib} sin la expectativa de la EEE, de modo que la inflación del trimestre anterior aproxima las expectativas, y se estima solo en las muestras que comienzan en 1997, en las que la encuesta no cubre todo el periodo. $^{***}p<0{,}01$;
$^{**}p<0{,}05$; $^{*}p<0{,}1$. Fuente: elaboración propia a partir del Banco Central de Chile e INE.
\end{minipage}
\end{table}
"""
(TEX / "tabla_macro_robustez.tex").write_text(tex, encoding="utf-8")

NOM = {"pi": "Inflación, $\\pi_t$", "pie_a": "Expectativa formada un año antes, $\\pi^e_{t\\mid t-4}$",
       "sorpresa": "Inflación no anticipada, $\\pi_t - \\pi^e_{t\\mid t-4}$", "acel": "Variación de la inflación, $\\pi_t - \\pi_{t-4}$",
       "u": "Tasa de desocupación, $U_t$", "brecha": "Brecha de desempleo, $\\tilde U_t$"}
L = []
for v, nm in NOM.items():
    o = res["raiz"][v]
    L.append(f"{nm} & {num(o['adf'], 2)}{'$^{**}$' if o['adf'] < o['adf_c5'] else ''} & {num(o['kpss'], 3)}{'$^{**}$' if o['kpss'] > o['kpss_c5'] else ''} \\\\")
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\caption{Pruebas de raíz unitaria (ADF) y de estacionariedad (KPSS) de las series trimestrales, 2010--2025}
\label{tab:macro_raiz}
\begin{tabular}{lcc}
\toprule
Serie & ADF & KPSS \\
\midrule
""" + "\n".join(L) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.8\linewidth}
\scriptsize \textit{Nota.} Ambas pruebas con constante; ADF con cuatro rezagos y valores críticos de \cite{mackinnon2010}; KPSS con
ventana de Bartlett. ADF: $H_0$ raíz unitaria; KPSS: $H_0$ estacionariedad. $^{**}$ indica rechazo de la hipótesis nula al 5\%.
Fuente: elaboración propia.
\end{minipage}
\end{table}
"""
(TEX / "tabla_macro_raiz.tex").write_text(tex, encoding="utf-8")

# Anexo: base trimestral de la Parte I
AN = q.loc[INI:FIN]
fmt = lambda x, d=2: (f"{x:.{d}f}".replace(".", ",").replace("-", "$-$")) if pd.notna(x) else "--"
filas = [" & ".join([str(i).replace("Q", "T"), fmt(r.pi), fmt(r.pi_sv), fmt(r.pie_a), fmt(r.pie_l1), fmt(r.u, 1), fmt(r.brecha),
                     fmt(r.dep, 1), fmt(r.pet, 1)]) + r" \\" for i, r in AN.iterrows()]
tex = r"""{\footnotesize
\setstretch{1.0}
\setlength{\tabcolsep}{5pt}
\begin{longtable}{lrrrrrrrr}
\caption{Base de datos trimestral de la Parte I, 2010--2025}\label{tab:anexo_macro}\\
\toprule
Trimestre & $\pi_t$ & $\pi^{SV}_t$ & $\pi^e_{t\mid t-4}$ & $\pi^e_{t-1}$ & $U_t$ & $\tilde U_t$ & \textit{dep} & \textit{pet} \\
\midrule
\endfirsthead
\multicolumn{9}{l}{\emph{\tablename\ \thetable\ (continuación)}} \\
\toprule
Trimestre & $\pi_t$ & $\pi^{SV}_t$ & $\pi^e_{t\mid t-4}$ & $\pi^e_{t-1}$ & $U_t$ & $\tilde U_t$ & \textit{dep} & \textit{pet} \\
\midrule
\endhead
\bottomrule
\multicolumn{9}{l}{\scriptsize \textit{Nota.} Todas las variables en porcentaje o en puntos porcentuales. $\pi_t$: inflación; $\pi^{SV}_t$: inflación sin volátiles;} \\
\multicolumn{9}{l}{\scriptsize $\pi^e$: expectativas de la EEE; $U_t$: tasa de desocupación desestacionalizada; $\tilde U_t$: brecha de desempleo;} \\
\multicolumn{9}{l}{\scriptsize \textit{dep}: depreciación del peso; \textit{pet}: variación del precio del petróleo. Fuente: elaboración propia a partir del Banco Central de Chile.} \\
\endlastfoot
""" + "\n".join(filas) + r"""
\end{longtable}
}
"""
(TEX / "anexo_macro.tex").write_text(tex, encoding="utf-8")

# ------------------------------------------------------------------ Figuras
C1, C2, C3 = "#2a78d6", "#eb6834", "#52514e"
INK, INK2, GRID, SOMBRA = "#0b0b0b", "#52514e", "#e4e3df", "#efeeea"
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.dpi": 100, "savefig.dpi": 300, "savefig.bbox": "tight"})
S = q.loc[INI:FIN]; x = S.index.year + (S.index.quarter - 0.5) / 4
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.3, 5.6), sharex=True)
for ax in (a1, a2):
    ax.axvspan(2020, 2022, color=SOMBRA, zorder=0)
a1.text(2021, 0.97, "Pandemia", transform=a1.get_xaxis_transform(), ha="center", va="top", fontsize=8, color=INK2)
a1.plot(x, S.pi, color=C1, lw=2, label="Inflación")
a1.plot(x, S.pie_a, color=C2, lw=1.6, ls=(0, (4, 2)), label="Expectativa formada un año antes (EEE)")
a1.axhline(3, color=INK2, lw=0.8, label="Meta del Banco Central (3%)")
a1.set_ylabel("Variación en 12 meses (%)"); a1.legend(frameon=False, fontsize=8.5, loc="upper left")
a2.plot(x, S.u, color=C1, lw=2, label="Tasa de desocupación (desestacionalizada)")
a2.plot(x, S.u - S.brecha, color=C3, lw=1.4, ls=(0, (4, 2)), label="Tendencia (Hodrick–Prescott)")
a2.set_ylabel("Porcentaje de la fuerza de trabajo"); a2.legend(frameon=False, fontsize=8.5, loc="upper left")
a2.set_xticks(range(2010, 2026, 2)); a2.set_xlim(2010, 2026)
fig.savefig(FIG / "fig7_inflacion_expectativas_desocupacion.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.3, 4.2))
pre = S[(S.index < pd.Period(POST)) & (S.d_pandemia == 0)]; pos = S[S.index >= pd.Period(POST)]
for dd, c, n, mk in [(pre, C1, "2010–2019", "o"), (pos, C2, "2022–2025", "D")]:
    ax.scatter(dd.brecha, dd.sorpresa, s=36, color=c, marker=mk, edgecolor="white", linewidth=1, label=n, zorder=3)
for dd, c in [(pre, C1), (pos, C2)]:
    b1, b0 = np.polyfit(dd.brecha, dd.sorpresa, 1); xs = np.linspace(dd.brecha.min(), dd.brecha.max(), 20)
    ax.plot(xs, b0 + b1 * xs, color=c, lw=1.4)
ax.axhline(0, color=INK2, lw=0.6); ax.axvline(0, color=INK2, lw=0.6)
ax.set_xlabel("Brecha de desempleo (puntos porcentuales)"); ax.set_ylabel("Inflación no anticipada (puntos porcentuales)")
ax.legend(frameon=False, fontsize=8.5, loc="upper right")
fig.savefig(FIG / "fig8_dispersion_brecha_sorpresa.png"); plt.close(fig)
print("Guardado:", RES / "macro.json", "y tablas/figuras de la Parte I")
