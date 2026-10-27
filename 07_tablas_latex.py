"""Paso 7. Tablas en LaTeX (decimales con coma) para el informe. Salida: resultados/tablas/*.tex"""
import json, pandas as pd
from config import RES, TEX

R = json.load(open(RES / "resultados_modelos.json"))

def num(x, d=3):
    s = f"{x:.{d}f}".replace(".", ",")
    return s.replace("-", "$-$") if s.startswith("-") else s

def est(p):
    if p is None: return ""
    return "$^{***}$" if p < .01 else "$^{**}$" if p < .05 else "$^{*}$" if p < .1 else ""

# ------------------------------------------------------------------ Tabla de regresiones
cols = ["M1", "M2", "M3", "M4", "M5", "M6"]
filas = [("Inflación", ["inflacion", "d_inflacion"]),
         ("Inflación $\\times$ Mujer", ["inflacion_x_sexo", "d_inflacion_x_sexo"]),
         ("Educación superior", ["educ_sup", "d_educ_sup"]),
         ("Educación superior $\\times$ Mujer", ["educ_sup_x_sexo", "d_educ_sup_x_sexo"]),
         ("Desempleo", ["desempleo"]),
         ("Desempleo $\\times$ Mujer", ["desempleo_x_sexo"]),
         ("Mujer", ["sexo"]),
         ("Constante", ["const"])]
L = []
for etq, vs in filas:
    a, s = [etq], [""]
    for m in cols:
        c = R[m]["coefs"]; v = next((x for x in vs if x in c), None)
        if v is None: a.append(""); s.append(""); continue
        a.append(num(c[v]["b"]) + est(c[v]["p_wb"])); s.append(f"({c[v]['se']:.3f})".replace(".", ","))
    L += [" & ".join(a) + r" \\", " & ".join(s) + r" \\[2pt]"]

def fila(etq, f): return " & ".join([etq] + [f(m) for m in cols]) + r" \\"
pie = [
    fila("Pendiente para mujeres", lambda m: num(R[m]["mujeres"]["b"]) + est(R[m]["mujeres"]["p_wb"])),
    fila("\\quad valor $p$ (bootstrap)", lambda m: f"[{R[m]['mujeres']['p_wb']:.3f}]".replace(".", ",")),
    r"\midrule",
    fila("Indicadora estallido (2019)", lambda m: "Sí"),
    fila("Indicadora(s) pandemia", lambda m: "Periodo" if m in ("M1", "M2") else "Por año"),
    fila("Indicadora 2022", lambda m: "Sí" if m == "M5" else "No"),
    r"\midrule",
    fila("Observaciones", lambda m: str(R[m]["n"])),
    fila("$R^2$ ajustado", lambda m: num(R[m]["r2a"])),
    fila("Durbin--Watson (H / M)", lambda m: f"{R[m]['dw'][0]:.2f} / {R[m]['dw'][1]:.2f}".replace(".", ",")),
]
tex = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{Estimación de la Curva de Phillips y de los determinantes del desempleo, Chile 2010--2024}
\label{tab:resultados}
\begin{tabular}{lcccccc}
\toprule
 & \multicolumn{2}{c}{Niveles} & \multicolumn{3}{c}{Primeras diferencias} & Salarial \\
\cmidrule(lr){2-3}\cmidrule(lr){4-6}\cmidrule(lr){7-7}
 & (1) & (2) & (3) & (4) & (5) & (6) \\
Variable dependiente & $U$ & $U$ & $\Delta U$ & $\Delta U$ & $\Delta U$ & $\Delta \ln w$ \\
\midrule
""" + "\n".join(L) + "\n\\midrule\n" + "\n".join(pie) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Mínimos cuadrados ordinarios sobre la base apilada por sexo. Mujer corresponde a la variable
dicotómica de sexo, $S_g$ (1 = mujeres; 0 = hombres). Entre paréntesis,
errores estándar agrupados por año. La significancia se evalúa mediante bootstrap salvaje restringido por año con
pesos de Rademacher y 9.999 réplicas \citep{cameron2008}, adecuado para un número reducido de grupos. En las columnas
(3) a (5) las variables explicativas continuas están en primeras diferencias; en la columna (6) la variable
dependiente es la variación porcentual del ingreso real por hora. Los coeficientes de Inflación, Educación superior y
Desempleo corresponden a la pendiente para los hombres; la fila ``Pendiente para mujeres'' reporta la suma con su
interacción. Los coeficientes de las variables indicadoras se omiten por identificarse con un solo año.
$^{***}p<0{,}01$; $^{**}p<0{,}05$; $^{*}p<0{,}1$. Fuente: elaboración propia a partir de INE y Banco Central de Chile.
\end{minipage}
\end{table}
"""
open(TEX / "tabla_resultados.tex", "w", encoding="utf-8").write(tex)

# ------------------------------------------------------------------ Tabla: curva salarial alternativa
S = R["salarial"]
def fs(etq, r):
    return " & ".join([etq, num(r["h"]), "[" + num(r["p_h"]) + "]", num(r["m"]), "[" + num(r["p_m"]) + "]",
                       num(r["r2a"]), str(r["n"])]) + r" \\"
filas_s = [r"\multicolumn{7}{l}{\textit{Desempleo en nivel, $U_{g,t}$}} \\",
           fs(r"\quad (1) Ingreso real por hora (Tabla \ref{tab:resultados}, col. 6)", S["real_nivel"]),
           fs(r"\quad (2) Ingreso nominal por hora", S["nominal_nivel"]),
           r"\midrule",
           r"\multicolumn{7}{l}{\textit{Desempleo en primeras diferencias, $\Delta U_{g,t}$}} \\",
           fs(r"\quad (3) Ingreso real por hora", S["real_dif"]),
           fs(r"\quad (4) Ingreso nominal por hora", S["nominal_dif"])]
tex_s = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{5pt}
\caption{Curva de Phillips salarial: especificaciones alternativas}
\label{tab:salarial}
\begin{tabular}{lcccccc}
\toprule
 & \multicolumn{2}{c}{Hombres} & \multicolumn{2}{c}{Mujeres} & & \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Variable dependiente & Pendiente & Valor $p$ & Pendiente & Valor $p$ & $R^2$ ajust. & $n$ \\
\midrule
""" + "\n".join(filas_s) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} Mínimos cuadrados ordinarios sobre la base apilada por sexo, con errores agrupados por año; valores $p$
(entre corchetes) obtenidos mediante bootstrap salvaje restringido por año con 9.999 réplicas. La variable dependiente es la
variación porcentual del ingreso por hora (diferencia del logaritmo multiplicada por 100). La pendiente corresponde al coeficiente
de la tasa de desocupación, o de su variación, para cada sexo. Todas las especificaciones incluyen la variable de sexo y las
indicadoras de 2019, 2020 y 2021. Elaboración propia a partir de los microdatos de la ESI (INE) y del Banco Central de Chile.
\end{minipage}
\end{table}
"""
open(TEX / "tabla_salarial.tex", "w", encoding="utf-8").write(tex_s)

# ------------------------------------------------------------------ Tabla de estacionariedad
e = pd.read_csv(RES / "estacionariedad.csv")
nom = {"desempleo": "Tasa de desempleo", "inflacion": "Inflación", "educ_sup": "Educación superior",
       "salario_real_hora": "Log. ingreso real por hora"}
sx = {"H": "Hombres", "M": "Mujeres", "común": "--"}
Ls = []
for _, r in e.iterrows():
    Ls.append(" & ".join([
        nom[r.serie], sx[r.sexo], "Nivel" if r.forma == "nivel" else "Diferencia",
        "C + T" if r.det == "ct" else "C",
        num(r.adf, 2) + ("$^{**}$" if r.adf_rechaza else ""), num(r.crit5, 2),
        f"{r.kpss:.3f}".replace(".", ",") + ("$^{**}$" if r.kpss_rechaza else ""), f"{r.kpss_crit5:.3f}".replace(".", ",")]) + r" \\")
tex2 = r"""\begin{table}[htbp]
\centering
\footnotesize
\caption{Pruebas de raíz unitaria (ADF) y de estacionariedad (KPSS), 2010--2024}
\label{tab:estacionariedad}
\begin{tabular}{lllcrrrr}
\toprule
Serie & Sexo & Forma & Det. & ADF & Crít. 5\% & KPSS & Crít. 5\% \\
\midrule
""" + "\n".join(Ls) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.92\linewidth}
\scriptsize \textit{Nota.} ADF: prueba de Dickey--Fuller aumentada sin rezagos, cuya hipótesis nula es la presencia
de raíz unitaria; valores críticos para muestra finita según \citet{mackinnon2010}. KPSS: hipótesis nula de
estacionariedad \citep{kpss1992}. Det.: componentes determinísticos (C = constante; T = tendencia lineal).
$^{**}$ indica rechazo de la hipótesis nula al 5\%. Con quince observaciones, ambas pruebas tienen baja potencia.
\end{minipage}
\end{table}
"""
open(TEX / "tabla_estacionariedad.tex", "w", encoding="utf-8").write(tex2)

# ------------------------------------------------------------------ Tabla de comparación (objetivo 5)
C = json.load(open(RES / "comparacion_modelos.json"))

def fc(x, d=3):
    if x is None or x != x: return "--"
    s = f"{x:.{d}f}".replace(".", ",")
    return s.replace("-", "$-$") if s.startswith("-") else s

def panel(clave, titulo):
    L = [r"\multicolumn{7}{l}{\textit{" + titulo + r" (n = " + str(C[clave]["n"]) + r")}} \\"]
    for r in C[clave]["tabla"]:
        F = None if r["F"] is None or r["F"] != r["F"] else r["F"]
        p = None if r["p"] is None or r["p"] != r["p"] else r["p"]
        L.append(" & ".join([r"\quad " + r["modelo"], fc(r["r2a"]), fc(r["aic"], 1), fc(r["bic"], 1),
                             "--" if r["r2_parcial"] == 0 else fc(r["r2_parcial"]),
                             fc(F, 2), fc(p)]) + r" \\")
    for nom, (F, p) in C[clave]["extra"].items():
        etq = "Aporte de la inflación, dada la educación" if nom.startswith("inflación") \
              else "Aporte de la educación, dada la inflación"
        L.append(" & ".join([r"\quad " + etq, "", "", "", "", fc(F, 2), fc(p)]) + r" \\")
    return L

tex3 = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{5pt}
\caption{Comparación de la capacidad explicativa de la inflación y de la composición educacional sobre la tasa de desocupación}
\label{tab:comparacion}
\begin{tabular}{lcccccc}
\toprule
Modelo & $R^2$ ajust. & AIC & BIC & $R^2$ parcial & $F$ & Valor $p$ \\
\midrule
""" + "\n".join(panel("dif", "Panel A. Primeras diferencias")) + "\n\\midrule\n" + \
"\n".join(panel("niv", "Panel B. Niveles")) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.95\linewidth}
\scriptsize \textit{Nota.} El modelo base incluye la variable de sexo y las indicadoras de eventos (en diferencias, 2019, 2020 y 2021;
en niveles, 2019 y el periodo 2020--2021). Cada bloque explicativo incluye la variable y su interacción con sexo. El $R^2$ parcial
corresponde a la proporción de la variación no explicada por el modelo base que explica el bloque agregado. AIC y BIC: criterios de
información de Akaike y de Schwarz; valores menores indican un mejor modelo. El estadístico $F$ contrasta la nulidad conjunta del
bloque; su valor $p$ se obtiene mediante bootstrap salvaje restringido por año con 9.999 réplicas \citep{cameron2008}.
Fuente: elaboración propia a partir de INE y Banco Central de Chile.
\end{minipage}
\end{table}
"""
open(TEX / "tabla_comparacion.tex", "w", encoding="utf-8").write(tex3)
print("Tablas guardadas en", TEX)

# ------------------------------------------------------------------ Estadísticas descriptivas
import numpy as np
b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
VARS_D = [("Tasa de desocupación (\\%)", "desempleo", 2), ("Inflación (\\%)", "inflacion", 2),
          ("Fuerza de trabajo con educación superior (\\%)", "educ_sup", 2),
          ("Ingreso real por hora (pesos de 2023)", "salario_real_hora", 0),
          ("Variación de la tasa de desocupación (pp.)", "d_desempleo", 2),
          ("Variación de la inflación (pp.)", "d_inflacion", 2),
          ("Variación de la educación superior (pp.)", "d_educ_sup", 2),
          ("Variación del ingreso real por hora (\\%)", "d_salario", 2)]


def nd(x, d):
    if d == 0:
        return f"{x:,.0f}".replace(",", ".")
    return num(x, d)


Ld = []
for i, (etq, v, d) in enumerate(VARS_D):
    if i == 4:
        Ld.append(r"\midrule")
        Ld.append(r"\multicolumn{9}{l}{\textit{Primeras diferencias}} \\")
    fila = [etq]
    for s in (0, 1):
        x = b[b.sexo == s][v].dropna()
        fila += [nd(x.mean(), d), nd(x.std(ddof=1), d), nd(x.min(), d), nd(x.max(), d)]
    Ld.append(" & ".join(fila) + r" \\")
tex4 = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{Estadísticas descriptivas de las variables por sexo, 2010--2024}
\label{tab:descriptivas}
\begin{tabular}{>{\raggedright\arraybackslash}p{4.6cm}rrrrrrrr}
\toprule
 & \multicolumn{4}{c}{Hombres} & \multicolumn{4}{c}{Mujeres} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-9}
Variable & Media & D.E. & Mín. & Máx. & Media & D.E. & Mín. & Máx. \\
\midrule
\multicolumn{9}{l}{\textit{Niveles}} \\
""" + "\n".join(Ld) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} D.E.: desviación estándar; pp.: puntos porcentuales. Niveles: 15 observaciones por sexo; primeras
diferencias: 14. La inflación es común a ambos sexos. La variación del ingreso real corresponde a la diferencia de su logaritmo
multiplicada por 100. Fuente: elaboración propia a partir de INE y Banco Central de Chile.
\end{minipage}
\end{table}
"""
open(TEX / "tabla_descriptivas.tex", "w", encoding="utf-8").write(tex4)

# ------------------------------------------------------------------ Anexo: resultados de la ESI por año
esi = pd.read_csv(RES / "esi_por_anio.csv")
La = []
for a in sorted(esi.anio.unique()):
    h, m = esi[(esi.anio == a) & (esi.sexo == 0)].iloc[0], esi[(esi.anio == a) & (esi.sexo == 1)].iloc[0]
    jor = r"\texttt{" + h.variables_jornada.replace("_", r"\_").replace(" x ", r"} $\times$ \texttt{") + "}"
    La.append(" & ".join([str(a), f"{h.n_asalariados:,}".replace(",", "."), f"{m.n_asalariados:,}".replace(",", "."),
                          nd(h.salario_hora_nominal, 0), nd(m.salario_hora_nominal, 0),
                          num(h.educ_sup, 1), num(m.educ_sup, 1), jor]) + r" \\")
tex5 = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{7pt}
\caption{Resultados del procesamiento de los microdatos de la ESI por año}
\label{tab:anexo_esi}
\begin{tabular}{lrrrrrrl}
\toprule
 & \multicolumn{2}{c}{Asalariados} & \multicolumn{2}{c}{Ingreso por hora} & \multicolumn{2}{c}{Ed. superior} & \\
 & \multicolumn{2}{c}{(n)} & \multicolumn{2}{c}{(\$ nominales)} & \multicolumn{2}{c}{(\%)} & \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}
Año & H & M & H & M & H & M & Variables de jornada \\
\midrule
""" + "\n".join(La) + r"""
\bottomrule
\end{tabular}

\vspace{4pt}
\begin{minipage}{0.97\linewidth}
\scriptsize \textit{Nota.} H: hombres; M: mujeres. Asalariados: número de observaciones muestrales (sin expandir) que cumplen los
filtros. El ingreso por hora y el porcentaje de la fuerza de trabajo con educación superior están ponderados por el factor de
expansión. Fuente: elaboración propia a partir de los microdatos de la ESI (INE).
\end{minipage}
\end{table}
"""
open(TEX / "anexo_esi.tex", "w", encoding="utf-8").write(tex5)

# ------------------------------------------------------------------ Anexo: base de datos completa
Lb = []
for _, r in b.sort_values(["anio", "sexo"]).iterrows():
    Lb.append(" & ".join([str(int(r.anio)), "Mujer" if r.sexo == 1 else "Hombre", num(r.desempleo, 2), num(r.inflacion, 2),
                          num(r.educ_sup, 2), nd(r.salario_real_hora, 0),
                          "--" if pd.isna(r.d_desempleo) else num(r.d_desempleo, 2),
                          "--" if pd.isna(r.d_salario) else num(r.d_salario, 2)]) + r" \\")
tex6 = r"""{\footnotesize
\setlength{\tabcolsep}{6pt}
\begin{longtable}{llrrrrrr}
\caption{Base de datos apilada por sexo, 2010--2024}\label{tab:anexo_base}\\
\toprule
Año & Sexo & $U$ (\%) & $\pi$ (\%) & $E$ (\%) & $w$ (\$ de 2023) & $\Delta U$ (pp.) & $\Delta \ln w$ (\%) \\
\midrule
\endfirsthead
\multicolumn{8}{l}{\emph{\tablename\ \thetable\ (continuación)}} \\
\toprule
Año & Sexo & $U$ (\%) & $\pi$ (\%) & $E$ (\%) & $w$ (\$ de 2023) & $\Delta U$ (pp.) & $\Delta \ln w$ (\%) \\
\midrule
\endhead
\bottomrule
\multicolumn{8}{l}{\scriptsize \textit{Nota.} $U$: tasa de desocupación; $\pi$: inflación; $E$: fuerza de trabajo con educación superior;} \\
\multicolumn{8}{l}{\scriptsize $w$: ingreso real por hora. Fuente: elaboración propia a partir de INE y Banco Central de Chile.} \\
\endlastfoot
""" + "\n".join(Lb) + r"""
\end{longtable}
}
"""
open(TEX / "anexo_base.tex", "w", encoding="utf-8").write(tex6)
print("Tablas adicionales guardadas.")
