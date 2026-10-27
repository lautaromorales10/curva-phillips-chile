"""Paso 3. Figuras descriptivas (PNG, 300 dpi). Salida: resultados/figuras/"""
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import RES, FIG

b = pd.read_csv(RES / "base_apilada_2010_2024.csv")
H, M = b[b.sexo == 0].set_index("anio"), b[b.sexo == 1].set_index("anio")

C_H, C_M, C_1 = "#2a78d6", "#eb6834", "#2a78d6"     # slots 1 y 2 de la paleta validada
INK, INK2, GRID, SOMBRA = "#0b0b0b", "#52514e", "#e4e3df", "#efeeea"
plt.rcParams.update({
    "font.family": "DejaVu Serif", "font.size": 10, "axes.edgecolor": INK2,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.dpi": 100, "savefig.dpi": 300, "savefig.bbox": "tight"})
anios = list(range(2010, 2025))

def marco(ax, ylabel, shocks=True):
    if shocks:
        ax.axvspan(2019.5, 2021.5, color=SOMBRA, zorder=0)
        ax.text(2020.5, 0.98, "Pandemia", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8, color=INK2)
        ax.axvline(2019, color=INK2, lw=0.8, ls=(0, (2, 2)), zorder=1)
        ax.text(2019, 0.98, "Estallido ", transform=ax.get_xaxis_transform(),
                ha="right", va="top", fontsize=8, color=INK2)
    ax.set_xticks(anios[::2]); ax.set_xlim(2009.6, 2024.6); ax.set_ylabel(ylabel)

def dos_series(col, ylabel, archivo, loc):
    fig, ax = plt.subplots(figsize=(6.3, 3.4))
    marco(ax, ylabel)
    for d, c, n in [(H, C_H, "Hombres"), (M, C_M, "Mujeres")]:
        ax.plot(d.index, d[col], color=c, lw=2, marker="o", ms=4,
                markeredgecolor="white", markeredgewidth=1, label=n)
        ax.annotate(n, (2024, d[col].iloc[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8.5, **loc)
    fig.savefig(FIG / archivo); plt.close(fig)

dos_series("desempleo", "Tasa de desocupación (%)", "fig1_desempleo_sexo.png", dict(loc="upper left", bbox_to_anchor=(0.22, 1.0)))
dos_series("salario_real_hora", "Ingreso real por hora (pesos de 2023)", "fig3_salario_real_sexo.png", dict(loc="upper left"))
dos_series("educ_sup", "Fuerza de trabajo con ed. superior (%)", "fig4_educacion_sexo.png", dict(loc="lower right"))

# Fig 2: inflación (una serie, sin leyenda)
fig, ax = plt.subplots(figsize=(6.3, 3.0))
marco(ax, "Inflación anual (%)")
ax.plot(H.index, H.inflacion, color=C_1, lw=2, marker="o", ms=4, markeredgecolor="white")
ax.axhline(3, color=INK2, lw=0.8)
ax.text(2009.8, 3.25, "Meta del Banco Central (3%)", fontsize=8, color=INK2)
ax.annotate(f"{H.inflacion[2022]:.1f}%", (2022, H.inflacion[2022]), xytext=(8, 0),
            textcoords="offset points", ha="left", va="center", fontsize=8.5, color=INK)
fig.savefig(FIG / "fig2_inflacion.png"); plt.close(fig)

# Fig 5: dispersión desempleo vs inflación
fig, ax = plt.subplots(figsize=(6.3, 4.0))
for d, c, n in [(H, C_H, "Hombres"), (M, C_M, "Mujeres")]:
    ax.scatter(d.inflacion, d.desempleo, s=42, color=c, edgecolor="white", linewidth=1.5,
               label=n, zorder=3)
    for a in (2020, 2022):
        ax.annotate(str(a), (d.inflacion[a], d.desempleo[a]), xytext=(6, -3),
                    textcoords="offset points", fontsize=8, color=INK2)
ax.set_xlabel("Inflación anual (%)"); ax.set_ylabel("Tasa de desocupación (%)")
ax.legend(frameon=False, fontsize=8.5)
fig.savefig(FIG / "fig5_dispersion_desempleo_inflacion.png"); plt.close(fig)
print("Figuras guardadas en", FIG)

# Figura conceptual (Marco Teórico): curvas de Phillips de corto plazo aumentadas por expectativas
# y curva de largo plazo vertical en la tasa natural. No usa datos.
import numpy as np
plt.rcParams["mathtext.fontset"] = "dejavuserif"
U = np.linspace(4.2, 11, 200); Un = 7.0
fig, ax = plt.subplots(figsize=(6.3, 3.8))
for pe, tono in [(0, "#9cc3ef"), (3, "#5b9be3"), (6, "#1f5fa8")]:
    pi = pe + 9.0 * (Un / U - 1.0)          # curva convexa que pasa por (Un, pe)
    ax.plot(U, pi, color=tono, lw=2)
    ax.annotate(rf"$\pi^e = {pe}\%$", (11, pe + 9.0 * (Un / 11 - 1.0)), xytext=(6, 0),
                textcoords="offset points", va="center", fontsize=9, color=INK, annotation_clip=False)
ax.axvline(Un, color=INK, lw=1.4)
ax.annotate("Curva de largo plazo", (Un, 10.6), xytext=(6, 0), textcoords="offset points",
            va="center", fontsize=9, color=INK)
ax.set_xticks([Un]); ax.set_xticklabels([r"$U^n$"])
ax.set_yticks([0, 3, 6]); ax.set_yticklabels(["0", "3", "6"])
ax.set_xlim(4.2, 11); ax.set_ylim(-4, 11.5)
ax.set_xlabel("Tasa de desempleo"); ax.set_ylabel("Inflación (%)")
ax.grid(False)
fig.savefig(FIG / "fig0_curvas_phillips_teoricas.png"); plt.close(fig)
