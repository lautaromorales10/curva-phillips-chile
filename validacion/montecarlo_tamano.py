"""Simulación de Monte Carlo del tamaño de las pruebas con el diseño del estudio (modelo de la columna 3
de la Tabla 7: 14 años x 2 sexos, con las mismas variables e indicadoras). Se generan datos en que la
inflación no afecta al desempleo y se cuenta cuántas veces cada prueba rechaza la hipótesis nula al 5%.
Requiere haber ejecutado antes ejecutar_todo.py.   Uso:  python validacion/montecarlo_tamano.py"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))
from config import RES, SEMILLA
from econometria import ols, wild_t

REPLICAS, B = 2000, 399
b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
for a in (2020, 2021):
    b[f"d{a}"] = (b.anio == a).astype(int)
xs = ["d_inflacion", "d_inflacion_x_sexo", "sexo", "d2019", "d2020", "d2021"]
d = b.dropna(subset=["d_desempleo"] + xs)
X = np.column_stack([np.ones(len(d))] + [d[v].values for v in xs]); g = d.anio.values
G = np.unique(g); idx = np.searchsorted(G, g)
R = np.zeros(X.shape[1]); R[1] = 1                       # H0: pendiente de los hombres igual a cero
extremo = np.abs(X[:, 1]) > 3                            # años con variaciones de la inflación mayores a 3 puntos

print(f"Rechazo de una hipótesis nula verdadera al 5% ({REPLICAS} simulaciones; error de simulación ~±1 punto)")
for heterocedastico in (False, True):
    for rho in (0.0, 0.5):                               # correlación de los errores de hombres y mujeres en un año
        rng = np.random.default_rng(SEMILLA)
        rech_conv = rech_boot = 0
        for _ in range(REPLICAS):
            u = np.sqrt(rho) * rng.standard_normal(len(G))[idx] + np.sqrt(1 - rho) * rng.standard_normal(len(g))
            if heterocedastico:
                u *= np.where(extremo, 3.0, 1.0)         # perturbaciones más volátiles en los años de inflación extrema
            r = ols(u, X, cluster=g)
            rech_conv += r["p"][1] < 0.05
            rech_boot += wild_t(X, u, g, R, rng, B)[1] < 0.05
        etq = "errores más volátiles en años extremos" if heterocedastico else "errores de igual varianza"
        print(f"  {etq:<40s} | correlación {rho:.1f} | agrupados convencionales: {100 * rech_conv / REPLICAS:5.1f}%"
              f" | bootstrap salvaje: {100 * rech_boot / REPLICAS:5.1f}%")
