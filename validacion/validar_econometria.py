"""Validación de econometria.py: cada función se compara con una implementación independiente
(fórmulas explícitas y bucles en lugar de cálculo vectorizado), usando los datos del estudio.
Requiere haber ejecutado antes ejecutar_todo.py.   Uso:  python validacion/validar_econometria.py"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))
from config import RES
from econometria import ols, wald, durbin_watson, wild_t, wild_F, adf, kpss

resultados = []
def revisar(condicion, texto):
    resultados.append(bool(condicion))
    print(("OK    " if condicion else "FALLA ") + texto)

b = pd.read_csv(RES / "base_apilada_2010_2024.csv").sort_values(["sexo", "anio"])
for a in (2020, 2021):
    b[f"d{a}"] = (b.anio == a).astype(int)
xs = ["d_inflacion", "d_inflacion_x_sexo", "sexo", "d2019", "d2020", "d2021"]   # columna 3 de la Tabla 7
d = b.dropna(subset=["d_desempleo"] + xs)
X = np.column_stack([np.ones(len(d))] + [d[v].values for v in xs]); Y = d.d_desempleo.values; g = d.anio.values
n, k = X.shape; G = np.unique(g); idx = np.searchsorted(G, g)

# 1-3. MCO, varianza agrupada CR1, R2, AIC y BIC
r = ols(Y, X, cluster=g)
b_ref = np.linalg.lstsq(X, Y, rcond=None)[0]; e = Y - X @ b_ref; XtXi = np.linalg.inv(X.T @ X)
revisar(np.allclose(r["b"], b_ref), "1. Coeficientes MCO iguales a los de numpy.linalg.lstsq")
meat = np.zeros((k, k))
for c in G:
    s = X[g == c].T @ e[g == c]; meat += np.outer(s, s)
V_ref = len(G) / (len(G) - 1) * (n - 1) / (n - k) * XtXi @ meat @ XtXi
revisar(np.allclose(r["V"], V_ref), "2. Varianza agrupada CR1 recalculada año por año")
ll = -n / 2 * (np.log(2 * np.pi) + np.log(e @ e / n) + 1)
revisar(np.isclose(r["aic"], -2 * ll + 2 * k) and np.isclose(r["bic"], -2 * ll + k * np.log(n))
        and np.isclose(r["r2"], 1 - e @ e / ((Y - Y.mean()) ** 2).sum()), "3. R2, AIC y BIC")

# 4-5. Durbin-Watson y Wald
e0 = e[d.sexo.values == 0]
revisar(np.isclose(durbin_watson(e0), sum((e0[t] - e0[t - 1]) ** 2 for t in range(1, len(e0))) / sum(e0 ** 2)),
        "4. Durbin-Watson con la fórmula de sumatorias")
R1 = np.zeros(k); R1[1] = 1
F, pF = wald(r, R1)
revisar(np.isclose(F, (r["b"][1] / r["se"][1]) ** 2) and np.isclose(pF, r["p"][1]), "5. Wald con una restricción igual a t al cuadrado")

# 6-8. Bootstrap salvaje: versión vectorizada frente a un bucle réplica por réplica con los mismos pesos
B = 999
W = np.random.default_rng(7).choice([-1.0, 1.0], size=(B, len(G)))
def p_bucle(R):
    """Modelo restringido por la fórmula de mínimos cuadrados restringidos (sin usar la SVD del código)."""
    br = b_ref - XtXi @ R * (R @ b_ref) / (R @ XtXi @ R); fit = X @ br; er = Y - fit
    t0 = (R @ b_ref) / np.sqrt(R @ r["V"] @ R); tb = []
    for w in W:
        rb = ols(fit + er * w[idx], X, cluster=g); tb.append((R @ rb["b"]) / np.sqrt(R @ rb["V"] @ R))
    return np.mean(np.abs(np.array(tb)) >= abs(t0))
t0, p_vec = wild_t(X, Y, g, R1, np.random.default_rng(7), B)
revisar(np.isclose(p_vec, p_bucle(R1)), f"6. Bootstrap salvaje: versión vectorizada = bucle (p = {p_vec:.4f})")
R2 = np.zeros(k); R2[1] = R2[2] = 1
revisar(np.isclose(wild_t(X, Y, g, R2, np.random.default_rng(7), B)[1], p_bucle(R2)),
        "7. Bootstrap salvaje para la pendiente de mujeres (combinación lineal)")
F0, pFb = wild_F(X, Y, g, R1[None, :], np.random.default_rng(7), B)
revisar(np.isclose(F0, t0 ** 2) and np.isclose(pFb, p_vec), "8. Bootstrap F con una restricción = bootstrap t bilateral")

# 9-10. ADF y KPSS
u = b[b.sexo == 0].desempleo.values
t_adf, cr = adf(u, "c")
dy = np.diff(u); Z = np.column_stack([u[:-1], np.ones(len(dy))]); bz = np.linalg.solve(Z.T @ Z, Z.T @ dy)
ez = dy - Z @ bz; se = np.sqrt(ez @ ez / (len(dy) - 2) * np.linalg.inv(Z.T @ Z)[0, 0]); T = len(dy)
revisar(np.isclose(t_adf, bz[0] / se) and np.isclose(cr[5], -2.86154 - 2.8903 / T - 4.234 / T ** 2 - 40.040 / T ** 3),
        "9. ADF: estadístico t de Dickey-Fuller y valor crítico de MacKinnon (2010, Tabla 2)")
eta, kc = kpss(u, "c"); ee = u - u.mean(); Tn = len(u); L = int(4 * (Tn / 100) ** 0.25)
s2 = sum(ee ** 2) / Tn + 2 / Tn * sum((1 - l / (L + 1)) * sum(ee[t] * ee[t - l] for t in range(l, Tn)) for l in range(1, L + 1))
revisar(np.isclose(eta, sum(np.cumsum(ee) ** 2) / (Tn ** 2 * s2)) and kc[5] == 0.463,
        "10. KPSS con la fórmula de Kwiatkowski et al. (1992)")

# 11. Filtro HP (se toma la función desde 08_robustez.py)
src = (RAIZ / "08_robustez.py").read_text(encoding="utf-8"); ns = {"np": np}
exec(src[src.index("def filtro_hp"):src.index("def brecha")], ns); hp = ns["filtro_hp"]
tau = hp(u, 6.25); K = np.diff(np.eye(Tn), 2, axis=0)
lineal = np.polyval(np.polyfit(np.arange(Tn), u, 1), np.arange(Tn))
revisar(np.allclose((u - tau) - 6.25 * K.T @ (K @ tau), 0) and np.allclose(hp(u, 1e-9), u)
        and np.allclose(hp(u, 1e8), lineal, atol=1e-5), "11. Filtro HP: condición de primer orden y casos límite")

print(f"\n{sum(resultados)} de {len(resultados)} comprobaciones correctas.")
sys.exit(0 if all(resultados) else 1)
