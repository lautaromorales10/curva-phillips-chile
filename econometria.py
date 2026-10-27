"""Herramientas econométricas (NumPy/SciPy): MCO con errores agrupados, bootstrap salvaje
restringido por grupos (Cameron, Gelbach y Miller, 2008) y pruebas ADF y KPSS."""
import numpy as np
from scipy import stats


# --------------------------------------------------------------------------- MCO
def ols(y, X, cluster=None):
    """MCO. Errores HC1, o agrupados (CR1) si se entrega 'cluster'."""
    y, X = np.asarray(y, float), np.asarray(X, float)
    n, k = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    if cluster is None:
        meat = (X * e[:, None]).T @ (X * e[:, None])
        V = n / (n - k) * XtXi @ meat @ XtXi
        df = n - k
    else:
        g = np.asarray(cluster); G = np.unique(g)
        meat = sum(np.outer(X[g == c].T @ e[g == c], X[g == c].T @ e[g == c]) for c in G)
        V = (len(G) / (len(G) - 1)) * ((n - 1) / (n - k)) * XtXi @ meat @ XtXi
        df = len(G) - 1
    se = np.sqrt(np.diag(V)); t = b / se
    p = 2 * stats.t.sf(np.abs(t), df)
    r2 = 1 - (e @ e) / ((y - y.mean()) @ (y - y.mean()))
    r2a = 1 - (1 - r2) * (n - 1) / (n - k)
    ll = -n / 2 * (np.log(2 * np.pi) + np.log(e @ e / n) + 1)
    return dict(b=b, se=se, t=t, p=p, V=V, df=df, n=n, k=k, r2=r2, r2a=r2a,
                aic=-2 * ll + 2 * k, bic=-2 * ll + k * np.log(n), resid=e)


def wald(res, R):
    """Prueba F de R b = 0 con la matriz de covarianza robusta."""
    R = np.atleast_2d(R); d = R @ res["b"]
    F = d @ np.linalg.inv(R @ res["V"] @ R.T) @ d / R.shape[0]
    return F, stats.f.sf(F, R.shape[0], res["df"])


def durbin_watson(e):
    return float(np.sum(np.diff(e) ** 2) / np.sum(e ** 2))


# --------------------------------------------------------------------------- Bootstrap salvaje
def wild_t(X, Y, g, R, rng, B, unilateral=False):
    """Valor p del bootstrap salvaje restringido para H0: R'b = 0 (una restricción).
    Bilateral por defecto; con unilateral=True contrasta H1: R'b < 0.
    Pesos de Rademacher por grupo; estadístico t con varianza agrupada."""
    full = ols(Y, X, cluster=g)
    t0 = (R @ full["b"]) / np.sqrt(R @ full["V"] @ R)
    # Modelo restringido: proyección de X en el espacio ortogonal a R
    Rn = R / np.linalg.norm(R); P = np.eye(len(R)) - np.outer(Rn, Rn)
    U, s, _ = np.linalg.svd(P); Z = X @ U[:, s > 1e-9]
    fit = Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]; er = Y - fit
    G = np.unique(g); idx = np.searchsorted(G, g)
    W = rng.choice([-1.0, 1.0], size=(B, len(G)))
    # MCO vectorizado para todas las réplicas
    XtXi = np.linalg.inv(X.T @ X); n, k = X.shape
    Yb = fit[None, :] + er[None, :] * W[:, idx]
    bb = Yb @ X @ XtXi
    E = Yb - bb @ X.T
    c = len(G) / (len(G) - 1) * (n - 1) / (n - k)
    xa = X @ (XtXi @ R)
    s_g = np.column_stack([E[:, idx == j] @ xa[idx == j] for j in range(len(G))])
    tb = (bb @ R) / np.sqrt(c * (s_g ** 2).sum(1))
    if unilateral:      # H1: R'b < 0
        return float(t0), float(np.mean(tb <= t0))
    return float(t0), float(np.mean(np.abs(tb) >= abs(t0)))


def _F_cluster(X, Y, g, R):
    r = ols(Y, X, cluster=g); d = R @ r["b"]
    return float(d @ np.linalg.inv(R @ r["V"] @ R.T) @ d / R.shape[0])


def wild_F(X, Y, g, R, rng, B):
    """Bootstrap salvaje restringido para H0: R b = 0 con varias restricciones de exclusión."""
    F0 = _F_cluster(X, Y, g, R)
    keep = [j for j in range(X.shape[1]) if not R[:, j].any()]
    Z = X[:, keep]; fit = Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]; e = Y - fit
    G = np.unique(g); idx = np.searchsorted(G, g)
    Fb = [_F_cluster(X, fit + e * rng.choice([-1., 1.], len(G))[idx], g, R) for _ in range(B)]
    return F0, float(np.mean(np.array(Fb) >= F0))


# --------------------------------------------------------------------------- Raíz unitaria
# Valores críticos de MacKinnon (2010) para muestra finita: b_inf + b1/T + b2/T^2 + b3/T^3
_MK = {"c":  {1: (-3.43035, -6.5393, -16.786, -79.433), 5: (-2.86154, -2.8903, -4.234, -40.040),
              10: (-2.56677, -1.5384, -2.809, 0.0)},
       "ct": {1: (-3.95877, -9.0531, -28.428, -134.155), 5: (-3.41049, -4.3904, -9.036, -45.374),
              10: (-3.12705, -2.5856, -3.925, -22.380)}}


def adf(x, reg="c", lags=0):
    """Dickey-Fuller aumentada. reg='c' (constante) o 'ct' (constante y tendencia)."""
    x = np.asarray(x, float); dx = np.diff(x)
    y = dx[lags:]; T = len(y)
    cols = [x[lags:-1]] + [dx[lags - j:-j] for j in range(1, lags + 1)] + [np.ones(T)]
    if reg == "ct":
        cols.append(np.arange(T, dtype=float))
    X = np.column_stack(cols)
    b = np.linalg.lstsq(X, y, rcond=None)[0]; e = y - X @ b
    s2 = e @ e / (T - X.shape[1]); se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[0, 0])
    crit = {lvl: c[0] + c[1] / T + c[2] / T ** 2 + c[3] / T ** 3 for lvl, c in _MK[reg].items()}
    return b[0] / se, crit


_KPSS = {"c": {10: 0.347, 5: 0.463, 1: 0.739}, "ct": {10: 0.119, 5: 0.146, 1: 0.216}}


def kpss(x, reg="c"):
    """KPSS (H0: estacionariedad), ventana de Bartlett con rezagos 'cortos'."""
    x = np.asarray(x, float); T = len(x)
    X = np.ones((T, 1)) if reg == "c" else np.column_stack([np.ones(T), np.arange(T)])
    e = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]
    L = int(np.floor(4 * (T / 100) ** 0.25))
    s2 = e @ e / T + 2 * sum((1 - l / (L + 1)) * (e[l:] @ e[:-l]) / T for l in range(1, L + 1))
    return (np.cumsum(e) ** 2).sum() / (T ** 2 * s2), _KPSS[reg]
