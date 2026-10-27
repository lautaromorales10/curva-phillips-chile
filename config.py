"""Rutas del proyecto. Todo se resuelve respecto de la carpeta donde están los scripts,
por lo que basta con respetar la estructura de carpetas descrita en README.md."""
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
DATOS = RAIZ / "datos"              # series del Banco Central (.xlsx)
ESI = DATOS / "esi"                 # microdatos ESI (.dta), un archivo por año
RES = RAIZ / "resultados"           # todo lo que generan los scripts
FIG = RES / "figuras"
TEX = RES / "tablas"
for d in (RES, FIG, TEX):
    d.mkdir(parents=True, exist_ok=True)

ANIOS = range(2010, 2025)
B = 9999                            # réplicas del bootstrap salvaje
SEMILLA = 2026
