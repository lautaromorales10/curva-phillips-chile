"""Ejecuta todo el análisis en orden. Uso:  python ejecutar_todo.py
(Si los microdatos de la ESI ya fueron procesados, puede omitirse el paso 1:  python ejecutar_todo.py --sin-esi)"""
import subprocess, sys, time
from pathlib import Path

PASOS = ["01_procesar_esi.py", "02_armar_base.py", "03_graficos.py", "04_estacionariedad.py",
         "05_estimar_modelos.py", "06_comparar_modelos.py", "07_tablas_latex.py", "08_robustez.py"]
if "--sin-esi" in sys.argv:
    PASOS = PASOS[1:]
aqui = Path(__file__).resolve().parent
for p in PASOS:
    print(f"\n{'=' * 70}\n>>> {p}\n{'=' * 70}", flush=True)
    t = time.time()
    if subprocess.run([sys.executable, str(aqui / p)], cwd=aqui).returncode != 0:
        sys.exit(f"\nError en {p}. Revisa el mensaje de arriba.")
    print(f"({time.time() - t:.0f} s)")
print("\nListo. Revisa la carpeta 'resultados'.")
