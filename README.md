# Curva de Phillips en Chile, 2010–2024

Código y datos de la memoria de título *Validez Empírica de la Curva de Phillips en Chile: Contraste y
Determinantes del Desempleo según Sexo y Composición Educacional* (Ingeniería Civil Industrial,
Universidad Técnica Federico Santa María, 2026). Autor: Lautaro Alen Morales Vásquez.

El análisis se hace íntegramente en Python. El script `ejecutar_todo.py` va desde los datos originales
hasta las tablas en LaTeX y las figuras de la memoria, en un par de minutos, y reproduce todas sus cifras.

## 1. Instalación

Se necesita Python 3. Desde una terminal, en la carpeta del repositorio:

```
pip install -r requirements.txt
```

(En Windows, si `pip` no se reconoce, usar `py -m pip install -r requirements.txt`.)

## 2. Datos

Las cuatro series del Banco Central de Chile (Base de Datos Estadísticos) ya vienen en `datos/`.
Los microdatos de la Encuesta Suplementaria de Ingresos (ESI) del INE no se incluyen en el repositorio:
deben descargarse desde el sitio del INE (https://www.ine.gob.cl/estadisticas-por-tema/mercado-laboral/encuesta-suplementaria-de-ingresos),
en formato Stata (`.dta`), uno por año entre 2010 y 2024, y copiarse en `datos/esi/`. Los archivos pueden
tener cualquier nombre, siempre que incluya el año (por ejemplo, `esi_2021.dta`), y debe haber un solo
archivo por año. Para 2010 a 2012 se utiliza la versión de la base de personas que incluye las becas de
estudio entregadas por el empleador (en la memoria se usaron `esi-2010-personas-cb.dta`,
`esi-2011-personas-cb.dta` y `esi-2012-personas-cb.dta`).

```
curva-phillips-chile/
├── datos/
│   ├── desocupacion_hombres.xlsx            (BDE Banco Central)
│   ├── desocupacion_mujeres.xlsx            (BDE Banco Central)
│   ├── ipc_general.xlsx                     (BDE Banco Central)
│   ├── poblacion_hombres_y_mujeres.xlsx     (BDE Banco Central)
│   └── esi/                                 (microdatos ESI, un .dta por año)
├── config.py, econometria.py
├── 01_procesar_esi.py ... 08_robustez.py
└── ejecutar_todo.py
```

## 3. Ejecución

```
python ejecutar_todo.py
```

Si la ESI ya fue procesada, `python ejecutar_todo.py --sin-esi` se salta ese paso. Cada script también
se puede ejecutar por separado (por ejemplo, `python 05_estimar_modelos.py`). Los resultados quedan en
la carpeta `resultados/`.

## 4. Qué hace cada script

| Script | Qué hace | Genera (en `resultados/`) |
|---|---|---|
| `01_procesar_esi.py` | Ingreso por hora de los asalariados y % de la fuerza de trabajo con educación superior, por sexo y año | `esi_por_anio.csv` |
| `02_armar_base.py` | Une ESI y Banco Central; calcula inflación, salario real (pesos de 2023) y nominal, interacciones, indicadoras y primeras diferencias | `base_apilada_2010_2024.csv` |
| `03_graficos.py` | Figuras descriptivas y figura conceptual de la memoria | `figuras/*.png` |
| `04_estacionariedad.py` | Pruebas ADF y KPSS | `estacionariedad.csv` |
| `05_estimar_modelos.py` | Modelos M1 a M6 (columnas 1 a 6 de la tabla de resultados), análisis de sensibilidad (incluidas indicadoras diferenciadas por sexo) y variantes de la curva salarial (ingreso real o nominal; desempleo en nivel o en diferencias); bootstrap salvaje por año | `resultados_modelos.json` |
| `06_comparar_modelos.py` | Inflación frente a composición educacional (R² ajustado, AIC, BIC, R² parcial, F) y correlación entre las series de hombres y mujeres | `comparacion_modelos.json` |
| `07_tablas_latex.py` | Tablas de la memoria en LaTeX | `tablas/*.tex` |
| `08_robustez.py` | Robustez: rezagos de la inflación, exclusión de un año a la vez y brecha de desempleo (filtro HP, con la tendencia estimada tras interpolar 2020 y 2021) | `robustez.json`, `tablas/tabla_robustez.tex`, `figuras/fig6_exclusion_anual.png` |

`config.py` define las rutas, el número de réplicas del bootstrap (9.999) y la semilla (2026), que hace
que los valores p sean siempre los mismos. `econometria.py` contiene las funciones de estimación
(MCO con errores agrupados, bootstrap salvaje restringido, ADF con valores críticos de MacKinnon (2010) y KPSS).
