# Series mensuales del Banco Central de Chile (BDE) para la parte macro

Archivo: `bcch_mensual_1995_2026.csv` (una fila por mes, ene. 1995 – oct. 2026; celdas vacías = sin dato).
Fuente: Base de Datos Estadísticos del Banco Central de Chile, https://si3.bcentral.cl/siete. Consulta: 2 de noviembre de 2026.
Cada serie se leyó del cuadro de la BDE y se verificó contra la página con una suma de control (n.º de datos, suma y suma ponderada).

| Columna | Código BDE | Cuadro (ruta en la BDE) | Cobertura | Descripción |
|---|---|---|---|---|
| u_ene_antigua | F049.DES.TAS.INE1.10.M (+ .10N.M para ene–feb 2010) | Empleo > Información histórica > Tasas de desocupación nacional y por región > Período 1986–2009 (ED_TDARM) y 2008–2010 (ED_TDNRM) | feb 1986 – feb 2010 | Tasa de desocupación nacional, antigua ENE, trimestre móvil terminado en el mes, % |
| u_ine | F049.DES.TAS.INE9.10.M | Empleo > Tasas de desocupación nacional y por región (ED_TDNRM2) | mar 2010 – sep 2026 | Tasa de desocupación nacional, ENE (proyecciones Censo 2017), trimestre móvil, sin desestacionalizar, % |
| u_ine_ae | F049.DES.TAS.INE.10.M | Empleo > Tasas de desocupación nacional y por género, cifras ajustadas (ED_TASAJ) | mar 2010 – sep 2026 | Ídem, desestacionalizada, % |
| eee_11m / eee_23m | F089.IPC.V12.14.M / F089.IPC.V12.15.M | Expectativas Económicas > Encuesta de Expectativas Económicas > Inflación (EXE_BCCH_01) | 2001 – oct 2026 (mensual desde sep 2001) | Mediana de la inflación esperada en 12 meses, a 11 y 23 meses plazo, % |
| ipc_var_m | F074.IPC.VAR.Z.Z.C.M | Precios > Información histórica > IPC, variación mensual (IPC_VAR_MEN1_HIST_NEW) | ene 1995 – sep 2026 | Variación mensual del IPC, % (una decimal) |
| inflacion_ine | F074.IPC.IND.Z.EP09.C.M | Precios > Información histórica > IPC, variación anual (IPC_VAR_ANUAL_HIST_NEW) | ene 1995 – sep 2026 | Variación en 12 meses del IPC según las bases oficiales del INE, % |
| inflacion_sae | F074.IPCSAE.VIN.Z.EP09.C.M | Ídem | ene 1995 – dic 2023 | IPC sin alimentos ni energía (histórico), variación en 12 meses, % |
| inflacion_bcch | G073.IPC.V12.2023.M | Precios > IPC y analíticos (empalme BCCh), variación anual (PEM_VAR12_IPC_2023) | dic 1999 – sep 2026 | IPC general empalmado por el BCCh (base 2023), variación en 12 meses, % |
| inflacion_sv | G073.IPCSV.V12.2023.M | Ídem | dic 1999 – sep 2026 | IPC sin volátiles, variación en 12 meses, % |
| inflacion_energia | G073.IPCEN.V12.2023.M | Ídem | dic 1999 – sep 2026 | Componente energía (volátiles) del IPC, variación en 12 meses, % |
| dolar_obs | F073.TCO.PRE.Z.D | Dólar Observado y Tipos de Cambio > Dólar observado (DOLAR_OBS_ADO), frecuencia mensual | ene 1995 – oct 2026 | Promedio mensual del dólar observado, pesos por dólar |
| brent / wti | F019.PPB.PRE.41AB.M / F019.PPB.PRE.41B.M | Economía Internacional > Precios de productos básicos, mensual (EI_PROD_BAS) | ene 1995 – sep 2026 | Precio del petróleo, dólares por barril |

## Notas para el uso de las series

1. **Quiebre metodológico de la desocupación en 2010.** La antigua ENE termina en el trimestre dic.–feb. 2010 (8,5%) y la ENE actual parte en ene.–mar. 2010 (9,2%). Para la robustez desde 1997 conviene incluir una variable indicadora desde marzo de 2010 (cambio de nivel) o estimar por separado.
2. **Inflación: empalme BCCh vs. bases oficiales del INE.** Coinciden salvo en los años de cambio de base (2009, 2013, 2018), donde el empalme del Banco Central usa la serie referencial con la canasta nueva; la mayor diferencia es 2009 (hasta 1,3 pp). Para 2010–2025 se recomienda `inflacion_bcch` (es la misma fuente del IPC usado en la memoria); para 1997–1999, `inflacion_ine`.
3. **Expectativas.** La EEE es mensual desde septiembre de 2001; antes hay solo observaciones aisladas. La robustez desde 1997 debe usar expectativas adaptativas (rezagos de la inflación) o la meta de inflación anunciada hasta 2001.
4. El valor de abril de 2026 de `eee_11m` (4,0%) se aparta de los meses vecinos (3,0–3,3%); queda fuera de la muestra 2010–2025, pero conviene revisarlo si se extiende la muestra.
5. Las variaciones de IPC de la BDE vienen con una decimal.
