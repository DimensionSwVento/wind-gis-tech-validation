## Prueba de Concepto PoC - Resultados y Marco de Comparación

### Objetivo general
- **Definir con mayor precisión qué tecnologías podrían adoptarse en el proyecto.**
- **Generar insumos para los resultados del artículo.**

### Alcance de la PoC estado actual
- Flujo de cómputo WSI (Python/rasterio) y generación de insumos:
  - Ráster WSI: `outputs/rasters/wsi.tif`
  - Sitios candidatos: `outputs/vectors/candidate_sites.gpkg`
  - Métricas: `outputs/reports/metrics.json`
- Mapas interactivos Folium/Leaflet:
  - Mapa WSI estándar: `outputs/maps/wsi_map.html`
  - Capas añadidas: múltiples basemaps, Hillshade, WMS OSM/SRTM30, overlay ráster de viento, sitios candidatos cluster y búsqueda, restricciones y red eléctrica si existen, selector regional, medidas, geolocalización, exportar PNG, buffers rápidos 5/10 km y reset mejorado.
- Mapas por límites administrativos (cuando haya estadísticas disponibles):
  - Coroplético por admin: `outputs/maps/wsi_by_<nivel>.html`
  - Flujo OCHA COD-AB: `outputs/maps/wsi_by_<nivel>_ocha.html`

### Cómo reproducir (rápido)
1) Activar entorno y dependencias (una vez):
   - `cd wind-gis-tech-validation`
   - `./venv/Scripts/activate`
   - `pip install -r requirements.txt`
2) Ejecutar WSI + mapa interactivo:
   - `python -m src.interface.cli compute-wsi configs/wsi.yaml -v`
3) Estadísticas y mapa por admin (si aplica):
   - `python -m src.interface.cli compute-wsi-by-admin configs/wsi.yaml -v`
   - `python -m src.interface.cli map-wsi-by-admin configs/wsi.yaml -v`
   - Alternativa OCHA directa: `python -m src.interface.cli map-ocha-wsi configs/wsi.yaml -v`

---

## Marco de comparación de tecnologías

### Tecnologías a considerar centrado en lo instalado y usado
- Python + rasterio + Folium/Leaflet (PoC actual; instalado según `requirements.txt`)

### Escala de puntuación
- 1: muy bajo / 2: bajo / 3: medio / 4: alto / 5: muy alto

### Criterios adicionales
- **Facilidad de instalación y configuración**
- **Curva de aprendizaje e intuitividad de la interfaz**
- **Disponibilidad de documentación y soporte**
- **Facilidad de integración con otras herramientas**
- **Velocidad de procesamiento (GIS)**
- **Precisión espacial (GIS)**
- **Compatibilidad con datos meteorológicos (GIS)**
- **Compatibilidad de formatos GIS (raster/vector/WMS/WFS)**
- **Escalabilidad y automatización (pipelines/CI/CD)**
- **Licenciamiento y costos**

### Evaluación por biblioteca

| Criterio | Numpy | Pandas | GeoPandas | Rasterio | Shapely | Folium | Matplotlib | Typer | PyYAML | Psutil | Requests |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Facilidad de instalación | 5 | 5 | 4 | 4 | 5 | 5 | 5 | 5 | 5 | 5 | 5 |
| Curva de aprendizaje | 3 | 4 | 4 | 3 | 4 | 4 | 3 | 5 | 5 | 4 | 5 |
| Documentación y soporte | 5 | 5 | 4 | 4 | 4 | 4 | 5 | 5 | 4 | 4 | 5 |
| Integración con herramientas | 5 | 5 | 4 | 4 | 5 | 4 | 4 | 5 | 5 | 4 | 5 |
| Velocidad de procesamiento GIS | 4 | 3 | 3 | 4 | 4 | 2 | 2 | 1 | 1 | 1 | 1 |
| Precisión espacial GIS | 2 | 2 | 4 | 5 | 4 | 3 | 2 | 1 | 1 | 1 | 1 |
| Compatibilidad con datos meteorológicos | 2 | 2 | 3 | 3 | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
| Compatibilidad de formatos GIS | 2 | 3 | 4 | 5 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
| Escalabilidad y automatización | 5 | 4 | 3 | 4 | 4 | 3 | 4 | 4 | 4 | 4 | 4 |
| Licenciamiento y costos | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 |

Notas de la tabla por biblioteca:
- Puntuaciones basadas en el uso real en esta PoC y el rol de cada librería.
- "Velocidad/Precisión GIS" favorecen a librerías geoespaciales núcleo (Rasterio, GeoPandas, Shapely).
- Varias librerías no GIS puras mantienen puntaje bajo en filas específicas por no aplicar a ese dominio.

Notas rápidas:
- Python+Folium: instalación vía pip estable en Windows; buen ecosistema (rasterio, geopandas, folium). Procesamiento moderado en CPU; exactitud sólida gestionando CRS y `mask(all_touched=True)`.
- PyQGIS/ArcPy: mayor precisión y tooling GIS, pero instalación/licenciamiento más exigentes (no instalados aquí).
- Deck.gl/Kepler: excelentes para visualización web; no procesan ráster pesado; formatos más limitados.
- GeoServer: ideal para publicar/servir (WMS/WFS/WCS); requiere despliegue Java.

Notas para calificación:
- Indicar versión y SO usados, y cualquier dependencia crítica (GDAL/PROJ/conda/QGIS/ArcGIS).
- Para rendimiento, medir con datasets de referencia y registrar tiempos (I/O, normalización, cálculo, exportación).
- Para precisión, documentar CRS, reproyecciones, resampling y errores conocidos.

---

## Evidencias y mediciones

### Entorno y artefactos
- SO: Windows 10 (win32 10.0.19045)
- Python: venv del proyecto (3.9.x)
- Librerías clave: rasterio, geopandas, shapely, folium, typer, matplotlib
- Artefactos generados:
  - `outputs/rasters/wsi.tif`
  - `outputs/vectors/candidate_sites.gpkg`
  - `outputs/maps/wsi_map.html`
  - `outputs/reports/metrics.json`

### Resultados medidos PoC
- Tiempo total de cómputo WSI: ~1.02 s
- Uso de memoria: ~128.5 MB
- Observaciones:
  - El flujo por admin ahora reproyecta geometrías al CRS del ráster y valida solape (evita errores de “no overlap”).
  - Mapas incluyen: basemaps (OSM/Carto/Stamen), WMS (OSM/SRTM30), Hillshade, ráster de viento opcional, restricciones/red eléctrica (si existen), clustering y búsqueda de sitios, medición, geolocalización, dibujar/editar, buffers 5/10 km, exportar PNG, reset mejorado.

### Run log fragmentos relevantes
- `vento.log` y consola: registrar tiempos, advertencias y errores (e.g., solapes, reproyecciones, geometrías inválidas).

### Métricas de rendimiento
- Tiempo total (s): `outputs/reports/metrics.json`
- Uso de memoria (MB): `outputs/reports/metrics.json`
- Tiempo por etapa: carga → normalización → WSI → exportación → mapas

### Calidad geoespacial
- CRS de entrada/salida, transformaciones aplicadas.
- Validación de solape entre geometrías y ráster.
- Parámetros de `mask()` (`all_touched`) y efectos en bordes.

### Integración y UX
- Capas y controles habilitados en mapas:
  - Basemaps (OSM, Carto, Stamen), WMS (OSM/SRTM30), Hillshade
  - Overlays: WSI points, restricciones, red eléctrica, ráster de viento (opcional)
  - Herramientas: medir, geolocalizar, dibujar/editar, buscar sitios, buffers 5/10 km, exportar PNG, reset
- Facilidad de uso: pasos para cargar datos, filtrar regiones, exportar resultados.

---

## Recomendaciones preliminares (espacio para conclusiones)
- Casos de uso donde cada tecnología destaca.
- Riesgos y mitigaciones (instalación, licencias, soporte, desempeño en datasets grandes).
- Siguientes pasos (pruebas adicionales, datasets reales, integración con pipeline de datos meteo).

---

## Bibliotecas utilizadas en la PoC

| Librería | ¿Usada? | Dónde/para qué |
|---|---|---|
| numpy | Sí | Cálculos y manejo de arreglos en WSI y estadísticas |
| pandas | Sí | Tablas de estadísticas por admin y lectura JSON |
| geopandas | Sí | Lectura/transformación de límites administrativos y AOIs |
| rasterio | Sí | Lectura/escritura de ráster, `mask`/cortes, bounds/CRS |
| shapely | Sí | Geometrías y reparación `buffer(0)` |
| folium | Sí | Mapa interactivo, capas, plugins, controles |
| matplotlib | Sí | Render en memoria para overlay de ráster de viento |
| typer | Sí | CLI `python -m src.interface.cli ...` |
| pyyaml | Sí | Carga de `configs/wsi.yaml` |
| psutil | Sí | Métricas de uso de memoria/CPU |
| requests | Sí | Utilidad de geolocalización fallback en mapa |
| scipy | No | No referenciado en el código actual |
| rioxarray | No | No referenciado en el código actual |
| fiona | No directo | Usado indirectamente por GeoPandas/OGR |
| pyproj | No directo | Usado indirectamente por reproyección en GeoPandas |
| seaborn | No | No referenciado en el código actual |
| jinja2 | No | No referenciado en el código actual |
| reportlab | No | No referenciado en el código actual |
| pytest/flake8/mypy/black | No en runtime | Solo desarrollo/pruebas |

Notas:
- Algunas dependencias se usan de forma transitiva vía GeoPandas/OGR/PROJ.
- La tabla refleja el uso real observado en `src/` en esta PoC.

 


