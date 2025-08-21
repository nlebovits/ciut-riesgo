# CIUT Riesgo

Un código base para analizar el riesgo de inundación en el Partido de La Plata.

Instrucciones de instalación están disponibles en [SETUP.md](/project-docs/SETUP.md).
El plan de trabajo está disponiblen en [ROADMAP.md](/project-docs/ROADMAP.md).

## Estructura del Proyecto

### Archivos Raíz

El directorio raíz contiene archivos Quarto sin procesar (`.qmd`) que se renderizan para generar el sitio web en Markdown desplegado en GitHub Pages.

### Notebooks de Análisis

Los archivos `.qmd` se basan en análisis realizados en el subdirectorio `notebooks/`, que contiene notebooks de Jupyter. Estos notebooks se limpian y procesan para crear los archivos `.qmd` finales.

### Documentación del Proyecto

La documentación técnica, incluyendo configuración y planificación, se encuentra en el directorio `project-docs/`.

## Recursos para Aprender GIS en Python

Para aprender más sobre el análisis geoespacial en Python, te recomendamos los siguientes recursos (en orden de complejidad):

1. [Introduction to Geospatial Raster and Vector Data with Python](https://carpentries-incubator.github.io/geospatial-python/) - Un excelente punto de partida que cubre los conceptos básicos de datos raster y vectoriales.

2. [PyGIS](https://pygis.io/docs/a_intro.html) - Una guía completa pero accesible para principiantes sobre programación espacial en Python.

3. [Geocomputation with Python](https://py.geocompx.org/) - Un recurso más técnico y detallado que PyGIS, ideal para profundizar en el análisis espacial.

4. [Geographic Data Science with Python](https://geographicdata.science/book/intro.html) - El recurso más avanzado técnicamente, excelente para análisis espacial avanzado y ciencia de datos geoespaciales.

*Nota: Actualmente estoy trabajando en compilar recursos en español.*

## Licencia

Este proyecto está licenciado bajo [GNU General Public License](LICENSE).
