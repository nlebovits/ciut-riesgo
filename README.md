# CIUT Riesgo

Un código base para analizar el riesgo de inundación en el Partido de La Plata.

## Software Requerido

Instala el software en el siguiente orden:

1. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) y [VSCode](https://code.visualstudio.com/download) - Para clonar y editar el repositorio
2. [Python](https://www.python.org/downloads/) - Necesario para las herramientas de Python
3. [`uv`](https://docs.astral.sh/uv/) y [`ruff`](https://docs.astral.sh/ruff/) - Herramientas de Python
4. [Google Cloud CLI](https://cloud.google.com/sdk/docs/install) - Para operaciones en la nube

La instalación de todo el software debería tomar no más de 15-20 minutos, asumiendo que no tienes ninguno instalado. (Python y Git ya vienen incluidos en sistemas macOS o Linux.)

### Para Windows

Para usuarios de Windows, la instalación puede ser un poco más compleja. Específicamente:

1. Descarga e instala [Git for Windows](https://github.com/git-for-windows/git/releases/latest) desde el enlace oficial.

2. La instalación de Python puede ser complicada porque no se agrega automáticamente al PATH del sistema. Durante la instalación de Python, asegúrate de marcar la casilla "Add Python to PATH". Si ya instalaste Python sin esta opción, puedes agregar Python al PATH usando PowerShell:

```powershell
# Obtener la última versión de Python instalada
$pythonVersion = (Get-ChildItem "$env:USERPROFILE\AppData\Local\Programs\Python" -Directory | Sort-Object Name -Descending | Select-Object -First 1).Name
$pythonPath = "$env:USERPROFILE\AppData\Local\Programs\Python\$pythonVersion"
$pythonScriptsPath = "$pythonPath\Scripts"

# Agregar Python al PATH del sistema
[System.Environment]::SetEnvironmentVariable('PATH', $env:PATH + ";$pythonPath;$pythonScriptsPath", "User")

# Verificar la instalación
python --version
pip --version
```

## Instalación

1. Clona este repositorio:

```bash
git clone https://github.com/nlebovits/ciut-tablero.git
cd ciut-tablero
```

2. Instala las dependencias usando uv:

```bash
uv sync
```

Para actualizar el archivo de bloqueo:

```bash
uv lock
```

3. Configura las credenciales:

Coloca tu clave de servicio de Google Cloud Storage en el directorio `credentials/`.

4. Configura las variables de entorno:

```bash
cp .env.example .env
```

Edita el archivo `.env` con tus credenciales. Consulta `.env.example` para ver la estructura requerida.

## Hoja de ruta

1.  **Automatizar el análisis manual de riesgo a nivel de cuenca utilizando Python para facilitar su actualización.**
    1.  Crear un **Producto Mínimo Viable (MVP)** combinando datos de archivos personales.
    2.  Actualizar el análisis para que sea 100% replicable y automatizado, utilizando **APIs y Geoservicios (incluido REDATAM)**.
2.  **Incorporar un análisis de barrios populares.**
    1.  Calcular el riesgo a nivel de barrio, empleando datos de **RENABAP** y datos censales.
    2.  Integrar datos de edificios para reducir la escala del análisis de exposición.
    3.  (Opcional) Calcular la expansión de barrios populares mediante datos de superficie impermeable y análisis de tendencias.

## Recursos para Aprender GIS en Python

Para aprender más sobre el análisis geoespacial en Python, te recomendamos los siguientes recursos (en orden de complejidad):

1. [Introduction to Geospatial Raster and Vector Data with Python](https://carpentries-incubator.github.io/geospatial-python/) - Un excelente punto de partida que cubre los conceptos básicos de datos raster y vectoriales.

2. [PyGIS](https://pygis.io/docs/a_intro.html) - Una guía completa pero accesible para principiantes sobre programación espacial en Python.

3. [Geocomputation with Python](https://py.geocompx.org/) - Un recurso más técnico y detallado que PyGIS, ideal para profundizar en el análisis espacial.

4. [Geographic Data Science with Python](https://geographicdata.science/book/intro.html) - El recurso más avanzado técnicamente, excelente para análisis espacial avanzado y ciencia de datos geoespaciales.

*Nota: Actualmente estoy trabajando en compilar recursos en español.*

## Licencia

Este proyecto está licenciado bajo [GNU General Public License](LICENSE).