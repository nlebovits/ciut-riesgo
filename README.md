# CIUT Riesgo

Un código base para analizar el riesgo de inundación en el Partido de La Plata.

## Software Requerido

Cambio!

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

## Licencia

Este proyecto está licenciado bajo [GNU General Public License](LICENSE).