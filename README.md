# CIUT Riesgo

Un código base para analizar el riesgo de inundación en el Partido de La Plata.

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