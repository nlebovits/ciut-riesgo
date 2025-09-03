import matplotlib.pyplot as plt
import contextily as ctx


from io import BytesIO
from owslib.wfs import WebFeatureService
from shapely.geometry import box
import geopandas as gpd
import pandas as pd
import numpy as np
import boto3
import duckdb
import s2sphere
from botocore.config import Config



from matplotlib_map_utils import north_arrow, ScaleBar





WEB_MERCATOR_CRS = "EPSG:3857"  # visualization
WGS84_CRS = "EPSG:4326"  # for API calls


# Basic visualization settings (only for repeated values)
DEFAULT_FIGSIZE = (12, 10)
MAP_PADDING = 500
PLASMA_CMAP = plt.cm.plasma

# Color schemes for visualization
PELIGROSIDAD_COLORS = {
    "alta": PLASMA_CMAP(0.8),
    "media": PLASMA_CMAP(0.5),
}


def setup_base_map(
    use_crs, figsize=None, bounds=None, boundary_gdf=None, padding_x=None, padding_y=None
):
    """Create figure and set up basic map boundaries with padding."""
    if figsize is None:
        figsize = DEFAULT_FIGSIZE
    if padding_x is None:
        padding_x = MAP_PADDING
    if padding_y is None:
        padding_y = MAP_PADDING

    if bounds is None and boundary_gdf is not None:
        bounds = boundary_gdf.total_bounds

    # Convert bounds to Web Mercator for basemap compatibility
    if bounds is not None:
        # Create a temporary GeoDataFrame with the bounds to reproject
        temp_bounds = gpd.GeoDataFrame(
            geometry=[box(bounds[0], bounds[1], bounds[2], bounds[3])], crs=use_crs
        )
        bounds_3857 = temp_bounds.to_crs(WEB_MERCATOR_CRS).total_bounds
    else:
        bounds_3857 = bounds

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(bounds_3857[0] - padding_x, bounds_3857[2] + padding_x)
    ax.set_ylim(bounds_3857[1] - padding_y, bounds_3857[3] + padding_y)
    return fig, ax


def add_basemap(ax, zoom=13):
    """Add CartoDB basemap to the axes."""

    ctx.add_basemap(
        ax,
        source=ctx.providers.CartoDB.PositronNoLabels,
        zorder=0,
        zoom=zoom,
    )

    return ax


ScaleBar.set_size("xs")


def add_scale_bar_and_north_arrow(
    ax, location="upper right", scale_color="black", arrow_color="black"
):
    """Add a scale bar and north arrow to the map using matplotlib_map_utils."""
    # Add scale bar using matplotlib_map_utils ScaleBar class with ticks style
    scalebar = ScaleBar(
        location="upper left",
        style="ticks",
        bar={
            "projection": "EPSG:3857",
            "tickcolors": scale_color,
            "basecolors": scale_color,
            "minor_type": "none",
        },
        labels={"style": "first_last"},
    )
    ax.add_artist(scalebar)

    # Add north arrow using matplotlib_map_utils
    north_arrow(
        ax,
        location=location,
        scale=0.3,  # Small size
        rotation={"degrees": 0},
        base={"facecolor": "none", "edgecolor": arrow_color, "linewidth": 1},
        fancy=True,
        shadow=True,
        label=False,  # Hide the "N" text
    )


def add_boundary_outline(ax, boundary_gdf, crs="EPSG:3857"):
    """Add the outline of a boundary geodataframe to a map."""
    boundary_3857 = boundary_gdf.to_crs(crs)
    boundary_3857.plot(
        ax=ax,
        facecolor="none",
        edgecolor="black",
        linewidth=0.5,
        linestyle="--",
        legend=False,
        zorder=5,
    )


def create_consistent_map(title, crs, boundary_gdf=None, bounds=None):
    """Create a map with consistent styling and basemap."""
    fig, ax = setup_base_map(crs, bounds=bounds, boundary_gdf=boundary_gdf)

    add_basemap(ax)

    add_scale_bar_and_north_arrow(ax)

    add_boundary_outline(ax, boundary_gdf)

    ax.set_title(title, fontsize=16, fontweight="bold", pad=20)

    ax.set_axis_off()

    return fig, ax


def wfs_to_gdf(
    wfs_url: str, layer_name: str, srs: str = "EPSG:4326"
) -> gpd.GeoDataFrame:
    """
    Descarga una capa WFS y la devuelve como GeoDataFrame.

    Args:
        wfs_url (str): URL del servicio WFS.
        layer_name (str): Nombre de la capa (typename).
        srs (str): Código EPSG del sistema de referencia de coordenadas.

    Returns:
        gpd.GeoDataFrame: Capa descargada como GeoDataFrame.
    """
    wfs = WebFeatureService(url=wfs_url, version="2.0.0")
    response = wfs.getfeature(typename=layer_name, srsname=srs)
    gdf = gpd.read_file(BytesIO(response.read()))
    return gdf


def fetch_buildings(geodataframe, temp_file="buildings_filtered.parquet"):
    """Fetch building data for a given GeoDataFrame region"""

    # Get S2 cell and bounds
    center = geodataframe.to_crs(WEB_MERCATOR_CRS).union_all().centroid
    center_wgs84 = (
        gpd.GeoDataFrame(geometry=[center], crs=WEB_MERCATOR_CRS)
        .to_crs(WGS84_CRS)
        .geometry.iloc[0]
    )
    cell = s2sphere.CellId.from_lat_lng(
        s2sphere.LatLng.from_degrees(center_wgs84.y, center_wgs84.x)
    ).parent(10)
    bounds = geodataframe.to_crs(WGS84_CRS).total_bounds

    # Find matching S2 partition
    s3 = boto3.client(
        "s3",
        endpoint_url="https://data.source.coop",
        aws_access_key_id="",
        aws_secret_access_key="",
        config=Config(s3={"addressing_style": "path"}),
    )

    partitions = {
        obj["Key"].split("/")[-1].replace(".parquet", "")
        for obj in s3.list_objects_v2(
            Bucket="vida",
            Prefix="google-microsoft-osm-open-buildings/geoparquet/by_country_s2/country_iso=ARG/",
        ).get("Contents", [])
    }

    parent_id = next(
        str(cell.parent(level).id())
        for level in range(10, 0, -1)
        if str(cell.parent(level).id()) in partitions
    )

    # Setup DuckDB and query
    con = duckdb.connect()
    for cmd in [
        "INSTALL spatial",
        "LOAD spatial",
        "INSTALL httpfs",
        "LOAD httpfs",
        "SET s3_region='us-east-1'",
        "SET s3_endpoint='data.source.coop'",
        "SET s3_use_ssl=true",
        "SET s3_url_style='path'",
    ]:
        con.execute(cmd)

    # Export and read back
    query = f"""
    COPY (SELECT * FROM 's3://vida/google-microsoft-osm-open-buildings/geoparquet/by_country_s2/country_iso=ARG/{parent_id}.parquet'
          WHERE bbox.xmax >= {bounds[0]} AND bbox.xmin <= {bounds[2]} AND
                bbox.ymax >= {bounds[1]} AND bbox.ymin <= {bounds[3]}
    ) TO '{temp_file}' (FORMAT PARQUET);
    """

    con.execute(query)
    df = pd.read_parquet(temp_file)
    df["geometry"] = gpd.GeoSeries.from_wkb(df["geometry"])

    return gpd.GeoDataFrame(df, geometry="geometry", crs=WGS84_CRS)