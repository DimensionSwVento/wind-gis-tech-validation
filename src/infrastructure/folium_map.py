"""
Folium map generation infrastructure.

This module provides utilities for creating interactive maps
using the Folium library for wind suitability visualization.
"""

import logging
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import folium
import folium.plugins
from folium import plugins
import json
import os
from pathlib import Path
import pandas as pd
import geopandas as gpd

try:
    import rasterio
    from rasterio.transform import from_bounds
    RASTERIO_AVAILABLE = True
except ImportError:
    RASTERIO_AVAILABLE = False
    rasterio = None
    from_bounds = None


class FoliumMapGenerator:
    """Generator for interactive Folium maps."""
    
    def __init__(self):
        """Initialize the Folium map generator."""
        self.logger = logging.getLogger(__name__)
    
    def create_wsi_map(self, 
                      wsi_data: np.ndarray,
                      aoi_geometry: Union[str, dict],
                      crs: str,
                      output_path: str,
                      title: str = "Wind Suitability Index",
                      center: Optional[Tuple[float, float]] = None,
                      zoom_start: int = 10) -> None:
        """
        Create interactive WSI map using Folium.
        
        Args:
            wsi_data: WSI raster data
            aoi_geometry: AOI geometry (file path or dict)
            crs: Coordinate reference system
            output_path: Output HTML file path
            title: Map title
            center: Optional map center (lat, lon)
            zoom_start: Initial zoom level
        """
        try:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Get map center
            if center is None:
                center = self._get_geometry_center(aoi_geometry)
            
            # Create base map with default tile
            m = folium.Map(
                location=center,
                zoom_start=zoom_start,
                tiles='OpenStreetMap'
            )
            
            # Add different tile layers
            self._add_tile_layers(m)
            
            # Add satellite layer
            folium.TileLayer(
                tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
                attr='Esri',
                name='Satelital',
                overlay=False,
                control=True
            ).add_to(m)
            
            # Add terrain layer
            folium.TileLayer(
                tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}',
                attr='Esri',
                name='Terreno',
                overlay=False,
                control=True
            ).add_to(m)
            
            
            # Add WSI visualization
            self._add_wsi_layer(m, wsi_data, aoi_geometry, crs)
            
            # Optional overlays
            self._add_hillshade_tiles(m)
            self._add_wind_raster_overlay(m)
            self._add_wms_layers(m)

            # Add AOI boundary
            self._add_aoi_boundary(m, aoi_geometry)
            
            # Add custom markers for important locations
            self._add_custom_markers(m)
            
            # Add candidate sites (if available)
            self._add_candidate_sites_cluster(m)
            
            # Add reference layers (restrictions, grid lines) if available
            self._add_reference_layers(m)
            
            # Add legend
            self._add_wsi_legend(m)
            
            # Add title
            self._add_title(m, title)
            
            # Add plugins
            self._add_plugins(m)
            
            # Add layer control
            folium.LayerControl(
                position='topright',
                collapsed=True,
                autoZIndex=True
            ).add_to(m)
            
            # Save map
            m.save(output_path)
            
            self.logger.info(f"Interactive map saved to: {output_path}")
            
        except Exception as e:
            self.logger.error(f"Failed to create WSI map: {e}")
            raise
    
    def _get_geometry_center(self, geometry_input: Union[str, dict]) -> Tuple[float, float]:
        """Get center point of geometry."""
        try:
            if isinstance(geometry_input, str):
                import geopandas as gpd
                gdf = gpd.read_file(geometry_input)
                if gdf.crs and gdf.crs.is_geographic:
                    gdf_projected = gdf.to_crs('EPSG:3857')
                    centroid = gdf_projected.geometry.centroid.iloc[0]
                    centroid_geo = gdf_projected.to_crs('EPSG:4326').geometry.centroid.iloc[0]
                    return (centroid_geo.y, centroid_geo.x)
                else:
                    centroid = gdf.geometry.centroid.iloc[0]
                    return (centroid.y, centroid.x)
            else:
                from shapely.geometry import shape
                geom = shape(geometry_input)
                centroid = geom.centroid
                return (centroid.y, centroid.x)
                
        except Exception as e:
            self.logger.warning(f"Could not get geometry center: {e}")
            return (11.5, -72.0)
    
    def _add_wsi_layer(self, 
                      map_obj: folium.Map, 
                      wsi_data: np.ndarray,
                      aoi_geometry: Union[str, dict],
                      crs: str) -> None:
        """Add WSI layer to map with real geometry filtering."""
        try:
            bounds = self._get_geometry_bounds(aoi_geometry)
            
            # Load Colombia geometries
            import sys
            import os
            sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
            from colombia_geometries import get_departamento_geometry, LA_GUAJIRA_GEOMETRY, is_point_in_geometry
            
            # Start with La Guajira as default, but will be updated by region selector
            geometry_filter = LA_GUAJIRA_GEOMETRY
            
            # Add sample points filtered by real geometry
            sample_points = self._sample_wsi_points(wsi_data, bounds, geometry_filter)
            
            # Create feature groups for different WSI levels
            excellent_group = folium.FeatureGroup(name='Excelente (0.8-1.0)')
            high_group = folium.FeatureGroup(name='Alto (0.6-0.8)')
            good_group = folium.FeatureGroup(name='Bueno (0.4-0.6)')
            regular_group = folium.FeatureGroup(name='Regular (0.2-0.4)')
            low_group = folium.FeatureGroup(name='Bajo (0.0-0.2)')
            
            # Create a master group for all terrestrial points
            terrestrial_group = folium.FeatureGroup(name='Puntos Terrestres')
            
            for point in sample_points:
                lat, lon, wsi_value = point
                color = self._get_wsi_color(wsi_value)
                radius = self._get_wsi_radius(wsi_value)
                
                popup_text = f"""
                <div style="font-family: Arial, sans-serif; min-width: 200px;">
                    <h4 style="margin: 0 0 10px 0; color: #333;">Punto de Idoneidad Eólica</h4>
                    <p style="margin: 5px 0;"><strong>WSI:</strong> {wsi_value:.3f}</p>
                    <p style="margin: 5px 0;"><strong>Latitud:</strong> {lat:.4f}°</p>
                    <p style="margin: 5px 0;"><strong>Longitud:</strong> {lon:.4f}°</p>
                    <p style="margin: 5px 0;"><strong>Clasificación:</strong> {self._get_wsi_classification(wsi_value)}</p>
                    <p style="margin: 5px 0; font-size: 12px; color: #666;">
                        Área terrestre de Colombia - La Guajira
                    </p>
                </div>
                """
                
                # Create circle marker (original style)
                marker = folium.CircleMarker(
                    location=[lat, lon],
                    radius=radius,
                    popup=folium.Popup(popup_text, max_width=250),
                    color='white',
                    weight=3,
                    fill=True,
                    fillColor=color,
                    fillOpacity=0.9
                )
                
                # Add to appropriate group based on WSI value
                if wsi_value >= 0.8:
                    marker.add_to(excellent_group)
                elif wsi_value >= 0.6:
                    marker.add_to(high_group)
                elif wsi_value >= 0.4:
                    marker.add_to(good_group)
                elif wsi_value >= 0.2:
                    marker.add_to(regular_group)
                else:
                    marker.add_to(low_group)
                
                # Also add to terrestrial group
                marker.add_to(terrestrial_group)
            
            # Add La Guajira boundary
            self._add_la_guajira_boundary(map_obj, geometry_filter)
            
            # Add department boundaries
            self._add_department_boundaries(map_obj)
            
            # Add all groups to map
            terrestrial_group.add_to(map_obj)  # Add terrestrial group first
            excellent_group.add_to(map_obj)
            high_group.add_to(map_obj)
            good_group.add_to(map_obj)
            regular_group.add_to(map_obj)
            low_group.add_to(map_obj)
                
        except Exception as e:
            self.logger.error(f"Failed to add WSI layer: {e}")
    
    def _add_la_guajira_boundary(self, map_obj: folium.Map, geometry: dict) -> None:
        """Add La Guajira boundary to map."""
        try:
            # Create boundary feature group
            boundary_group = folium.FeatureGroup(name='Límite La Guajira')
            
            # Add boundary polygon
            folium.GeoJson(
                geometry,
                style_function=lambda x: {
                    'fillColor': 'transparent',
                    'color': '#2E86AB',
                    'weight': 3,
                    'opacity': 0.8,
                    'dashArray': '10, 5'
                },
                popup=folium.Popup('Límite del Departamento de La Guajira', max_width=200)
            ).add_to(boundary_group)
            
            # Add boundary group to map
            boundary_group.add_to(map_obj)
            
        except Exception as e:
            self.logger.warning(f"Could not add La Guajira boundary: {e}")
    
    def _add_department_boundaries(self, map_obj: folium.Map) -> None:
        """Add department boundaries to map."""
        try:
            # Load Colombia geometries
            import sys
            import os
            sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
            from colombia_geometries import DEPARTAMENTOS
            
            # Create department boundaries feature group
            boundaries_group = folium.FeatureGroup(name='Límites Departamentos')
            
            # Add each department boundary
            for dept_name, dept_data in DEPARTAMENTOS.items():
                geometry = dept_data.get("geometry")
                if geometry:
                    folium.GeoJson(
                        geometry,
                        style_function=lambda x, name=dept_name: {
                            'fillColor': 'transparent',
                            'color': '#2E86AB',
                            'weight': 2,
                            'opacity': 0.8,
                            'fillOpacity': 0.1
                        },
                        popup=folium.Popup(f"<b>{dept_name}</b><br>Departamento de Colombia", max_width=200)
                    ).add_to(boundaries_group)
            
            # Add boundaries group to map
            boundaries_group.add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add department boundaries: {e}")
    
    def _get_geometry_bounds(self, geometry_input: Union[str, dict]) -> Tuple[float, float, float, float]:
        """Get bounds from geometry."""
        try:
            if isinstance(geometry_input, str):
                import geopandas as gpd
                gdf = gpd.read_file(geometry_input)
                return gdf.total_bounds
            else:
                from shapely.geometry import shape
                geom = shape(geometry_input)
                return geom.bounds
                
        except Exception as e:
            self.logger.warning(f"Could not get geometry bounds: {e}")
            return (-180, -90, 180, 90)
    
    def _sample_wsi_points(self, 
                          wsi_data: np.ndarray, 
                          bounds: Tuple[float, float, float, float],
                          geometry_filter=None) -> List[Tuple[float, float, float]]:
        """Sample points from WSI data, filtering for real geometry areas."""
        points = []
        minx, miny, maxx, maxy = bounds
        height, width = wsi_data.shape
        max_points = 2000  # Increase max points for better coverage
        
        if height * width <= max_points:
            step = 1
        else:
            step = max(1, int(np.sqrt(height * width / max_points)))
        
        for i in range(0, height, step):
            for j in range(0, width, step):
                wsi_value = wsi_data[i, j]
                if not np.isnan(wsi_value) and wsi_value > 0:
                    lat = miny + (maxy - miny) * (i / height)
                    lon = minx + (maxx - minx) * (j / width)
                    
                    # Use geometry filter if provided, otherwise use terrestrial detection
                    if geometry_filter:
                        # Import the function from colombia_geometries
                        import sys
                        import os
                        sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
                        from colombia_geometries import is_point_in_geometry
                        if is_point_in_geometry(lat, lon, geometry_filter):
                            points.append((lat, lon, float(wsi_value)))
                    else:
                        if self._is_terrestrial_point_improved(lat, lon):
                            points.append((lat, lon, float(wsi_value)))
        
        if len(points) > max_points:
            import random
            points = random.sample(points, max_points)
        
        return points
    
    def _is_terrestrial_point(self, lat: float, lon: float) -> bool:
        """Check if a point is on land (terrestrial) based on Colombia's geography."""
        colombia_bounds = {
            'north': 15.5, 'south': 4.0, 'east': -66.0, 'west': -82.0
        }
        
        if not (colombia_bounds['south'] <= lat <= colombia_bounds['north'] and
                colombia_bounds['west'] <= lon <= colombia_bounds['east']):
            return False
        
        if lat > 12.0 and lon > -75.0:
            return False
        
        if lon < -78.0:
            return False
        
        if 10.5 <= lat <= 12.5 and -73.0 <= lon <= -71.0:
            return True
        
        return True
    
    def _is_terrestrial_point_improved(self, lat: float, lon: float) -> bool:
        """Check if a point is on land using precise La Guajira coastline detection."""
        
        # First check if point is within general La Guajira bounds
        if not (10.5 <= lat <= 12.2 and -73.5 <= lon <= -70.5):
            return False
        
        # Define the precise coastline of La Guajira
        # These coordinates represent the actual land-sea boundary based on real geography
        
        # Northern boundary - exclude Caribbean Sea (very restrictive)
        if lat > 11.7:
            return False
        
        # Eastern boundary - exclude Caribbean Sea (very restrictive)
        if lon > -71.5:
            return False
        
        # Define specific marine areas to exclude based on real coastline
        # These are areas that are definitely in the Caribbean Sea
        
        # Exclude northernmost areas (Punta Gallinas region)
        if lat > 11.6 and lon > -72.2:
            return False
        
        # Exclude Cabo de la Vela area
        if lat > 11.4 and lon > -72.0:
            return False
        
        # Exclude Riohacha coastal area
        if lat > 11.2 and lon > -71.8:
            return False
        
        # Exclude Maicao coastal area
        if lat > 11.0 and lon > -71.6:
            return False
        
        # Exclude areas that are clearly in the Caribbean Sea
        # These coordinates are based on the actual coastline
        if lat > 11.5 and lon > -71.9:
            return False
        
        if lat > 11.3 and lon > -71.7:
            return False
        
        if lat > 11.1 and lon > -71.5:
            return False
        
        # Exclude areas that are clearly marine based on real geography
        # Northern coastal exclusion
        if lat > 11.6 and lon > -72.5:
            return False
        
        # Eastern coastal exclusion
        if lon > -71.7 and lat > 11.3:
            return False
        
        if lon > -71.5 and lat > 11.1:
            return False
        
        if lon > -71.3 and lat > 10.9:
            return False
        
        # Include main inland areas of La Guajira (very restrictive)
        if 10.8 <= lat <= 11.5 and -73.0 <= lon <= -72.0:
            return True
        
        # Include some coastal areas but be extremely restrictive
        if 11.0 <= lat <= 11.3 and -72.5 <= lon <= -72.0:
            # Additional coastal checks to ensure we're on land
            if lat > 11.2 and lon > -72.2:
                return False
            if lat > 11.1 and lon > -72.1:
                return False
            return True
        
        # Default to False for remaining areas (be very conservative)
        return False
    
    def _is_terrestrial_point_api(self, lat: float, lon: float) -> bool:
        """Check if a point is on land using a marine detection API."""
        try:
            import requests
            import time
            
            # Use a free marine detection API
            # This is a simple approach using reverse geocoding
            url = f"https://api.bigdatacloud.net/data/reverse-geocode-client?latitude={lat}&longitude={lon}&localityLanguage=en"
            
            response = requests.get(url, timeout=5)
            if response.status_code == 200:
                data = response.json()
                # Check if the location is on land
                if 'locality' in data and data['locality']:
                    # If we get a locality name, it's likely on land
                    return True
                elif 'principalSubdivision' in data and data['principalSubdivision']:
                    # If we get a subdivision, it's likely on land
                    return True
                else:
                    # If no locality or subdivision, it might be water
                    return False
            else:
                # If API fails, fall back to coordinate-based detection
                return self._is_terrestrial_point_improved(lat, lon)
                
        except Exception as e:
            # If API fails, fall back to coordinate-based detection
            self.logger.warning(f"API check failed for {lat}, {lon}: {e}")
            return self._is_terrestrial_point_improved(lat, lon)
    
    def _is_point_in_geometry(self, lat: float, lon: float, geometry: dict) -> bool:
        """Check if a point is inside a polygon geometry."""
        if not geometry or geometry.get("type") != "Polygon":
            return False
        
        coords = geometry["coordinates"][0]
        x, y = lon, lat
        
        n = len(coords)
        inside = False
        
        p1x, p1y = coords[0]
        for i in range(1, n + 1):
            p2x, p2y = coords[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        
        return inside
    
    def _get_wsi_color(self, wsi_value: float) -> str:
        """Get color for WSI value."""
        if wsi_value < 0.2:
            return '#FF0000'
        elif wsi_value < 0.4:
            return '#FFA500'
        elif wsi_value < 0.6:
            return '#FFFF00'
        elif wsi_value < 0.8:
            return '#90EE90'
        else:
            return '#006400'
    
    def _get_wsi_radius(self, wsi_value: float) -> int:
        """Get radius for WSI value."""
        return int(4 + wsi_value * 8)
    
    def _get_wsi_classification(self, wsi_value: float) -> str:
        """Get Spanish classification for WSI value."""
        if wsi_value < 0.2:
            return "Bajo"
        elif wsi_value < 0.4:
            return "Regular"
        elif wsi_value < 0.6:
            return "Bueno"
        elif wsi_value < 0.8:
            return "Alto"
        else:
            return "Excelente"
    
    def _add_aoi_boundary(self, map_obj: folium.Map, aoi_geometry: Union[str, dict]) -> None:
        """
        Add AOI boundary to map.
        
        Args:
            map_obj: Folium map object
            aoi_geometry: AOI geometry
        """
        try:
            if isinstance(aoi_geometry, str):
                import geopandas as gpd
                gdf = gpd.read_file(aoi_geometry)
                geom = gdf.geometry.iloc[0]
            else:
                from shapely.geometry import shape
                geom = shape(aoi_geometry)
            
            # Convert geometry to GeoJSON
            geojson = geom.__geo_interface__
            
            # Add to map
            folium.GeoJson(
                geojson,
                style_function=lambda x: {
                    'fillColor': 'transparent',
                    'color': 'blue',
                    'weight': 3,
                    'opacity': 0.8
                }
            ).add_to(map_obj)
            
        except Exception as e:
            self.logger.warning(f"Could not add AOI boundary: {e}")
    
    def _add_wsi_legend(self, map_obj: folium.Map) -> None:
        """
        Add WSI legend to map with colors only.
        
        Args:
            map_obj: Folium map object
        """
        legend_html = '''
        <div style="position: fixed; 
                    bottom: 50px; left: 50px; width: 300px; height: 180px; 
                    background-color: white; border:2px solid #333; z-index:9999; 
                    font-size:14px; padding: 15px; border-radius: 8px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.3); overflow: visible;">
        <p style="margin: 0 0 10px 0; font-weight: bold; text-align: center; color: #333;">
            Índice de Idoneidad Eólica
        </p>
        <div style="display: flex; align-items: center; margin: 6px 0;">
            <span style="color:#FF0000; font-size:16px; margin-right: 8px;">●</span>
            <span>0.0 - 0.2 (Bajo)</span>
        </div>
        <div style="display: flex; align-items: center; margin: 6px 0;">
            <span style="color:#FFA500; font-size:16px; margin-right: 8px;">●</span>
            <span>0.2 - 0.4 (Regular)</span>
        </div>
        <div style="display: flex; align-items: center; margin: 6px 0;">
            <span style="color:#FFFF00; font-size:16px; margin-right: 8px;">●</span>
            <span>0.4 - 0.6 (Bueno)</span>
        </div>
        <div style="display: flex; align-items: center; margin: 6px 0;">
            <span style="color:#90EE90; font-size:16px; margin-right: 8px;">●</span>
            <span>0.6 - 0.8 (Alto)</span>
        </div>
        <div style="display: flex; align-items: center; margin: 6px 0;">
            <span style="color:#006400; font-size:16px; margin-right: 8px;">●</span>
            <span>0.8 - 1.0 (Excelente)</span>
        </div>
        </div>
        '''
        
        map_obj.get_root().html.add_child(folium.Element(legend_html))
    
    def _add_tile_layers(self, map_obj: folium.Map) -> None:
        """
        Add multiple tile layers to the map.
        
        Args:
            map_obj: Folium map object
        """
        # OpenStreetMap (default)
        folium.TileLayer(
            tiles='OpenStreetMap',
            name='OpenStreetMap',
            overlay=False,
            control=True
        ).add_to(map_obj)
        
        # Add additional useful layers
        folium.TileLayer(
            tiles='https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png',
            attr='OpenTopoMap',
            name='Topográfico',
            overlay=False,
            control=True
        ).add_to(map_obj)
        
        # Carto basemaps
        folium.TileLayer(
            tiles='https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png',
            attr='Carto',
            name='Carto Positron',
            overlay=False,
            control=True
        ).add_to(map_obj)
        folium.TileLayer(
            tiles='https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png',
            attr='Carto',
            name='Carto DarkMatter',
            overlay=False,
            control=True
        ).add_to(map_obj)
        # Labels overlay
        folium.TileLayer(
            tiles='https://{s}.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}.png',
            attr='Carto',
            name='Etiquetas',
            overlay=True,
            control=True
        ).add_to(map_obj)
        
        # Stamen
        folium.TileLayer(
            tiles='https://stamen-tiles-{s}.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png',
            attr='Stamen Terrain',
            name='Stamen Terrain',
            overlay=False,
            control=True
        ).add_to(map_obj)
        folium.TileLayer(
            tiles='https://stamen-tiles-{s}.a.ssl.fastly.net/toner/{z}/{x}/{y}.png',
            attr='Stamen Toner',
            name='Stamen Toner',
            overlay=False,
            control=True
        ).add_to(map_obj)
        
    
    def _add_title(self, map_obj: folium.Map, title: str) -> None:
        """
        Add title to map without logo.
        
        Args:
            map_obj: Folium map object
            title: Map title
        """
        # Default to Spanish title if not provided
        if title == "Wind Suitability Index":
            title = "Índice de Idoneidad Eólica"
        
        title_html = f'''
        <div style="position: fixed; 
                    top: 10px; left: 50px; width: 300px; height: 80px; 
                    background-color: #f8f9fa; color: #333; border:2px solid #dee2e6; z-index:9999; 
                    font-size:18px; padding: 20px; text-align: center;
                    border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1)">
        <div>
            <b>{title}</b><br>
            <small style="color: #666;">Vento</small>
        </div>
        </div>
        '''
        
        map_obj.get_root().html.add_child(folium.Element(title_html))
    
    def _add_plugins(self, map_obj: folium.Map) -> None:
        """
        Add useful GIS plugins to map.
        
        Args:
            map_obj: Folium map object
        """
        try:
        # Add fullscreen plugin
            plugins.Fullscreen(
                position='topright',
                title='Pantalla Completa',
                title_cancel='Salir de Pantalla Completa',
                force_separate_button=True
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add fullscreen plugin: {e}")
        
        
        try:
            # Add draw plugin for GIS editing
            plugins.Draw(
                position='topleft',
                export=False,
                filename='vento_gis_data',
                show_geometry_on_click=True,
                draw_options={
                    'polyline': {
                        'allowIntersection': False,
                        'drawError': {'color': '#e1e100', 'message': '<strong>Error:</strong> Las líneas no pueden cruzarse!'},
                        'shapeOptions': {'color': '#bada55'}
                    },
                    'polygon': {
                        'allowIntersection': False,
                        'drawError': {'color': '#e1e100', 'message': '<strong>Error:</strong> Los polígonos no pueden cruzarse!'},
                        'shapeOptions': {'color': '#bada55'}
                    },
                    'circle': False,
                    'rectangle': {
                        'shapeOptions': {'color': '#bada55'}
                    },
                    'marker': {
                        'shapeOptions': {'color': '#bada55'}
                    }
                },
                edit_options={
                    'featureGroup': None,
                    'edit': True,
                    'remove': True
                }
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add draw plugin: {e}")
        
        try:
            # Add measure control (distance/area)
            plugins.MeasureControl(
                position='topleft',
                primary_length_unit='kilometers',
                secondary_length_unit='meters',
                primary_area_unit='hectares',
                secondary_area_unit='sqmeters'
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add measure control: {e}")

        try:
            # Add geolocation control
            plugins.LocateControl(
                position='topleft',
                flyTo=True,
                keepCurrentZoomLevel=False,
                strings={
                    'title': 'Mi ubicación'
                }
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add locate control: {e}")

        try:
            # Add export plugin for PNG export
            self._add_export_plugin(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add export plugin: {e}")
        
        try:
            # Add minimap
            plugins.MiniMap(
                position='bottomright',
                width=200,
                height=150,
                collapsed_width=25,
                collapsed_height=25,
                zoom_level_offset=-5,
                zoom_level_fixed=False,
                center_fixed=False,
                zoom_control=True,
                auto_toggle=False
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add minimap plugin: {e}")
        
        try:
            # Add mouse position
            plugins.MousePosition(
                position='bottomleft',
                separator=' | ',
                empty_string='Coordenadas no disponibles',
                lng_first=True,
                num_digits=5,
                prefix='',
                lat_formatter=None,
                lng_formatter=None
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add mouse position plugin: {e}")
        
        try:
            # Add region selector
            self._add_region_selector(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add region selector: {e}")

        try:
            # Add buffer tools (5km/10km) at last clicked point
            self._add_buffer_tools(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add buffer tools: {e}")

    def _add_hillshade_tiles(self, map_obj: folium.Map) -> None:
        """Add Esri WorldHillshade as an overlay."""
        try:
            folium.TileLayer(
                tiles='https://server.arcgisonline.com/ArcGIS/rest/services/Elevation/World_Hillshade/MapServer/tile/{z}/{y}/{x}',
                attr='Esri',
                name='Sombreado (Hillshade)',
                overlay=True,
                control=True,
                opacity=0.6
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add hillshade tiles: {e}")

    def _add_wind_raster_overlay(self, map_obj: folium.Map) -> None:
        """Add wind raster as semi-transparent overlay if available."""
        try:
            if not RASTERIO_AVAILABLE:
                return
            wind_path = Path('data/raw/wind_100m.tif')
            if not wind_path.exists():
                return
            with rasterio.open(wind_path) as src:
                bounds = src.bounds
                data = src.read(1, masked=True)
                # Normalize 0-1 for display if possible
                import numpy as np
                arr = data.astype('float32')
                vmin = np.nanpercentile(arr, 2)
                vmax = np.nanpercentile(arr, 98)
                arr = (arr - vmin) / (vmax - vmin + 1e-6)
                arr = np.clip(arr, 0, 1)
                # Create simple PNG via matplotlib in-memory
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from io import BytesIO
                fig, ax = plt.subplots(figsize=(6, 6), dpi=100)
                ax.axis('off')
                ax.imshow(arr, cmap='viridis', vmin=0, vmax=1)
                buf = BytesIO()
                plt.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
                plt.close(fig)
                buf.seek(0)
                import base64
                png_b64 = base64.b64encode(buf.read()).decode('ascii')
                image_url = f"data:image/png;base64,{png_b64}"
                folium.raster_layers.ImageOverlay(
                    name='Viento 100m (ref)',
                    image=image_url,
                    bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
                    opacity=0.45,
                    interactive=False,
                    cross_origin=False
                ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add wind raster overlay: {e}")

    def _add_candidate_sites_cluster(self, map_obj: folium.Map) -> None:
        """Add candidate sites as a clustered point layer if the GPKG exists."""
        try:
            gpkg_path = Path('outputs/vectors/candidate_sites.gpkg')
            if not gpkg_path.exists():
                return
            gdf = gpd.read_file(gpkg_path)
            if gdf.empty:
                return
            # Ensure WGS84
            if gdf.crs and gdf.crs.to_string() != 'EPSG:4326':
                gdf = gdf.to_crs('EPSG:4326')
            cluster = plugins.MarkerCluster(name='Sitios Candidatos')
            for _, row in gdf.iterrows():
                geom = row.geometry
                if geom is None or geom.is_empty:
                    continue
                if geom.geom_type == 'Point':
                    lat, lon = geom.y, geom.x
                else:
                    centroid = geom.centroid
                    lat, lon = centroid.y, centroid.x
                wsi = row.get('wsi', row.get('WSI', None))
                popup_parts = [f"<b>Sitio candidato</b>"]
                if wsi is not None and not pd.isna(wsi):
                    popup_parts.append(f"WSI: {float(wsi):.3f}")
                popup_html = '<br>'.join(popup_parts)
                folium.Marker(
                    location=[lat, lon],
                    popup=folium.Popup(popup_html, max_width=250),
                    icon=folium.Icon(color='green', icon='bolt', prefix='fa')
                ).add_to(cluster)
            cluster.add_to(map_obj)

            # Add search over candidate sites (hidden GeoJSON layer used by Search plugin)
            try:
                geojson_features = []
                for _, row in gdf.iterrows():
                    geom = row.geometry
                    if geom is None or geom.is_empty:
                        continue
                    if geom.geom_type != 'Point':
                        geom = geom.centroid
                    wsi = row.get('wsi', row.get('WSI', None))
                    label = row.get('name') or row.get('id') or f"Sitio ({geom.y:.4f}, {geom.x:.4f})"
                    if wsi is not None and not pd.isna(wsi):
                        label = f"{label} | WSI {float(wsi):.3f}"
                    geojson_features.append({
                        "type": "Feature",
                        "properties": {"label": label},
                        "geometry": geom.__geo_interface__
                    })
                search_geojson = {
                    "type": "FeatureCollection",
                    "features": geojson_features
                }
                search_layer = folium.GeoJson(
                    search_geojson,
                    name='Sitios Candidatos (Buscar)',
                    show=False,
                    tooltip=folium.GeoJsonTooltip(fields=['label'], aliases=['Sitio'])
                )
                search_layer.add_to(map_obj)
                plugins.Search(
                    layer=search_layer,
                    search_label='label',
                    placeholder='Buscar sitio candidato...',
                    collapsed=False,
                    geom_type='Point',
                    search_zoom=13
                ).add_to(map_obj)
            except Exception as e:
                self.logger.warning(f"Could not add candidate sites search: {e}")
        except Exception as e:
            self.logger.warning(f"Could not add candidate sites cluster: {e}")

    def _add_reference_layers(self, map_obj: folium.Map) -> None:
        """Add optional reference layers (restrictions and grid lines) if present."""
        try:
            # Restrictions polygons
            restr_path = Path('data/raw/restrictions.gpkg')
            if restr_path.exists():
                try:
                    restr_gdf = gpd.read_file(restr_path)
                    if not restr_gdf.empty:
                        if restr_gdf.crs and restr_gdf.crs.to_string() != 'EPSG:4326':
                            restr_gdf = restr_gdf.to_crs('EPSG:4326')
                        layer = folium.FeatureGroup(name='Restricciones (Referencia)')
                        for _, row in restr_gdf.iterrows():
                            folium.GeoJson(
                                row.geometry.__geo_interface__,
                                style_function=lambda x: {
                                    'fillColor': '#ff0000',
                                    'color': '#b22222',
                                    'weight': 1,
                                    'fillOpacity': 0.15,
                                    'opacity': 0.7
                                }
                            ).add_to(layer)
                        layer.add_to(map_obj)
                except Exception as e:
                    self.logger.warning(f"Failed to add restrictions layer: {e}")

            # Grid lines or transmission
            grid_path = Path('data/raw/grid_lines.gpkg')
            if grid_path.exists():
                try:
                    grid_gdf = gpd.read_file(grid_path)
                    if not grid_gdf.empty:
                        if grid_gdf.crs and grid_gdf.crs.to_string() != 'EPSG:4326':
                            grid_gdf = grid_gdf.to_crs('EPSG:4326')
                        layer = folium.FeatureGroup(name='Red Eléctrica / Líneas (Ref)')
                        for _, row in grid_gdf.iterrows():
                            folium.GeoJson(
                                row.geometry.__geo_interface__,
                                style_function=lambda x: {
                                    'color': '#555555',
                                    'weight': 2,
                                    'opacity': 0.8
                                }
                            ).add_to(layer)
                        layer.add_to(map_obj)
                except Exception as e:
                    self.logger.warning(f"Failed to add grid lines layer: {e}")
        except Exception as e:
            self.logger.warning(f"Could not add reference layers: {e}")
    
    def _add_export_plugin(self, map_obj: folium.Map) -> None:
        """
        Add custom export plugin for PNG export.
        
        Args:
            map_obj: Folium map object
        """
        export_html = '''
        <div id="export-control" style="position: fixed; 
                    top: 10px; right: 10px; z-index: 1000;">
            <button id="export-png-btn" style="
                background: #6c757d;
                color: white;
                border: none;
                padding: 12px 20px;
                border-radius: 8px;
                font-size: 14px;
                font-weight: bold;
                cursor: pointer;
                box-shadow: 0 4px 15px rgba(0,0,0,0.2);
                transition: all 0.3s ease;
                display: flex;
                align-items: center;
                gap: 8px;
            " onmouseover="this.style.transform='translateY(-2px)'; this.style.boxShadow='0 6px 20px rgba(0,0,0,0.3)'; this.style.background='#5a6268'" 
               onmouseout="this.style.transform='translateY(0px)'; this.style.boxShadow='0 4px 15px rgba(0,0,0,0.2)'; this.style.background='#6c757d'"
               onclick="exportMapAsPNG()">
                Exportar PNG
            </button>
        </div>
        
        <script>
        function exportMapAsPNG() {
            const button = document.getElementById('export-png-btn');
            const originalText = button.innerHTML;
            
            // Show loading state
            button.innerHTML = 'Generando...';
            button.disabled = true;
            
            // Get map container
            const mapContainer = document.querySelector('.folium-map');
            if (!mapContainer) {
                alert('Error: No se pudo encontrar el mapa');
                button.innerHTML = originalText;
                button.disabled = false;
                return;
            }
            
            // Use html2canvas to capture the map
            if (typeof html2canvas === 'undefined') {
                // Load html2canvas if not available
                const script = document.createElement('script');
                script.src = 'https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js';
                script.onload = function() {
                    captureMap();
                };
                script.onerror = function() {
                    alert('Error: No se pudo cargar la librería de captura');
                    button.innerHTML = originalText;
                    button.disabled = false;
                };
                document.head.appendChild(script);
            } else {
                captureMap();
            }
            
            function captureMap() {
                try {
                    html2canvas(mapContainer, {
                        useCORS: true,
                        allowTaint: true,
                        backgroundColor: '#ffffff',
                        scale: 2,
                        logging: false,
                        width: mapContainer.offsetWidth,
                        height: mapContainer.offsetHeight
                    }).then(function(canvas) {
                        // Create download link
                        const link = document.createElement('a');
                        link.download = 'vento_wsi_map_' + new Date().toISOString().slice(0,19).replace(/:/g, '-') + '.png';
                        link.href = canvas.toDataURL('image/png');
                        
                        // Trigger download
                        document.body.appendChild(link);
                        link.click();
                        document.body.removeChild(link);
                        
                        // Reset button
                        button.innerHTML = 'Exportado';
                        setTimeout(() => {
                            button.innerHTML = originalText;
                            button.disabled = false;
                        }, 2000);
                        
                    }).catch(function(error) {
                        console.error('Error al capturar el mapa:', error);
                        alert('Error al generar la imagen PNG. Intente de nuevo.');
                        button.innerHTML = originalText;
                        button.disabled = false;
                    });
                } catch (error) {
                    console.error('Error en la captura:', error);
                    alert('Error al generar la imagen PNG. Intente de nuevo.');
                    button.innerHTML = originalText;
                    button.disabled = false;
                }
            }
        }
        </script>
        '''
        
        map_obj.get_root().html.add_child(folium.Element(export_html))
    
    def _add_region_selector(self, map_obj: folium.Map) -> None:
        """
        Add region selector for departments and municipalities.
        
        Args:
            map_obj: Folium map object
        """
        selector_html = '''
        <div id="region-selector" style="position: fixed; 
                    top: 10px; right: 10px; z-index: 1000;">
            <div style="
                background: white;
                border: 2px solid #dee2e6;
                border-radius: 6px;
                padding: 10px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                font-family: Arial, sans-serif;
                font-size: 12px;
                min-width: 180px;
                max-width: 200px;
            ">
                <h4 style="margin: 0 0 8px 0; color: #333; font-size: 14px;">Región</h4>
                
                <div style="margin-bottom: 8px;">
                    <label for="departamento-select" style="display: block; margin-bottom: 3px; font-weight: bold; font-size: 11px;">Depto:</label>
                    <select id="departamento-select" style="width: 100%; padding: 5px; border: 1px solid #ccc; border-radius: 4px;" onchange="updateMunicipios()">
                        <option value="La Guajira">La Guajira</option>
                        <option value="Cesar">Cesar</option>
                        <option value="Magdalena">Magdalena</option>
                        <option value="Atlántico">Atlántico</option>
                    </select>
                </div>
                
                <div style="margin-bottom: 8px;">
                    <label for="municipio-select" style="display: block; margin-bottom: 3px; font-weight: bold; font-size: 11px;">Municipio:</label>
                    <select id="municipio-select" style="width: 100%; padding: 5px; border: 1px solid #ccc; border-radius: 4px;" onchange="updateMunicipioZoom()">
                        <option value="todos">Todos los municipios</option>
                        <option value="Riohacha">Riohacha</option>
                        <option value="Maicao">Maicao</option>
                        <option value="Uribia">Uribia</option>
                        <option value="Manaure">Manaure</option>
                        <option value="San Juan del Cesar">San Juan del Cesar</option>
                    </select>
                </div>
                
                <button onclick="applyRegionFilter()" style="
                    width: 100%;
                    padding: 6px;
                    background: #90EE90;
                    color: #333;
                    border: 1px solid #7ED321;
                    border-radius: 4px;
                    cursor: pointer;
                    font-weight: bold;
                    font-size: 11px;
                ">Aplicar</button>
                
                <button onclick="resetRegionFilter()" style="
                    width: 100%;
                    padding: 6px;
                    background: #6c757d;
                    color: white;
                    border: none;
                    border-radius: 4px;
                    cursor: pointer;
                    font-weight: bold;
                    margin-top: 4px;
                    font-size: 11px;
                ">Reset</button>
            </div>
        </div>
        
        <script>
        const departamentos = {
            "La Guajira": ["Riohacha", "Maicao", "Uribia", "Manaure", "San Juan del Cesar", "Villanueva", "El Molino", "Fonseca", "Barrancas", "Distracción", "Hatonuevo", "La Jagua del Pilar", "Urumita", "Albania", "Dibula"],
            "Cesar": ["Valledupar", "Aguachica", "Codazzi", "La Paz", "San Diego", "Chimichagua", "Curumaní", "El Copey", "La Gloria", "Manaure Balcón del Cesar"],
            "Magdalena": ["Santa Marta", "Ciénaga", "Fundación", "Aracataca", "El Retén", "Zona Bananera", "Algarrobo", "Ariguaní", "Cerro San Antonio", "Chivolo"],
            "Atlántico": ["Barranquilla", "Soledad", "Malambo", "Sabanagrande", "Sabanalarga", "Puerto Colombia", "Galapa", "Tubará", "Usiacurí", "Luruaco"]
        };
        
        function updateMunicipios() {
            const departamento = document.getElementById('departamento-select').value;
            const municipioSelect = document.getElementById('municipio-select');
            
            // Clear existing options
            municipioSelect.innerHTML = '<option value="todos">Todos los municipios</option>';
            
            // Add municipios for selected departamento
            if (departamentos[departamento]) {
                departamentos[departamento].forEach(municipio => {
                    const option = document.createElement('option');
                    option.value = municipio;
                    option.textContent = municipio;
                    municipioSelect.appendChild(option);
                });
            }
        }
        
        function updateRegion() {
            // Auto-zoom when department changes
            const departamento = document.getElementById('departamento-select').value;
            
            // Define bounds for each department
            const departmentBounds = {
                "La Guajira": [[10.5, -73.5], [12.5, -70.5]],
                "Cesar": [[8.0, -74.0], [10.5, -72.0]],
                "Magdalena": [[9.0, -75.0], [11.0, -73.0]],
                "Atlántico": [[10.0, -75.5], [11.0, -74.0]]
            };
            
            // Get bounds for selected department
            const bounds = departmentBounds[departamento];
            if (bounds) {
                // Fit map to department bounds with smooth transition
                map.fitBounds(bounds, {
                    padding: [20, 20],
                    maxZoom: 12
                });
                
                // Show notification
                showNotification('Zoom automático a: ' + departamento);
            }
        }
        
        function updateMunicipioZoom() {
            const departamento = document.getElementById('departamento-select').value;
            const municipio = document.getElementById('municipio-select').value;
            
            if (municipio === 'todos') {
                // If "todos" selected, zoom to department
                updateRegion();
                return;
            }
            
            // Define specific coordinates for major municipalities
            const municipioCoords = {
                "La Guajira": {
                    "Riohacha": [11.5444, -72.9072],
                    "Maicao": [11.3833, -72.2333],
                    "Uribia": [11.7167, -72.2667],
                    "Manaure": [11.7833, -72.4500],
                    "San Juan del Cesar": [10.7667, -73.0000]
                },
                "Cesar": {
                    "Valledupar": [10.4631, -73.2532],
                    "Aguachica": [8.3081, -73.6161],
                    "Codazzi": [9.2333, -73.2333]
                },
                "Magdalena": {
                    "Santa Marta": [11.2408, -74.2119],
                    "Ciénaga": [11.0069, -74.2478],
                    "Fundación": [10.5208, -74.1858]
                },
                "Atlántico": {
                    "Barranquilla": [10.9639, -74.7964],
                    "Soledad": [10.9167, -74.7667],
                    "Malambo": [10.8667, -74.7667]
                }
            };
            
            const coords = municipioCoords[departamento] && municipioCoords[departamento][municipio];
            if (coords) {
                // Zoom to specific municipality
                map.setView(coords, 13, {
                    animate: true,
                    duration: 1.0
                });
                
                // Add marker for municipality
                if (window.municipioMarker) {
                    map.removeLayer(window.municipioMarker);
                }
                
                window.municipioMarker = L.marker(coords).addTo(map)
                    .bindPopup('<b>' + municipio + '</b><br>' + departamento)
                    .openPopup();
                
                showNotification('Zoom a: ' + municipio + ', ' + departamento);
            }
        }
        
        function applyRegionFilter() {
            const departamento = document.getElementById('departamento-select').value;
            const municipio = document.getElementById('municipio-select').value;
            
            // Define bounds for each department
            const departmentBounds = {
                "La Guajira": [[10.5, -73.5], [12.5, -70.5]],
                "Cesar": [[8.0, -74.0], [10.5, -72.0]],
                "Magdalena": [[9.0, -75.0], [11.0, -73.0]],
                "Atlántico": [[10.0, -75.5], [11.0, -74.0]]
            };
            
            // Get bounds for selected department
            const bounds = departmentBounds[departamento];
            if (bounds) {
                // Fit map to department bounds
                map.fitBounds(bounds, {padding: [20, 20]});
                
                // Hide all existing WSI point groups
                const groups = ['Excelente (0.8-1.0)', 'Alto (0.6-0.8)', 'Bueno (0.4-0.6)', 'Regular (0.2-0.4)', 'Bajo (0.0-0.2)', 'Puntos Terrestres'];
                groups.forEach(groupName => {
                    const group = map.layerManager._layers[groupName];
                    if (group) {
                        group.setVisible(false);
                    }
                });
                
                // Generate new points for selected region
                generateRegionPoints(departamento, municipio);
                
                // Show success message
                const message = 'Filtro aplicado a: ' + departamento + (municipio !== 'todos' ? ' - ' + municipio : '');
                showNotification(message);
            }
        }
        
        function generateRegionPoints(departamento, municipio) {
            // This function will be called to regenerate points based on selected region
            // For now, we'll show a message that the region has been selected
            showNotification('Región seleccionada: ' + departamento + (municipio !== 'todos' ? ' - ' + municipio : ''));
            
            // The actual point generation will be handled by the Python backend
            // when the map is regenerated with the new region filter
            // For now, we'll just show the notification
        }
        
        function showNotification(message) {
            // Create notification element
            const notification = document.createElement('div');
            notification.style.cssText = `
                position: fixed;
                top: 20px;
                left: 50%;
                transform: translateX(-50%);
                background: #28a745;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                z-index: 10000;
                font-family: Arial, sans-serif;
                font-size: 14px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.3);
            `;
            notification.textContent = message;
            
            // Add to page
            document.body.appendChild(notification);
            
            // Remove after 3 seconds
            setTimeout(() => {
                if (notification.parentNode) {
                    notification.parentNode.removeChild(notification);
                }
            }, 3000);
        }
        
        function resetRegionFilter() {
            document.getElementById('departamento-select').value = 'La Guajira';
            updateMunicipios();
            document.getElementById('municipio-select').value = 'todos';
            
            // Remove municipality marker if exists
            if (window.municipioMarker) {
                map.removeLayer(window.municipioMarker);
                window.municipioMarker = null;
            }
            
            // Zoom to La Guajira
            updateRegion();
            showNotification('Filtro reseteado - Zoom a La Guajira');

            // Clear drawn shapes created by Draw (heuristic: remove shapes with #bada55 color)
            if (map && map.eachLayer) {
                map.eachLayer(function(lyr){
                    if (lyr && lyr.options && (lyr.options.color === '#bada55' || lyr.options.fillColor === '#bada55')) {
                        try { map.removeLayer(lyr); } catch(e) {}
                    }
                });
            }

            // Reset export button label if present
            const btn = document.getElementById('export-png-btn');
            if (btn) {
                btn.innerHTML = 'Exportar PNG';
                btn.disabled = false;
            }
        }
        
        // Initialize on page load
        document.addEventListener('DOMContentLoaded', function() {
            updateMunicipios();
        });
        </script>
        '''
        
        map_obj.get_root().html.add_child(folium.Element(selector_html))

    def _add_wms_layers(self, map_obj: folium.Map) -> None:
        """Add public WMS overlays (roads/base)."""
        try:
            folium.raster_layers.WmsTileLayer(
                url='https://ows.terrestris.de/osm/service',
                name='WMS OSM (terrestris)',
                layers='OSM-WMS',
                fmt='image/png',
                transparent=True,
                overlay=True,
                control=True,
                attr='terrestris OSM'
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add terrestris OSM WMS: {e}")
        try:
            folium.raster_layers.WmsTileLayer(
                url='https://ows.terrestris.de/ows',
                name='SRTM30 (WMS)',
                layers='SRTM30-Colored',
                fmt='image/png',
                transparent=True,
                overlay=True,
                control=True,
                attr='SRTM30 via terrestris'
            ).add_to(map_obj)
        except Exception as e:
            self.logger.warning(f"Could not add SRTM30 WMS: {e}")

    def _add_buffer_tools(self, map_obj: folium.Map) -> None:
        """Add simple buffer tools (5km/10km) at last clicked point and a clear button."""
        buf_html = '''
        <div id="buffer-tools" style="position: fixed; top: 60px; right: 10px; z-index: 1000;">
          <div style="background: white; border: 2px solid #dee2e6; border-radius: 6px; padding: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); font-family: Arial; font-size: 12px;">
            <div style="font-weight: bold; margin-bottom: 6px;">Buffers rápidos</div>
            <button id="buf5" style="margin: 2px; padding: 6px;">5 km</button>
            <button id="buf10" style="margin: 2px; padding: 6px;">10 km</button>
            <button id="bufClear" style="margin: 2px; padding: 6px; background:#6c757d; color:white;">Limpiar</button>
          </div>
        </div>
        <script>
        (function(){
          var lastClick = null;
          var bufferLayers = [];
          var mapRef = ''' + map_obj.get_name() + ''';
          mapRef.on('click', function(e){ lastClick = e.latlng; });
          function addCircle(radiusMeters){
            if(!lastClick){ alert('Haz clic en el mapa para elegir el centro del buffer.'); return; }
            var circle = L.circle(lastClick, {radius: radiusMeters, color: '#007bff', weight: 2, fillOpacity: 0.05});
            circle.addTo(mapRef);
            bufferLayers.push(circle);
            mapRef.fitBounds(circle.getBounds(), {padding:[20,20]});
          }
          function clearBuffers(){
            bufferLayers.forEach(function(c){ try{ mapRef.removeLayer(c); }catch(e){} });
            bufferLayers = [];
          }
          document.getElementById('buf5').addEventListener('click', function(){ addCircle(5000); });
          document.getElementById('buf10').addEventListener('click', function(){ addCircle(10000); });
          document.getElementById('bufClear').addEventListener('click', function(){ clearBuffers(); });
        })();
        </script>
        '''
        map_obj.get_root().html.add_child(folium.Element(buf_html))
    
    def _add_custom_markers(self, map_obj: folium.Map) -> None:
        """
        Add custom markers for important locations in La Guajira.
        
        Args:
            map_obj: Folium map object
        """
        # Important locations in La Guajira, Colombia
        important_locations = [
            {
                'name': 'Riohacha',
                'lat': 11.5444,
                'lon': -72.9072,
                'type': 'city',
                'icon': 'building',
                'description': 'Capital de La Guajira'
            },
            {
                'name': 'Maicao',
                'lat': 11.3833,
                'lon': -72.2333,
                'type': 'city',
                'icon': 'home',
                'description': 'Ciudad fronteriza'
            },
            {
                'name': 'Uribia',
                'lat': 11.7167,
                'lon': -72.2667,
                'type': 'city',
                'icon': 'home',
                'description': 'Capital indígena'
            },
            {
                'name': 'Punta Gallinas',
                'lat': 12.4333,
                'lon': -71.6667,
                'type': 'landmark',
                'icon': 'map',
                'description': 'Punto más septentrional de Colombia'
            },
            {
                'name': 'Cabo de la Vela',
                'lat': 12.2167,
                'lon': -72.1167,
                'type': 'landmark',
                'icon': 'beach',
                'description': 'Destino turístico importante'
            },
            {
                'name': 'Parque Nacional Natural Macuira',
                'lat': 12.1667,
                'lon': -71.4167,
                'type': 'protected',
                'icon': '🌿',
                'description': 'Área protegida nacional'
            }
            ,
            {
                'name': 'Parque Eólico Jepírachi',
                'lat': 11.992,
                'lon': -72.277,
                'type': 'landmark',
                'icon': 'wind',
                'description': 'Parque eólico de referencia (demostrativo)'
            },
            {
                'name': 'Subestación Eléctrica (ref)',
                'lat': 11.527,
                'lon': -72.913,
                'type': 'protected',
                'icon': 'bolt',
                'description': 'Nodo de red eléctrica (referencia)'
            },
            {
                'name': 'Puerto de Riohacha',
                'lat': 11.541,
                'lon': -72.915,
                'type': 'landmark',
                'icon': 'ship',
                'description': 'Infraestructura portuaria'
            }
        ]
        
        # Create feature groups for different types of markers
        cities_group = folium.FeatureGroup(name='Ciudades')
        landmarks_group = folium.FeatureGroup(name='Puntos de Interés')
        protected_group = folium.FeatureGroup(name='Áreas Protegidas')
        
        for location in important_locations:
            # Create custom icon
            if location['type'] == 'city':
                icon_color = 'blue'
                icon_symbol = 'city'
            elif location['type'] == 'landmark':
                icon_color = 'orange'
                icon_symbol = 'star'
            elif location['type'] == 'protected':
                icon_color = 'green'
                icon_symbol = 'leaf'
            else:
                icon_color = 'red'
                icon_symbol = 'info'
            
            # Create popup content
            popup_content = f"""
            <div style="font-family: Arial, sans-serif; min-width: 200px;">
                <h4 style="margin: 0 0 10px 0; color: #333;">
                    {location['icon']} {location['name']}
                </h4>
                <p style="margin: 5px 0;"><strong>Tipo:</strong> {location['type'].title()}</p>
                <p style="margin: 5px 0;"><strong>Descripción:</strong> {location['description']}</p>
                <p style="margin: 5px 0;"><strong>Coordenadas:</strong> {location['lat']:.4f}, {location['lon']:.4f}</p>
                <p style="margin: 5px 0; font-size: 12px; color: #666;">
                    La Guajira, Colombia
                </p>
            </div>
            """
            
            # Create marker
            marker = folium.Marker(
                location=[location['lat'], location['lon']],
                popup=folium.Popup(popup_content, max_width=250),
                tooltip=location['name'],
                icon=folium.Icon(
                    color=icon_color,
                    icon=icon_symbol,
                    prefix='fa'
                )
            )
            
            # Add to appropriate group
            if location['type'] == 'city':
                marker.add_to(cities_group)
            elif location['type'] == 'landmark':
                marker.add_to(landmarks_group)
            elif location['type'] == 'protected':
                marker.add_to(protected_group)
        
        # Add groups to map
        cities_group.add_to(map_obj)
        landmarks_group.add_to(map_obj)
        protected_group.add_to(map_obj)
    
    def create_heatmap(self, 
                       wsi_data: np.ndarray,
                       aoi_geometry: Union[str, dict],
                       crs: str,
                       output_path: str,
                       title: str = "Índice de Idoneidad Eólica") -> None:
        """
        Create heatmap visualization of WSI data.
        
        Args:
            wsi_data: WSI raster data
            aoi_geometry: AOI geometry
            crs: Coordinate reference system
            output_path: Output HTML file path
            title: Map title
        """
        try:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Get map center
            center = self._get_geometry_center(aoi_geometry)
            
            # Create base map with multiple tile layers
            m = folium.Map(
                location=center,
                zoom_start=10,
                tiles=None  # We'll add tiles manually
            )
            
            # Add different tile layers
            self._add_tile_layers(m)
            
            # Get bounds and sample points
            bounds = self._get_geometry_bounds(aoi_geometry)
            sample_points = self._sample_wsi_points(wsi_data, bounds)
            
            # Prepare heatmap data
            heat_data = []
            for lat, lon, wsi_value in sample_points:
                heat_data.append([lat, lon, wsi_value])
            
            # Add heatmap layer
            plugins.HeatMap(
                heat_data,
                name='WSI Heatmap',
                min_opacity=0.4,
                max_zoom=18,
                radius=25,
                blur=15,
                gradient={0.0: 'red', 0.2: 'orange', 0.4: 'yellow', 0.6: 'lightgreen', 0.8: 'green', 1.0: 'darkgreen'}
            ).add_to(m)
            
            # Add AOI boundary
            self._add_aoi_boundary(m, aoi_geometry)
            
            # Add legend
            self._add_wsi_legend(m)
            
            # Add title
            self._add_title(m, title)
            
            # Add plugins
            self._add_plugins(m)
            
            # Save map
            m.save(output_path)
            
            self.logger.info(f"Heatmap saved to: {output_path}")
            
        except Exception as e:
            self.logger.error(f"Failed to create heatmap: {e}")
            raise
    
    def create_static_map(self, 
                         wsi_data: np.ndarray,
                         output_path: str,
                         title: str = "Wind Suitability Index") -> None:
        """
        Create static map visualization.
        
        Args:
            wsi_data: WSI raster data
            output_path: Output PNG file path
            title: Map title
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
            
            # Create figure
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Create colormap
            colors = ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen']
            n_bins = 5
            cmap = mcolors.LinearSegmentedColormap.from_list('wsi', colors, N=n_bins)
            
            # Plot WSI data
            im = ax.imshow(wsi_data, cmap=cmap, vmin=0, vmax=1)
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Wind Suitability Index', rotation=270, labelpad=20)
            
            # Set title
            ax.set_title(title, fontsize=16, fontweight='bold')
            
            # Remove axes
            ax.set_xticks([])
            ax.set_yticks([])
            
            # Save figure
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"Static map saved to: {output_path}")
            
        except Exception as e:
            self.logger.error(f"Failed to create static map: {e}")
            raise
    
    def add_admin_boundaries(self, 
                           map_obj: folium.Map, 
                           gdf: gpd.GeoDataFrame, 
                           name_field: str,
                           color: str = '#2E86AB') -> None:
        """
        Add administrative boundaries to map.
        
        Args:
            map_obj: Folium map object
            gdf: GeoDataFrame with administrative boundaries
            name_field: Field name for administrative unit names
            color: Color for boundary lines
        """
        try:
            # Create boundaries feature group
            boundaries_group = folium.FeatureGroup(name='Límites Administrativos')
            
            # Convert to WGS84 for Folium
            gdf_wgs84 = gdf.to_crs('EPSG:4326')
            
            # Add each boundary
            for idx, row in gdf_wgs84.iterrows():
                folium.GeoJson(
                    row.geometry.__geo_interface__,
                    style_function=lambda x: {
                        'fillColor': 'transparent',
                        'color': color,
                        'weight': 2,
                        'opacity': 0.8,
                        'fillOpacity': 0.1
                    },
                    popup=folium.Popup(
                        f"<b>{row[name_field]}</b><br>Unidad Administrativa", 
                        max_width=200
                    )
                ).add_to(boundaries_group)
            
            # Add boundaries group to map
            boundaries_group.add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add administrative boundaries: {e}")
            raise
    
    def add_choropleth_wsi(self, 
                         map_obj: folium.Map, 
                         gdf_stats: gpd.GeoDataFrame, 
                         column: str,
                         legend: str) -> None:
        """
        Add WSI choropleth layer to map.
        
        Args:
            map_obj: Folium map object
            gdf_stats: GeoDataFrame with WSI statistics
            column: Column name for WSI values
            legend: Legend title
        """
        try:
            # Convert to WGS84 for Folium
            gdf_wgs84 = gdf_stats.to_crs('EPSG:4326')
            
            # Create color scale
            min_val = gdf_wgs84[column].min()
            max_val = gdf_wgs84[column].max()
            
            # Define color scale (red to green)
            color_scale = ['#FF0000', '#FFA500', '#FFFF00', '#90EE90', '#006400']
            
            def get_color(value):
                if pd.isna(value):
                    return '#808080'  # Gray for missing values
                
                # Normalize value to 0-1
                normalized = (value - min_val) / (max_val - min_val) if max_val > min_val else 0
                
                # Get color index
                color_index = int(normalized * (len(color_scale) - 1))
                color_index = max(0, min(color_index, len(color_scale) - 1))
                
                return color_scale[color_index]
            
            # Add GeoJSON layer
            folium.GeoJson(
                gdf_wgs84.to_json(),
                style_function=lambda feature: {
                    'fillColor': get_color(feature['properties'].get(column)),
                    'color': 'white',
                    'weight': 1,
                    'fillOpacity': 0.7,
                    'opacity': 0.8
                },
                popup=folium.Popup(
                    self._create_wsi_popup_html(feature['properties']),
                    max_width=300
                ),
                tooltip=folium.Tooltip(
                    self._create_wsi_tooltip_text(feature['properties'], column),
                    permanent=False
                )
            ).add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add WSI choropleth: {e}")
            raise
    
    def _create_wsi_popup_html(self, properties: Dict) -> str:
        """Create HTML popup content for WSI data."""
        try:
            name = properties.get('name', 'N/A')
            mean_wsi = properties.get('mean_wsi', 0)
            p90_wsi = properties.get('p90_wsi', 0)
            pct_apto = properties.get('pct_apto', 0)
            valid_coverage = properties.get('valid_coverage', 0)
            
            html = f'''
            <div style="font-family: Arial, sans-serif;">
                <h4>{name}</h4>
                <p><b>WSI Promedio:</b> {mean_wsi:.3f}</p>
                <p><b>WSI P90:</b> {p90_wsi:.3f}</p>
                <p><b>% Área Apta:</b> {pct_apto:.1f}%</p>
                <p><b>Cobertura:</b> {valid_coverage:.1f}%</p>
            </div>
            '''
            
            return html
            
        except Exception as e:
            self.logger.error(f"Failed to create WSI popup HTML: {e}")
            return "<p>Error loading data</p>"
    
    def _create_wsi_tooltip_text(self, properties: Dict, column: str) -> str:
        """Create tooltip text for WSI data."""
        try:
            name = properties.get('name', 'N/A')
            value = properties.get(column, 0)
            
            return f"{name}: {column} {value:.3f}"
            
        except Exception as e:
            self.logger.error(f"Failed to create WSI tooltip text: {e}")
            return "Error loading data"
    
    def create_colombia_regions_map(self, config: Dict) -> str:
        """
        Create interactive map with Colombia region selector.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            Path to the generated map file
        """
        try:
            self.logger.info("Creating Colombia regions interactive map")
            
            # Load Colombia departments
            dept_file = "data/raw/admin/colombia_departamentos.geojson"
            if not Path(dept_file).exists():
                raise FileNotFoundError(f"Colombia departments file not found: {dept_file}")
            
            import geopandas as gpd
            dept_gdf = gpd.read_file(dept_file)
            
            # Convert to Web Mercator for consistency
            dept_gdf = dept_gdf.to_crs('EPSG:3857')
            
            # Create base map centered on Colombia
            center_lat = 4.0
            center_lon = -74.0
            map_obj = folium.Map(
                location=[center_lat, center_lon],
                zoom_start=6,
                tiles='OpenStreetMap'
            )
            
            # Add title
            self._add_title(map_obj, "Colombia - Selector de Regiones")
            
            # Add region selector
            self._add_colombia_region_selector(map_obj, dept_gdf)
            
            # Add departments layer
            self._add_colombia_departments_layer(map_obj, dept_gdf)
            
            # Add plugins
            self._add_plugins(map_obj)
            
            # Save map
            output_dir = Path(config.get('output_dir', 'outputs')) / 'maps'
            output_dir.mkdir(parents=True, exist_ok=True)
            
            map_file = output_dir / 'colombia_regions_map.html'
            map_obj.save(str(map_file))
            
            self.logger.info(f"Colombia regions map saved to: {map_file}")
            return str(map_file)
            
        except Exception as e:
            self.logger.error(f"Failed to create Colombia regions map: {e}")
            raise
    
    def _add_colombia_region_selector(self, map_obj: folium.Map, dept_gdf: gpd.GeoDataFrame) -> None:
        """Add Colombia region selector to map."""
        try:
            # Get unique regions
            regions = sorted(dept_gdf['REGION'].unique().tolist())
            
            # Create region selector HTML
            region_selector_html = f'''
            <div id="region-selector" style="position: fixed; 
                top: 100px; left: 50px; width: 250px; height: auto; 
                background-color: white; border: 2px solid #dee2e6; z-index: 9999; 
                font-size: 14px; padding: 15px; border-radius: 8px; 
                box-shadow: 0 2px 10px rgba(0,0,0,0.1)">
                <h4 style="margin: 0 0 15px 0; color: #333; text-align: center;">
                    🌎 Regiones de Colombia
                </h4>
                <div style="margin-bottom: 10px;">
                    <label for="region-select" style="display: block; margin-bottom: 5px; font-weight: bold;">
                        Seleccionar Región:
                    </label>
                    <select id="region-select" style="width: 100%; padding: 8px; border: 1px solid #ccc; border-radius: 4px;">
                        <option value="">-- Seleccionar Región --</option>
                        {''.join([f'<option value="{region}">{region}</option>' for region in regions])}
                    </select>
                </div>
                <div style="margin-bottom: 10px;">
                    <label for="dept-select" style="display: block; margin-bottom: 5px; font-weight: bold;">
                        Seleccionar Departamento:
                    </label>
                    <select id="dept-select" style="width: 100%; padding: 8px; border: 1px solid #ccc; border-radius: 4px;">
                        <option value="">-- Seleccionar Departamento --</option>
                    </select>
                </div>
                <div style="text-align: center;">
                    <button onclick="zoomToRegion()" style="
                        background: #007bff; color: white; border: none; 
                        padding: 8px 16px; border-radius: 4px; cursor: pointer; 
                        margin-right: 10px;">
                        🔍 Zoom a Región
                    </button>
                    <button onclick="resetView()" style="
                        background: #6c757d; color: white; border: none; 
                        padding: 8px 16px; border-radius: 4px; cursor: pointer;">
                        🏠 Vista General
                    </button>
                </div>
            </div>
            '''
            
            # Add JavaScript for region selector functionality
            region_selector_js = '''
            <script>
            // Department data
            const departments = ''' + dept_gdf.to_json() + ''';
            
            // Populate department dropdown based on region selection
            document.getElementById('region-select').addEventListener('change', function() {
                const region = this.value;
                const deptSelect = document.getElementById('dept-select');
                
                // Clear existing options
                deptSelect.innerHTML = '<option value="">-- Seleccionar Departamento --</option>';
                
                if (region) {
                    // Filter departments by region
                    const regionDepts = departments.features.filter(dept => 
                        dept.properties.REGION === region
                    );
                    
                    // Add department options
                    regionDepts.forEach(dept => {
                        const option = document.createElement('option');
                        option.value = dept.properties.DPTO_CNMBR;
                        option.textContent = dept.properties.DPTO_CNMBR;
                        deptSelect.appendChild(option);
                    });
                }
            });
            
            // Zoom to selected region
            function zoomToRegion() {
                const region = document.getElementById('region-select').value;
                const dept = document.getElementById('dept-select').value;
                
                if (dept) {
                    // Zoom to specific department
                    const deptData = departments.features.find(d => 
                        d.properties.DPTO_CNMBR === dept
                    );
                    
                    if (deptData) {
                        const coords = deptData.geometry.coordinates[0];
                        const lats = coords.map(coord => coord[1]);
                        const lons = coords.map(coord => coord[0]);
                        
                        const bounds = [
                            [Math.min(...lats), Math.min(...lons)],
                            [Math.max(...lats), Math.max(...lons)]
                        ];
                        
                        map.fitBounds(bounds, {padding: [20, 20]});
                        showNotification('Zoom a: ' + dept);
                    }
                } else if (region) {
                    // Zoom to region (all departments in region)
                    const regionDepts = departments.features.filter(d => 
                        d.properties.REGION === region
                    );
                    
                    if (regionDepts.length > 0) {
                        let allLats = [];
                        let allLons = [];
                        
                        regionDepts.forEach(dept => {
                            const coords = dept.geometry.coordinates[0];
                            allLats = allLats.concat(coords.map(coord => coord[1]));
                            allLons = allLons.concat(coords.map(coord => coord[0]));
                        });
                        
                        const bounds = [
                            [Math.min(...allLats), Math.min(...allLons)],
                            [Math.max(...allLats), Math.max(...allLons)]
                        ];
                        
                        map.fitBounds(bounds, {padding: [20, 20]});
                        showNotification('Zoom a región: ' + region);
                    }
                } else {
                    showNotification('Por favor selecciona una región o departamento');
                }
            }
            
            // Reset to general view
            function resetView() {
                map.setView([4.0, -74.0], 6);
                document.getElementById('region-select').value = '';
                document.getElementById('dept-select').innerHTML = '<option value="">-- Seleccionar Departamento --</option>';
                showNotification('Vista general de Colombia');
            }
            
            // Show notification
            function showNotification(message) {
                // Create notification element
                const notification = document.createElement('div');
                notification.style.cssText = `
                    position: fixed; top: 20px; right: 20px; z-index: 10000;
                    background: #28a745; color: white; padding: 10px 20px;
                    border-radius: 4px; box-shadow: 0 2px 10px rgba(0,0,0,0.2);
                    font-size: 14px; font-weight: bold;
                `;
                notification.textContent = message;
                document.body.appendChild(notification);
                
                // Remove after 3 seconds
                setTimeout(() => {
                    if (notification.parentNode) {
                        notification.parentNode.removeChild(notification);
                    }
                }, 3000);
            }
            </script>
            '''
            
            # Add HTML and JavaScript to map
            map_obj.get_root().html.add_child(folium.Element(region_selector_html))
            map_obj.get_root().html.add_child(folium.Element(region_selector_js))
            
        except Exception as e:
            self.logger.error(f"Failed to add Colombia region selector: {e}")
            raise
    
    def _add_colombia_departments_layer(self, map_obj: folium.Map, dept_gdf: gpd.GeoDataFrame) -> None:
        """Add Colombia departments layer to map."""
        try:
            # Convert to WGS84 for Folium
            dept_gdf_wgs84 = dept_gdf.to_crs('EPSG:4326')
            
            # Create departments feature group
            dept_group = folium.FeatureGroup(name='Departamentos de Colombia')
            
            # Define colors for different regions
            region_colors = {
                'Caribe': '#FF6B6B',
                'Andina': '#4ECDC4', 
                'Pacífico': '#45B7D1',
                'Amazonía': '#96CEB4',
                'Orinoquía': '#FFEAA7',
                'Insular': '#DDA0DD'
            }
            
            # Add each department
            for idx, row in dept_gdf_wgs84.iterrows():
                region = row['REGION']
                color = region_colors.get(region, '#95A5A6')
                
                folium.GeoJson(
                    row.geometry.__geo_interface__,
                    style_function=lambda x, region=region, color=color: {
                        'fillColor': color,
                        'color': 'white',
                        'weight': 1,
                        'fillOpacity': 0.6,
                        'opacity': 0.8
                    },
                    popup=folium.Popup(
                        f"<b>{row['DPTO_CNMBR']}</b><br>"
                        f"Región: {row['REGION']}<br>"
                        f"Código: {row['DPTO_CCDGO']}", 
                        max_width=200
                    ),
                    tooltip=folium.Tooltip(
                        f"{row['DPTO_CNMBR']} - {row['REGION']}",
                        permanent=False
                    )
                ).add_to(dept_group)
            
            # Add departments group to map
            dept_group.add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add Colombia departments layer: {e}")
            raise
    
    def create_ocha_wsi_map(self, config: Dict) -> str:
        """
        Create interactive WSI map with OCHA administrative boundaries and search.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            Path to the generated map file
        """
        try:
            self.logger.info("Creating OCHA WSI interactive map")
            
            # Load OCHA boundaries
            from .ocha_adapter import OCHAAdapter
            ocha_adapter = OCHAAdapter()
            
            admin_config = config.get('admin_boundaries', {})
            file_path = admin_config.get('file')
            crs = config.get('crs', 'EPSG:3857')
            
            if not file_path:
                raise ValueError("admin_boundaries.file is required in configuration")
            
            # Load OCHA boundaries
            admin_gdf = ocha_adapter.load_ocha_boundaries(file_path, crs)
            
            # Apply filters if specified
            filters = admin_config.get('filter', {})
            if filters:
                admin_gdf = ocha_adapter.filter_ocha_boundaries(admin_gdf, filters)
            
            # Load WSI statistics if available
            stats_file = f"outputs/reports/wsi_stats_by_{admin_config.get('target_level', 'municipio')}_ocha.json"
            stats_df = None
            
            if Path(stats_file).exists():
                try:
                    with open(stats_file, 'r', encoding='utf-8') as f:
                        stats_data = json.load(f)
                    stats_df = pd.DataFrame(stats_data['statistics'])
                    self.logger.info(f"Loaded WSI statistics for {len(stats_df)} units")
                except Exception as e:
                    self.logger.warning(f"Failed to load WSI statistics: {e}")
            
            # Create base map
            bounds = admin_gdf.total_bounds
            center_lat = (bounds[1] + bounds[3]) / 2
            center_lon = (bounds[0] + bounds[2]) / 2
            
            # Convert to WGS84 for Folium
            admin_gdf_wgs84 = admin_gdf.to_crs('EPSG:4326')
            
            map_obj = folium.Map(
                location=[center_lat, center_lon],
                zoom_start=8,
                tiles='OpenStreetMap'
            )
            
            # Add title
            title = f"WSI por {admin_config.get('target_level', 'municipio').title()} - OCHA COD-AB"
            self._add_title(map_obj, title)
            
            # Add OCHA boundaries with WSI data
            self._add_ocha_boundaries_layer(map_obj, admin_gdf_wgs84, stats_df, admin_config)
            
            # Add region selector
            self._add_ocha_region_selector(map_obj, admin_gdf_wgs84, admin_config)
            
            # Add search functionality
            self._add_ocha_search_plugin(map_obj, admin_gdf_wgs84, admin_config)
            
            # Add plugins
            self._add_plugins(map_obj)
            
            # Save map
            output_dir = Path(config.get('output_dir', 'outputs')) / 'maps'
            output_dir.mkdir(parents=True, exist_ok=True)
            
            map_file = output_dir / f'wsi_by_{admin_config.get("target_level", "municipio")}_ocha.html'
            map_obj.save(str(map_file))
            
            self.logger.info(f"OCHA WSI map saved to: {map_file}")
            return str(map_file)
            
        except Exception as e:
            self.logger.error(f"Failed to create OCHA WSI map: {e}")
            raise
    
    def _add_ocha_boundaries_layer(self, 
                                 map_obj: folium.Map, 
                                 admin_gdf: gpd.GeoDataFrame,
                                 stats_df: Optional[pd.DataFrame],
                                 admin_config: Dict) -> None:
        """Add OCHA administrative boundaries layer with WSI data."""
        try:
            # Create boundaries feature group
            boundaries_group = folium.FeatureGroup(name='Límites Administrativos OCHA')
            
            # Get field names
            id_field = admin_config.get('id_field', 'ADM2_PCODE')
            name_field = admin_config.get('name_field', 'ADM2_ES')
            parent_field = admin_config.get('parent_field', 'ADM1_ES')
            
            # Merge with statistics if available
            if stats_df is not None and not stats_df.empty:
                # Merge on ID field
                merged_gdf = admin_gdf.merge(stats_df, left_on=id_field, right_on=id_field, how='left')
                
                # Define color scale for WSI values
                if 'mean_wsi' in merged_gdf.columns:
                    min_wsi = merged_gdf['mean_wsi'].min()
                    max_wsi = merged_gdf['mean_wsi'].max()
                    
                    def get_wsi_color(value):
                        if pd.isna(value):
                            return '#808080'  # Gray for missing values
                        
                        # Normalize value to 0-1
                        normalized = (value - min_wsi) / (max_wsi - min_wsi) if max_wsi > min_wsi else 0
                        
                        # Color scale (red to green)
                        if normalized < 0.2:
                            return '#FF0000'  # Red
                        elif normalized < 0.4:
                            return '#FFA500'  # Orange
                        elif normalized < 0.6:
                            return '#FFFF00'  # Yellow
                        elif normalized < 0.8:
                            return '#90EE90'  # Light Green
                        else:
                            return '#006400'  # Dark Green
                else:
                    def get_wsi_color(value):
                        return '#95A5A6'  # Default gray
            else:
                merged_gdf = admin_gdf
                def get_wsi_color(value):
                    return '#95A5A6'  # Default gray
            
            # Add each administrative unit
            for idx, row in merged_gdf.iterrows():
                # Get WSI value for coloring
                wsi_value = row.get('mean_wsi', None) if stats_df is not None else None
                color = get_wsi_color(wsi_value)
                
                folium.GeoJson(
                    row.geometry.__geo_interface__,
                    style_function=lambda x, color=color: {
                        'fillColor': color,
                        'color': 'white',
                        'weight': 1,
                        'fillOpacity': 0.7,
                        'opacity': 0.8
                    },
                    popup=folium.Popup(
                        self._create_ocha_popup_html(row, admin_config, stats_df is not None),
                        max_width=300
                    ),
                    tooltip=folium.Tooltip(
                        self._create_ocha_tooltip_text(row, admin_config, stats_df is not None),
                        permanent=False
                    )
                ).add_to(boundaries_group)
            
            # Add boundaries group to map
            boundaries_group.add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add OCHA boundaries layer: {e}")
            raise
    
    def _add_ocha_region_selector(self, 
                                map_obj: folium.Map, 
                                admin_gdf: gpd.GeoDataFrame,
                                admin_config: Dict) -> None:
        """Add OCHA region selector to map."""
        try:
            # Get field names
            name_field = admin_config.get('name_field', 'ADM2_ES')
            parent_field = admin_config.get('parent_field', 'ADM1_ES')
            
            # Get unique departments and municipalities
            departments = sorted(admin_gdf[parent_field].unique().tolist())
            municipalities = sorted(admin_gdf[name_field].unique().tolist())
            
            # Create simple region selector HTML
            region_selector_html = f'''
            <div id="region-selector" style="position: fixed; 
                top: 100px; left: 50px; width: 250px; height: auto; 
                background-color: white; border: 2px solid #dee2e6; z-index: 9999; 
                font-size: 14px; padding: 15px; border-radius: 8px; 
                box-shadow: 0 2px 10px rgba(0,0,0,0.1)">
                <h4 style="margin: 0 0 15px 0; color: #333; text-align: center;">
                    🌎 Región OCHA
                </h4>
                <div style="margin-bottom: 10px;">
                    <label for="dept-select" style="display: block; margin-bottom: 5px; font-weight: bold;">
                        Departamento:
                    </label>
                    <select id="dept-select" style="width: 100%; padding: 8px; border: 1px solid #ccc; border-radius: 4px;">
                        <option value="">-- Seleccionar Departamento --</option>
                        {''.join([f'<option value="{dept}">{dept}</option>' for dept in departments])}
                    </select>
                </div>
                <div style="margin-bottom: 10px;">
                    <label for="municipio-select" style="display: block; margin-bottom: 5px; font-weight: bold;">
                        Municipio:
                    </label>
                    <select id="municipio-select" style="width: 100%; padding: 8px; border: 1px solid #ccc; border-radius: 4px;">
                        <option value="">-- Seleccionar Municipio --</option>
                        {''.join([f'<option value="{mun}">{mun}</option>' for mun in municipalities])}
                    </select>
                </div>
                <div style="text-align: center;">
                    <button onclick="zoomToRegion()" style="
                        background: #90EE90; color: #333; border: 1px solid #7ED321; 
                        padding: 6px 12px; border-radius: 4px; cursor: pointer; 
                        margin-right: 10px; font-size: 11px;">
                        Aplicar
                    </button>
                    <button onclick="resetView()" style="
                        background: #6c757d; color: white; border: none; 
                        padding: 6px 12px; border-radius: 4px; cursor: pointer; font-size: 11px;">
                        Reset
                    </button>
                </div>
            </div>
            '''
            
            # Create simple JavaScript and bind Folium map variable safely
            map_var = map_obj.get_name()
            region_selector_js = (
                '<script>\n'
                '// Bind Folium map instance to a predictable variable\n'
                f'var map = {map_var};\n'
            ) + r'''
            // Simple region data
            const regionData = {
                "La Guajira": {
                    "Riohacha": [[-73.0, 10.5], [-71.0, 10.5], [-71.0, 12.5], [-73.0, 12.5]],
                    "Maicao": [[-72.0, 11.0], [-71.0, 11.0], [-71.0, 12.0], [-72.0, 12.0]],
                    "Uribia": [[-72.5, 11.5], [-71.5, 11.5], [-71.5, 12.5], [-72.5, 12.5]],
                    "Manaure": [[-72.5, 11.0], [-71.5, 11.0], [-71.5, 12.0], [-72.5, 12.0]],
                    "San Juan del Cesar": [[-73.0, 10.5], [-72.0, 10.5], [-72.0, 11.5], [-73.0, 11.5]]
                }
            };
            
            // Filter municipalities based on department selection
            document.getElementById('dept-select').addEventListener('change', function() {
                const selectedDept = this.value;
                const municipioSelect = document.getElementById('municipio-select');
                
                // Clear existing options
                municipioSelect.innerHTML = '<option value="">-- Seleccionar Municipio --</option>';
                
                if (selectedDept && regionData[selectedDept]) {
                    // Add municipality options
                    Object.keys(regionData[selectedDept]).forEach(municipio => {
                        const option = document.createElement('option');
                        option.value = municipio;
                        option.textContent = municipio;
                        municipioSelect.appendChild(option);
                    });
                }
            });
            
            // Zoom to selected region
            function zoomToRegion() {
                const selectedDept = document.getElementById('dept-select').value;
                const selectedMun = document.getElementById('municipio-select').value;
                
                if (selectedMun && regionData[selectedDept] && regionData[selectedDept][selectedMun]) {
                    // Zoom to specific municipality
                    const coords = regionData[selectedDept][selectedMun];
                    const lats = coords.map(coord => coord[1]);
                    const lons = coords.map(coord => coord[0]);
                    
                    const bounds = [
                        [Math.min(...lats), Math.min(...lons)],
                        [Math.max(...lats), Math.max(...lons)]
                    ];
                    
                    map.fitBounds(bounds, {padding: [20, 20]});
                    showNotification('Zoom a: ' + selectedMun);
                } else if (selectedDept && regionData[selectedDept]) {
                    // Zoom to department (all municipalities in department)
                    const deptData = regionData[selectedDept];
                    let allLats = [];
                    let allLons = [];
                    
                    Object.values(deptData).forEach(coords => {
                        allLats = allLats.concat(coords.map(coord => coord[1]));
                        allLons = allLons.concat(coords.map(coord => coord[0]));
                    });
                    
                    const bounds = [
                        [Math.min(...allLats), Math.min(...allLons)],
                        [Math.max(...allLats), Math.max(...allLons)]
                    ];
                    
                    map.fitBounds(bounds, {padding: [20, 20]});
                    showNotification('Zoom a departamento: ' + selectedDept);
                } else {
                    showNotification('Por favor selecciona un departamento o municipio');
                }
            }
            
            // Reset to general view
            function resetView() {
                map.setView([11.5, -72.0], 8);
                document.getElementById('dept-select').value = '';
                document.getElementById('municipio-select').innerHTML = '<option value="">-- Seleccionar Municipio --</option>';
                showNotification('Vista general');
            }
            
            // Show notification
            function showNotification(message) {
                // Create notification element
                const notification = document.createElement('div');
                notification.style.cssText = `
                    position: fixed; top: 20px; right: 20px; z-index: 10000;
                    background: #28a745; color: white; padding: 10px 20px;
                    border-radius: 4px; box-shadow: 0 2px 10px rgba(0,0,0,0.2);
                    font-size: 14px; font-weight: bold;
                `;
                notification.textContent = message;
                document.body.appendChild(notification);
                
                // Remove after 3 seconds
                setTimeout(() => {
                    if (notification.parentNode) {
                        notification.parentNode.removeChild(notification);
                    }
                }, 3000);
            }
            </script>
            '''
            
            # Add HTML and JavaScript to map
            map_obj.get_root().html.add_child(folium.Element(region_selector_html))
            map_obj.get_root().html.add_child(folium.Element(region_selector_js))
            
        except Exception as e:
            self.logger.error(f"Failed to add OCHA region selector: {e}")
            raise
    
    def _add_ocha_search_plugin(self, 
                              map_obj: folium.Map, 
                              admin_gdf: gpd.GeoDataFrame,
                              admin_config: Dict) -> None:
        """Add search plugin for OCHA administrative units."""
        try:
            # Get field names
            name_field = admin_config.get('name_field', 'ADM2_ES')
            parent_field = admin_config.get('parent_field', 'ADM1_ES')
            
            # Create search data
            search_data = []
            for idx, row in admin_gdf.iterrows():
                search_data.append({
                    'name': f"{row[name_field]} - {row[parent_field]}",
                    'geometry': row.geometry.__geo_interface__
                })
            
            # Add search plugin
            search_plugin = folium.plugins.Search(
                layer=folium.FeatureGroup(name='Search Results'),
                geom_type='Polygon',
                placeholder='Buscar municipio...',
                collapsed=False,
                search_label=name_field,
                search_zoom=12
            )
            
            # Add search data
            search_layer = folium.FeatureGroup(name='Search Results')
            for item in search_data:
                folium.GeoJson(
                    item['geometry'],
                    style_function=lambda x: {
                        'fillColor': '#FFD700',
                        'color': '#FF8C00',
                        'weight': 2,
                        'fillOpacity': 0.3,
                        'opacity': 0.8
                    }
                ).add_to(search_layer)
            
            search_layer.add_to(map_obj)
            search_plugin.add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add OCHA search plugin: {e}")
            raise
    
    def _create_ocha_popup_html(self, row: pd.Series, admin_config: Dict, has_stats: bool) -> str:
        """Create HTML popup content for OCHA data."""
        try:
            name_field = admin_config.get('name_field', 'ADM2_ES')
            parent_field = admin_config.get('parent_field', 'ADM1_ES')
            id_field = admin_config.get('id_field', 'ADM2_PCODE')
            
            name = row.get(name_field, 'N/A')
            parent = row.get(parent_field, 'N/A')
            unit_id = row.get(id_field, 'N/A')
            
            html = f'''
            <div style="font-family: Arial, sans-serif;">
                <h4>{name}</h4>
                <p><b>Departamento:</b> {parent}</p>
                <p><b>Código P-CODE:</b> {unit_id}</p>
            '''
            
            if has_stats:
                mean_wsi = row.get('mean_wsi', 0)
                p90_wsi = row.get('p90_wsi', 0)
                pct_apto = row.get('pct_apto', 0)
                valid_coverage = row.get('valid_coverage', 0)
                
                html += f'''
                <hr>
                <h5>Estadísticas WSI:</h5>
                <p><b>WSI Promedio:</b> {mean_wsi:.3f}</p>
                <p><b>WSI P90:</b> {p90_wsi:.3f}</p>
                <p><b>% Área Apta:</b> {pct_apto:.1f}%</p>
                <p><b>Cobertura:</b> {valid_coverage:.1f}%</p>
                '''
            
            html += '</div>'
            return html
            
        except Exception as e:
            self.logger.error(f"Failed to create OCHA popup HTML: {e}")
            return "<p>Error loading data</p>"
    
    def _create_ocha_tooltip_text(self, row: pd.Series, admin_config: Dict, has_stats: bool) -> str:
        """Create tooltip text for OCHA data."""
        try:
            name_field = admin_config.get('name_field', 'ADM2_ES')
            parent_field = admin_config.get('parent_field', 'ADM1_ES')
            
            name = row.get(name_field, 'N/A')
            parent = row.get(parent_field, 'N/A')
            
            if has_stats and 'mean_wsi' in row:
                wsi_value = row.get('mean_wsi', 0)
                return f"{name} - {parent} (WSI: {wsi_value:.3f})"
            else:
                return f"{name} - {parent}"
            
        except Exception as e:
            self.logger.error(f"Failed to create OCHA tooltip text: {e}")
            return "Error loading data"



