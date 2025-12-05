"""
Map WSI by administrative boundaries use case.

This module handles the generation of interactive choropleth maps
showing WSI statistics by administrative units.
"""

import folium
import pandas as pd
import geopandas as gpd
from typing import Dict, List, Optional, Tuple, Union
import logging
from pathlib import Path
import json

from ..infrastructure.admin_boundaries import AdminBoundariesAdapter
from ..infrastructure.folium_map import FoliumMapGenerator

logger = logging.getLogger(__name__)


class MapWSIByAdminUseCase:
    """Use case for generating WSI choropleth maps by administrative boundaries."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.admin_adapter = AdminBoundariesAdapter()
        self.folium_adapter = FoliumMapGenerator()
    
    def execute(self, config: Dict) -> Dict:
        """
        Execute WSI choropleth map generation by administrative boundaries.
        
        Args:
            config: Configuration dictionary with admin_boundaries settings
            
        Returns:
            Dictionary with execution results
        """
        try:
            self.logger.info("Starting WSI choropleth map generation by administrative boundaries")
            
            # Extract configuration
            admin_config = config.get('admin_boundaries', {})
            output_dir = config.get('output_dir', 'outputs')
            target_level = admin_config.get('target_level', 'admin')
            crs = config.get('crs', 'EPSG:3857')
            
            # Load administrative boundaries
            admin_gdf = self._load_admin_boundaries(admin_config, crs)
            
            # Load statistics
            stats_df = self._load_statistics(output_dir, target_level)
            
            # Merge boundaries with statistics
            merged_gdf = self._merge_boundaries_with_stats(admin_gdf, stats_df, admin_config)
            
            # Generate choropleth map
            map_file = self._generate_choropleth_map(
                merged_gdf, 
                admin_config, 
                output_dir,
                target_level
            )
            
            self.logger.info("WSI choropleth map generation completed successfully")
            
            return {
                "success": True,
                "map_file": map_file,
                "total_units": len(merged_gdf),
                "statistics_loaded": len(stats_df)
            }
            
        except Exception as e:
            self.logger.error(f"Failed to generate WSI choropleth map: {e}")
            raise
    
    def _load_admin_boundaries(self, admin_config: Dict, crs: str) -> gpd.GeoDataFrame:
        """Load administrative boundaries."""
        try:
            file_path = admin_config.get('file')
            layer = admin_config.get('layer')
            
            if not file_path:
                raise ValueError("admin_boundaries.file is required in configuration")
            
            # Load boundaries
            admin_gdf = self.admin_adapter.load_admin_boundaries(file_path, layer, crs)
            
            # Apply filters if specified
            filters = admin_config.get('filter', {})
            if filters:
                admin_gdf = self.admin_adapter.filter_admin_boundaries(admin_gdf, filters)
            
            return admin_gdf
            
        except Exception as e:
            self.logger.error(f"Failed to load administrative boundaries: {e}")
            raise
    
    def _load_statistics(self, output_dir: str, target_level: str) -> pd.DataFrame:
        """Load WSI statistics from JSON file."""
        try:
            stats_file = Path(output_dir) / 'reports' / f'wsi_stats_by_{target_level}.json'
            
            if not stats_file.exists():
                raise FileNotFoundError(f"Statistics file not found: {stats_file}")
            
            with open(stats_file, 'r', encoding='utf-8') as f:
                stats_data = json.load(f)
            
            # Convert to DataFrame
            stats_df = pd.DataFrame(stats_data['statistics'])
            
            self.logger.info(f"Loaded statistics for {len(stats_df)} administrative units")
            return stats_df
            
        except Exception as e:
            self.logger.error(f"Failed to load statistics: {e}")
            raise
    
    def _merge_boundaries_with_stats(self, 
                                   admin_gdf: gpd.GeoDataFrame,
                                   stats_df: pd.DataFrame,
                                   admin_config: Dict) -> gpd.GeoDataFrame:
        """Merge administrative boundaries with statistics."""
        try:
            id_field = admin_config.get('id_field', 'id')
            
            # Merge on ID field
            merged_gdf = admin_gdf.merge(stats_df, on=id_field, how='left')
            
            # Check for missing statistics
            missing_stats = merged_gdf['mean_wsi'].isna().sum()
            if missing_stats > 0:
                self.logger.warning(f"Missing statistics for {missing_stats} administrative units")
            
            self.logger.info(f"Merged boundaries with statistics: {len(merged_gdf)} units")
            return merged_gdf
            
        except Exception as e:
            self.logger.error(f"Failed to merge boundaries with statistics: {e}")
            raise
    
    def _generate_choropleth_map(self, 
                               merged_gdf: gpd.GeoDataFrame,
                               admin_config: Dict,
                               output_dir: str,
                               target_level: str) -> str:
        """Generate choropleth map with WSI statistics."""
        try:
            # Create map centered on the data
            bounds = merged_gdf.total_bounds
            center_lat = (bounds[1] + bounds[3]) / 2
            center_lon = (bounds[0] + bounds[2]) / 2
            
            # Convert to WGS84 for Folium
            gdf_wgs84 = merged_gdf.to_crs('EPSG:4326')
            
            # Create base map
            map_obj = folium.Map(
                location=[center_lat, center_lon],
                zoom_start=8,
                tiles='OpenStreetMap'
            )
            
            # Add title
            title = f"Índice de Idoneidad Eólica por {target_level.title()}"
            self.folium_adapter._add_title(map_obj, title)
            
            # Add choropleth layer
            self._add_choropleth_layer(
                map_obj, 
                gdf_wgs84, 
                admin_config,
                'mean_wsi',
                'WSI Promedio'
            )
            
            # Add administrative boundaries
            self._add_admin_boundaries_layer(
                map_obj, 
                gdf_wgs84, 
                admin_config
            )
            
            # Add legend
            self._add_choropleth_legend(map_obj, 'mean_wsi', 'WSI Promedio')
            
            # Add plugins
            self.folium_adapter._add_plugins(map_obj)
            
            # Save map
            maps_dir = Path(output_dir) / 'maps'
            maps_dir.mkdir(parents=True, exist_ok=True)
            
            map_file = maps_dir / f'wsi_by_{target_level}.html'
            map_obj.save(str(map_file))
            
            self.logger.info(f"Choropleth map saved to {map_file}")
            return str(map_file)
            
        except Exception as e:
            self.logger.error(f"Failed to generate choropleth map: {e}")
            raise
    
    def _add_choropleth_layer(self, 
                            map_obj: folium.Map,
                            gdf: gpd.GeoDataFrame,
                            admin_config: Dict,
                            column: str,
                            name: str) -> None:
        """Add choropleth layer to map."""
        try:
            # Create color scale
            min_val = gdf[column].min()
            max_val = gdf[column].max()
            
            # Define color scale (green to red)
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
                gdf.to_json(),
                style_function=lambda feature: {
                    'fillColor': get_color(feature['properties'].get(column)),
                    'color': 'white',
                    'weight': 1,
                    'fillOpacity': 0.7,
                    'opacity': 0.8
                },
                popup=folium.Popup(
                    self._create_popup_html(feature['properties'], admin_config),
                    max_width=300
                ),
                tooltip=folium.Tooltip(
                    self._create_tooltip_text(feature['properties'], admin_config),
                    permanent=False
                )
            ).add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add choropleth layer: {e}")
            raise
    
    def _add_admin_boundaries_layer(self, 
                                  map_obj: folium.Map,
                                  gdf: gpd.GeoDataFrame,
                                  admin_config: Dict) -> None:
        """Add administrative boundaries layer."""
        try:
            # Add boundaries as separate layer
            boundaries_group = folium.FeatureGroup(name='Límites Administrativos')
            
            folium.GeoJson(
                gdf.to_json(),
                style_function=lambda x: {
                    'fillColor': 'transparent',
                    'color': '#2E86AB',
                    'weight': 2,
                    'opacity': 0.8,
                    'fillOpacity': 0.1
                }
            ).add_to(boundaries_group)
            
            boundaries_group.add_to(map_obj)
            
        except Exception as e:
            self.logger.error(f"Failed to add administrative boundaries layer: {e}")
            raise
    
    def _add_choropleth_legend(self, map_obj: folium.Map, column: str, title: str) -> None:
        """Add choropleth legend to map."""
        try:
            legend_html = f'''
            <div style="position: fixed; 
                        bottom: 50px; right: 50px; width: 200px; height: 120px; 
                        background-color: white; border: 2px solid grey; z-index:9999; 
                        font-size:14px; padding: 10px">
            <p><b>{title}</b></p>
            <p><i class="fa fa-square" style="color:#FF0000"></i> Muy Bajo (0.0-0.2)</p>
            <p><i class="fa fa-square" style="color:#FFA500"></i> Bajo (0.2-0.4)</p>
            <p><i class="fa fa-square" style="color:#FFFF00"></i> Regular (0.4-0.6)</p>
            <p><i class="fa fa-square" style="color:#90EE90"></i> Alto (0.6-0.8)</p>
            <p><i class="fa fa-square" style="color:#006400"></i> Muy Alto (0.8-1.0)</p>
            </div>
            '''
            
            map_obj.get_root().html.add_child(folium.Element(legend_html))
            
        except Exception as e:
            self.logger.error(f"Failed to add choropleth legend: {e}")
            raise
    
    def _create_popup_html(self, properties: Dict, admin_config: Dict) -> str:
        """Create HTML popup content."""
        try:
            name_field = admin_config.get('name_field', 'name')
            id_field = admin_config.get('id_field', 'id')
            
            name = properties.get(name_field, 'N/A')
            unit_id = properties.get(id_field, 'N/A')
            mean_wsi = properties.get('mean_wsi', 0)
            p90_wsi = properties.get('p90_wsi', 0)
            pct_apto = properties.get('pct_apto', 0)
            valid_coverage = properties.get('valid_coverage', 0)
            
            html = f'''
            <div style="font-family: Arial, sans-serif;">
                <h4>{name}</h4>
                <p><b>ID:</b> {unit_id}</p>
                <p><b>WSI Promedio:</b> {mean_wsi:.3f}</p>
                <p><b>WSI P90:</b> {p90_wsi:.3f}</p>
                <p><b>% Área Apta:</b> {pct_apto:.1f}%</p>
                <p><b>Cobertura:</b> {valid_coverage:.1f}%</p>
            </div>
            '''
            
            return html
            
        except Exception as e:
            self.logger.error(f"Failed to create popup HTML: {e}")
            return "<p>Error loading data</p>"
    
    def _create_tooltip_text(self, properties: Dict, admin_config: Dict) -> str:
        """Create tooltip text."""
        try:
            name_field = admin_config.get('name_field', 'name')
            name = properties.get(name_field, 'N/A')
            mean_wsi = properties.get('mean_wsi', 0)
            
            return f"{name}: WSI {mean_wsi:.3f}"
            
        except Exception as e:
            self.logger.error(f"Failed to create tooltip text: {e}")
            return "Error loading data"
