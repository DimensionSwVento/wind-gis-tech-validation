"""
Administrative boundaries infrastructure module.

This module handles loading and processing of administrative boundaries
(departments and municipalities) for WSI analysis by administrative units.
"""

import geopandas as gpd
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
import logging
from pathlib import Path
import json

logger = logging.getLogger(__name__)


class AdminBoundariesAdapter:
    """Adapter for handling administrative boundaries data."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def load_admin_boundaries(self, 
                            file_path: str, 
                            layer: Optional[str] = None,
                            crs: str = "EPSG:3857") -> gpd.GeoDataFrame:
        """
        Load administrative boundaries from vector file.
        
        Args:
            file_path: Path to the vector file (.gpkg, .geojson, .shp)
            layer: Layer name for multi-layer formats (e.g., .gpkg)
            crs: Target CRS for reprojection
            
        Returns:
            GeoDataFrame with administrative boundaries
        """
        try:
            self.logger.info(f"Loading administrative boundaries from {file_path}")
            
            # Load the vector file
            if layer:
                gdf = gpd.read_file(file_path, layer=layer)
            else:
                gdf = gpd.read_file(file_path)
            
            # Ensure it's a GeoDataFrame
            if not isinstance(gdf, gpd.GeoDataFrame):
                raise ValueError("Loaded data is not a GeoDataFrame")
            
            # Check if geometry column exists
            if not hasattr(gdf, 'geometry'):
                raise ValueError("No geometry column found in the data")
            
            # Reproject to target CRS
            if gdf.crs != crs:
                self.logger.info(f"Reprojecting from {gdf.crs} to {crs}")
                gdf = gdf.to_crs(crs)
            
            # Repair geometries
            gdf = self._repair_geometries(gdf)
            
            self.logger.info(f"Loaded {len(gdf)} administrative units")
            return gdf
            
        except Exception as e:
            self.logger.error(f"Failed to load administrative boundaries: {e}")
            raise
    
    def _repair_geometries(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Repair invalid geometries using buffer(0).
        
        Args:
            gdf: GeoDataFrame with potentially invalid geometries
            
        Returns:
            GeoDataFrame with repaired geometries
        """
        try:
            # Check for invalid geometries
            invalid_mask = ~gdf.geometry.is_valid
            invalid_count = invalid_mask.sum()
            
            if invalid_count > 0:
                self.logger.warning(f"Found {invalid_count} invalid geometries, repairing...")
                
                # Repair invalid geometries
                gdf.loc[invalid_mask, 'geometry'] = gdf.loc[invalid_mask, 'geometry'].buffer(0)
                
                # Verify repair
                still_invalid = ~gdf.geometry.is_valid
                if still_invalid.sum() > 0:
                    self.logger.warning(f"Could not repair {still_invalid.sum()} geometries, removing them")
                    gdf = gdf[~still_invalid]
            
            return gdf
            
        except Exception as e:
            self.logger.error(f"Failed to repair geometries: {e}")
            raise
    
    def filter_admin_boundaries(self, 
                              gdf: gpd.GeoDataFrame,
                              filters: Dict[str, List[str]]) -> gpd.GeoDataFrame:
        """
        Filter administrative boundaries based on criteria.
        
        Args:
            gdf: GeoDataFrame with administrative boundaries
            filters: Dictionary with field names as keys and lists of values to filter by
            
        Returns:
            Filtered GeoDataFrame
        """
        try:
            if not filters:
                return gdf
            
            self.logger.info(f"Applying filters: {filters}")
            
            filtered_gdf = gdf.copy()
            
            for field, values in filters.items():
                if field in filtered_gdf.columns:
                    # Filter by field values
                    mask = filtered_gdf[field].isin(values)
                    filtered_gdf = filtered_gdf[mask]
                    self.logger.info(f"Filtered by {field}: {mask.sum()} units remaining")
                else:
                    self.logger.warning(f"Field {field} not found in data, skipping filter")
            
            self.logger.info(f"Final filtered count: {len(filtered_gdf)} units")
            return filtered_gdf
            
        except Exception as e:
            self.logger.error(f"Failed to filter administrative boundaries: {e}")
            raise
    
    def validate_admin_data(self, 
                          gdf: gpd.GeoDataFrame,
                          required_fields: List[str]) -> bool:
        """
        Validate administrative boundaries data.
        
        Args:
            gdf: GeoDataFrame to validate
            required_fields: List of required field names
            
        Returns:
            True if valid, False otherwise
        """
        try:
            # Check if GeoDataFrame is not empty
            if len(gdf) == 0:
                self.logger.error("Administrative boundaries data is empty")
                return False
            
            # Check required fields
            missing_fields = [field for field in required_fields if field not in gdf.columns]
            if missing_fields:
                self.logger.error(f"Missing required fields: {missing_fields}")
                return False
            
            # Check for valid geometries
            invalid_geometries = ~gdf.geometry.is_valid
            if invalid_geometries.sum() > 0:
                self.logger.error(f"Found {invalid_geometries.sum()} invalid geometries")
                return False
            
            # Check for empty geometries
            empty_geometries = gdf.geometry.is_empty
            if empty_geometries.sum() > 0:
                self.logger.error(f"Found {empty_geometries.sum()} empty geometries")
                return False
            
            self.logger.info("Administrative boundaries data validation passed")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to validate administrative boundaries: {e}")
            return False
    
    def get_admin_bounds(self, gdf: gpd.GeoDataFrame) -> Tuple[float, float, float, float]:
        """
        Get bounding box of administrative boundaries.
        
        Args:
            gdf: GeoDataFrame with administrative boundaries
            
        Returns:
            Tuple of (minx, miny, maxx, maxy)
        """
        try:
            bounds = gdf.total_bounds
            return tuple(bounds)
        except Exception as e:
            self.logger.error(f"Failed to get administrative boundaries bounds: {e}")
            raise
    
    def save_admin_stats(self, 
                        stats_df: pd.DataFrame, 
                        output_path: str,
                        level: str) -> None:
        """
        Save administrative statistics to JSON file.
        
        Args:
            stats_df: DataFrame with statistics
            output_path: Output directory path
            level: Administrative level (municipio, departamento)
        """
        try:
            output_file = Path(output_path) / f"wsi_stats_by_{level}.json"
            
            # Convert DataFrame to dictionary
            stats_dict = {
                "metadata": {
                    "level": level,
                    "total_units": len(stats_df),
                    "generated_at": pd.Timestamp.now().isoformat()
                },
                "statistics": stats_df.to_dict('records')
            }
            
            # Save to JSON
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(stats_dict, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"Saved administrative statistics to {output_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to save administrative statistics: {e}")
            raise

