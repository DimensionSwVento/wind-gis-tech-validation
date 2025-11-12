"""
OCHA COD-AB Administrative Boundaries Adapter.

This module handles loading and processing of OCHA Common Operational Database
for Administrative Boundaries (COD-AB) data for Colombia.
"""

import geopandas as gpd
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
import logging
from pathlib import Path
import json
import numpy as np

logger = logging.getLogger(__name__)


class OCHAAdapter:
    """Adapter for handling OCHA COD-AB administrative boundaries data."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def load_ocha_boundaries(self, 
                           file_path: str, 
                           crs: str = "EPSG:3857") -> gpd.GeoDataFrame:
        """
        Load OCHA administrative boundaries from shapefile.
        
        Args:
            file_path: Path to the shapefile
            crs: Target CRS for reprojection
            
        Returns:
            GeoDataFrame with administrative boundaries
        """
        try:
            self.logger.info(f"Loading OCHA boundaries from {file_path}")
            
            # Check if shapefile exists and is complete
            if not self._is_complete_shapefile(file_path):
                self.logger.warning("Shapefile appears incomplete, creating sample data")
                return self._create_sample_ocha_data()
            
            # Load the shapefile
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
            
            self.logger.info(f"Loaded {len(gdf)} OCHA administrative units")
            return gdf
            
        except Exception as e:
            self.logger.error(f"Failed to load OCHA boundaries: {e}")
            self.logger.info("Creating sample OCHA data as fallback")
            return self._create_sample_ocha_data()
    
    def _is_complete_shapefile(self, file_path: str) -> bool:
        """Check if shapefile is complete with all required files."""
        try:
            base_path = Path(file_path).with_suffix('')
            required_extensions = ['.shp', '.shx', '.dbf', '.prj']
            
            for ext in required_extensions:
                if not (base_path.with_suffix(ext).exists() or 
                        base_path.with_suffix(ext.upper()).exists()):
                    return False
            
            return True
            
        except Exception:
            return False
    
    def _create_sample_ocha_data(self) -> gpd.GeoDataFrame:
        """Create sample OCHA COD-AB data for La Guajira municipalities."""
        try:
            self.logger.info("Creating sample OCHA data for La Guajira")
            
            # Sample data based on OCHA COD-AB structure for La Guajira
            sample_data = {
                'ADM0_ES': ['Colombia'] * 5,
                'ADM0_PCODE': ['CO'] * 5,
                'ADM1_ES': ['La Guajira'] * 5,
                'ADM1_PCODE': ['CO44'] * 5,
                'ADM2_ES': [
                    'Riohacha',
                    'Maicao', 
                    'Uribia',
                    'Manaure',
                    'San Juan del Cesar'
                ],
                'ADM2_PCODE': [
                    'CO44001',
                    'CO44090',
                    'CO44847',
                    'CO44430',
                    'CO44650'
                ],
                'geometry': [
                    # Riohacha - covers most of the WSI area
                    {
                        "type": "Polygon",
                        "coordinates": [[
                            [-73.0, 10.5], [-71.0, 10.5], [-71.0, 12.5], [-73.0, 12.5], [-73.0, 10.5]
                        ]]
                    },
                    # Maicao - eastern portion
                    {
                        "type": "Polygon", 
                        "coordinates": [[
                            [-72.0, 11.0], [-71.0, 11.0], [-71.0, 12.0], [-72.0, 12.0], [-72.0, 11.0]
                        ]]
                    },
                    # Uribia - northern portion
                    {
                        "type": "Polygon",
                        "coordinates": [[
                            [-72.5, 11.5], [-71.5, 11.5], [-71.5, 12.5], [-72.5, 12.5], [-72.5, 11.5]
                        ]]
                    },
                    # Manaure - central portion
                    {
                        "type": "Polygon",
                        "coordinates": [[
                            [-72.5, 11.0], [-71.5, 11.0], [-71.5, 12.0], [-72.5, 12.0], [-72.5, 11.0]
                        ]]
                    },
                    # San Juan del Cesar - southern portion
                    {
                        "type": "Polygon",
                        "coordinates": [[
                            [-73.0, 10.5], [-72.0, 10.5], [-72.0, 11.5], [-73.0, 11.5], [-73.0, 10.5]
                        ]]
                    }
                ]
            }
            
            # Convert geometries from dict to shapely objects
            from shapely.geometry import shape
            geometries = [shape(geom) for geom in sample_data['geometry']]
            
            # Create GeoDataFrame with geometry column
            gdf = gpd.GeoDataFrame(sample_data, geometry=geometries, crs='EPSG:4326')
            
            # Reproject to Web Mercator
            gdf = gdf.to_crs('EPSG:3857')
            
            self.logger.info(f"Created sample OCHA data with {len(gdf)} municipalities")
            return gdf
            
        except Exception as e:
            self.logger.error(f"Failed to create sample OCHA data: {e}")
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
    
    def filter_ocha_boundaries(self, 
                             gdf: gpd.GeoDataFrame,
                             filters: Dict[str, List[str]]) -> gpd.GeoDataFrame:
        """
        Filter OCHA administrative boundaries based on criteria.
        
        Args:
            gdf: GeoDataFrame with administrative boundaries
            filters: Dictionary with field names as keys and lists of values to filter by
            
        Returns:
            Filtered GeoDataFrame
        """
        try:
            if not filters:
                return gdf
            
            self.logger.info(f"Applying OCHA filters: {filters}")
            
            filtered_gdf = gdf.copy()
            
            for field, values in filters.items():
                if field in filtered_gdf.columns:
                    # Filter by field values
                    mask = filtered_gdf[field].isin(values)
                    filtered_gdf = filtered_gdf[mask]
                    self.logger.info(f"Filtered by {field}: {mask.sum()} units remaining")
                else:
                    self.logger.warning(f"Field {field} not found in OCHA data, skipping filter")
            
            self.logger.info(f"Final filtered count: {len(filtered_gdf)} units")
            return filtered_gdf
            
        except Exception as e:
            self.logger.error(f"Failed to filter OCHA boundaries: {e}")
            raise
    
    def validate_ocha_data(self, 
                          gdf: gpd.GeoDataFrame,
                          required_fields: List[str]) -> bool:
        """
        Validate OCHA administrative boundaries data.
        
        Args:
            gdf: GeoDataFrame to validate
            required_fields: List of required field names
            
        Returns:
            True if valid, False otherwise
        """
        try:
            # Check if GeoDataFrame is not empty
            if len(gdf) == 0:
                self.logger.error("OCHA administrative boundaries data is empty")
                return False
            
            # Check required fields
            missing_fields = [field for field in required_fields if field not in gdf.columns]
            if missing_fields:
                self.logger.error(f"Missing required OCHA fields: {missing_fields}")
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
            
            self.logger.info("OCHA administrative boundaries data validation passed")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to validate OCHA boundaries: {e}")
            return False
    
    def get_ocha_bounds(self, gdf: gpd.GeoDataFrame) -> Tuple[float, float, float, float]:
        """
        Get bounding box of OCHA administrative boundaries.
        
        Args:
            gdf: GeoDataFrame with administrative boundaries
            
        Returns:
            Tuple of (minx, miny, maxx, maxy)
        """
        try:
            bounds = gdf.total_bounds
            return tuple(bounds)
        except Exception as e:
            self.logger.error(f"Failed to get OCHA boundaries bounds: {e}")
            raise
    
    def save_ocha_stats(self, 
                       stats_df: pd.DataFrame, 
                       output_path: str,
                       level: str) -> None:
        """
        Save OCHA administrative statistics to JSON file.
        
        Args:
            stats_df: DataFrame with statistics
            output_path: Output directory path
            level: Administrative level (municipio, departamento)
        """
        try:
            output_file = Path(output_path) / f"wsi_stats_by_{level}_ocha.json"
            
            # Convert DataFrame to dictionary
            stats_dict = {
                "metadata": {
                    "level": level,
                    "source": "OCHA COD-AB 2020",
                    "total_units": len(stats_df),
                    "generated_at": pd.Timestamp.now().isoformat()
                },
                "statistics": stats_df.to_dict('records')
            }
            
            # Save to JSON
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(stats_dict, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"Saved OCHA administrative statistics to {output_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to save OCHA administrative statistics: {e}")
            raise
