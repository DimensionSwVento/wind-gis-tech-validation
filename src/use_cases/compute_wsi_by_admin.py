"""
Compute WSI statistics by administrative boundaries use case.

This module handles the calculation of Wind Suitability Index (WSI) statistics
for administrative units (departments and municipalities).
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from typing import Dict, List, Optional, Tuple, Union
import logging
from pathlib import Path
import json

from ..infrastructure.admin_boundaries import AdminBoundariesAdapter
from ..infrastructure.ocha_adapter import OCHAAdapter
from ..infrastructure.rasterio_adapter import RasterioAdapter

logger = logging.getLogger(__name__)


class ComputeWSIByAdminUseCase:
    """Use case for computing WSI statistics by administrative boundaries."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.admin_adapter = AdminBoundariesAdapter()
        self.ocha_adapter = OCHAAdapter()
        self.raster_adapter = RasterioAdapter()
    
    def execute(self, config: Dict) -> Dict:
        """
        Execute WSI computation by administrative boundaries.
        
        Args:
            config: Configuration dictionary with admin_boundaries settings
            
        Returns:
            Dictionary with execution results
        """
        try:
            self.logger.info("Starting WSI computation by administrative boundaries")
            
            # Extract configuration
            admin_config = config.get('admin_boundaries', {})
            wsi_raster_path = config.get('wsi_raster_path', 'outputs/rasters/wsi.tif')
            output_dir = config.get('output_dir', 'outputs')
            viability_threshold = config.get('thresholds', {}).get('viability_threshold', 0.5)
            crs = config.get('crs', 'EPSG:3857')
            
            # Load administrative boundaries (try OCHA first, fallback to generic)
            admin_gdf = self._load_admin_boundaries(admin_config, crs)
            
            # Load WSI raster
            wsi_raster = self._load_wsi_raster(wsi_raster_path)
            
            # Reproject admin geometries to raster CRS before analysis
            try:
                admin_gdf = admin_gdf.to_crs(self._get_raster_crs(wsi_raster))
            except Exception as e:
                self.logger.warning(f"Failed to reproject admin boundaries to raster CRS: {e}")
            
            # Compute statistics for each administrative unit
            stats_df = self._compute_admin_statistics(
                admin_gdf,
                wsi_raster, 
                admin_config,
                viability_threshold
            )
            
            # Save cropped rasters
            self._save_cropped_rasters(
                admin_gdf,
                wsi_raster, 
                admin_config,
                output_dir
            )
            
            # Save statistics
            is_ocha = 'col_admbnda' in admin_config.get('file', '') or 'ocha' in admin_config.get('file', '').lower()
            self._save_statistics(stats_df, output_dir, admin_config.get('target_level', 'admin'), is_ocha)
            
            self.logger.info("WSI computation by administrative boundaries completed successfully")
            
            return {
                "success": True,
                "total_units": len(stats_df),
                "statistics_file": f"{output_dir}/reports/wsi_stats_by_{admin_config.get('target_level', 'admin')}.json",
                "rasters_dir": f"{output_dir}/rasters/wsi_by_{admin_config.get('target_level', 'admin')}/"
            }
            
        except Exception as e:
            self.logger.error(f"Failed to compute WSI by administrative boundaries: {e}")
            raise
    
    def _load_admin_boundaries(self, admin_config: Dict, crs: str) -> gpd.GeoDataFrame:
        """Load and process administrative boundaries (OCHA or generic)."""
        try:
            file_path = admin_config.get('file')
            layer = admin_config.get('layer')
            
            if not file_path:
                raise ValueError("admin_boundaries.file is required in configuration")
            
            # Check if this is an OCHA file
            if 'col_admbnda' in file_path or 'ocha' in file_path.lower():
                self.logger.info("Detected OCHA COD-AB file, using OCHA adapter")
                
                # Load OCHA boundaries
                admin_gdf = self.ocha_adapter.load_ocha_boundaries(file_path, crs)
                
                # Apply filters if specified
                filters = admin_config.get('filter', {})
                if filters:
                    admin_gdf = self.ocha_adapter.filter_ocha_boundaries(admin_gdf, filters)
                
                # Validate required fields
                required_fields = [
                    admin_config.get('id_field', 'ADM2_PCODE'),
                    admin_config.get('name_field', 'ADM2_ES')
                ]
                
                if not self.ocha_adapter.validate_ocha_data(admin_gdf, required_fields):
                    raise ValueError("OCHA administrative boundaries data validation failed")
                
            else:
                self.logger.info("Using generic admin boundaries adapter")
                
                # Load generic boundaries
                admin_gdf = self.admin_adapter.load_admin_boundaries(file_path, layer, crs)
                
                # Apply filters if specified
                filters = admin_config.get('filter', {})
                if filters:
                    admin_gdf = self.admin_adapter.filter_admin_boundaries(admin_gdf, filters)
                
                # Validate required fields
                required_fields = [
                    admin_config.get('id_field', 'id'),
                    admin_config.get('name_field', 'name')
                ]
                
                if not self.admin_adapter.validate_admin_data(admin_gdf, required_fields):
                    raise ValueError("Administrative boundaries data validation failed")
            
            return admin_gdf
            
        except Exception as e:
            self.logger.error(f"Failed to load administrative boundaries: {e}")
            raise
    
    def _load_wsi_raster(self, wsi_raster_path: str) -> rasterio.DatasetReader:
        """Load WSI raster."""
        try:
            if not Path(wsi_raster_path).exists():
                raise FileNotFoundError(f"WSI raster not found: {wsi_raster_path}")
            
            return rasterio.open(wsi_raster_path)
            
        except Exception as e:
            self.logger.error(f"Failed to load WSI raster: {e}")
            raise
    
    def _compute_admin_statistics(self, 
                                admin_gdf: gpd.GeoDataFrame,
                                wsi_raster: rasterio.DatasetReader,
                                admin_config: Dict,
                                viability_threshold: float) -> pd.DataFrame:
        """Compute WSI statistics for each administrative unit."""
        try:
            self.logger.info(f"Computing statistics for {len(admin_gdf)} administrative units")
            
            id_field = admin_config.get('id_field', 'id')
            name_field = admin_config.get('name_field', 'name')
            parent_field = admin_config.get('parent_field')
            
            stats_list = []
            
            # Ensure same CRS as raster
            try:
                if admin_gdf.crs and admin_gdf.crs != wsi_raster.crs:
                    admin_gdf = admin_gdf.to_crs(wsi_raster.crs)
            except Exception as e:
                self.logger.warning(f"Could not align CRS before stats: {e}")

            raster_bounds = wsi_raster.bounds
            raster_bbox = gpd.GeoSeries([gpd.box(*raster_bounds)], crs=wsi_raster.crs)

            for idx, row in admin_gdf.iterrows():
                try:
                    # Get geometry
                    geometry = row.geometry
                    if geometry is None or geometry.is_empty:
                        continue
                    # Fix invalid geometries
                    if not geometry.is_valid:
                        geometry = geometry.buffer(0)
                    # Skip if no bbox overlap
                    if not gpd.GeoSeries([geometry], crs=admin_gdf.crs).intersects(raster_bbox.iloc[0]).iloc[0]:
                        self.logger.warning(f"No overlap with raster for {row.get(name_field, 'unknown')} - skipping")
                        continue
                    
                    # Mask WSI raster by geometry
                    try:
                        masked_data, transform = mask(
                            wsi_raster,
                            [geometry],
                            crop=True,
                            nodata=np.nan,
                            all_touched=True
                        )
                    except ValueError as ve:
                        # Common when slight misalignment occurs
                        self.logger.warning(f"Masking failed for {row.get(name_field, 'unknown')}: {ve}")
                        continue
                    
                    # Flatten and remove nodata values
                    wsi_values = masked_data[0].flatten()
                    valid_mask = ~np.isnan(wsi_values)
                    valid_values = wsi_values[valid_mask]
                    
                    if len(valid_values) == 0:
                        self.logger.warning(f"No valid WSI values for {row[name_field]}")
                        continue
                    
                    # Compute statistics
                    stats = {
                        id_field: row[id_field],
                        name_field: row[name_field],
                        'mean_wsi': float(np.mean(valid_values)),
                        'p90_wsi': float(np.percentile(valid_values, 90)),
                        'max_wsi': float(np.max(valid_values)),
                        'min_wsi': float(np.min(valid_values)),
                        'std_wsi': float(np.std(valid_values)),
                        'valid_coverage': float(valid_mask.sum() / len(wsi_values)) * 100,
                        'pct_apto': float(np.sum(valid_values >= viability_threshold) / len(valid_values)) * 100,
                        'total_pixels': int(len(valid_values))
                    }
                    
                    # Add parent field if specified
                    if parent_field and parent_field in row:
                        stats[parent_field] = row[parent_field]
                    
                    stats_list.append(stats)
                    
                except Exception as e:
                    self.logger.warning(f"Failed to compute statistics for {row.get(name_field, 'unknown')}: {e}")
                    continue
            
            # Create DataFrame
            stats_df = pd.DataFrame(stats_list)
            
            if len(stats_df) == 0:
                raise ValueError("No valid statistics computed")
            
            self.logger.info(f"Computed statistics for {len(stats_df)} administrative units")
            return stats_df
            
        except Exception as e:
            self.logger.error(f"Failed to compute administrative statistics: {e}")
            raise
    
    def _save_cropped_rasters(self, 
                            admin_gdf: gpd.GeoDataFrame,
                            wsi_raster: rasterio.DatasetReader,
                            admin_config: Dict,
                            output_dir: str) -> None:
        """Save cropped WSI rasters for each administrative unit."""
        try:
            target_level = admin_config.get('target_level', 'admin')
            id_field = admin_config.get('id_field', 'id')
            name_field = admin_config.get('name_field', 'name')
            
            # Create output directory
            rasters_dir = Path(output_dir) / 'rasters' / f'wsi_by_{target_level}'
            rasters_dir.mkdir(parents=True, exist_ok=True)
            
            self.logger.info(f"Saving cropped rasters to {rasters_dir}")
            
            # Ensure same CRS as raster
            try:
                if admin_gdf.crs and admin_gdf.crs != wsi_raster.crs:
                    admin_gdf = admin_gdf.to_crs(wsi_raster.crs)
            except Exception as e:
                self.logger.warning(f"Could not align CRS before saving crops: {e}")

            raster_bounds = wsi_raster.bounds
            raster_bbox = gpd.GeoSeries([gpd.box(*raster_bounds)], crs=wsi_raster.crs)

            for idx, row in admin_gdf.iterrows():
                try:
                    # Get geometry
                    geometry = row.geometry
                    if geometry is None or geometry.is_empty:
                        continue
                    if not geometry.is_valid:
                        geometry = geometry.buffer(0)
                    if not gpd.GeoSeries([geometry], crs=admin_gdf.crs).intersects(raster_bbox.iloc[0]).iloc[0]:
                        continue
                    
                    # Mask WSI raster by geometry
                    masked_data, transform = mask(
                        wsi_raster,
                        [geometry],
                        crop=True,
                        nodata=np.nan,
                        all_touched=True
                    )
                    
                    # Create output filename
                    unit_id = str(row[id_field]).replace(' ', '_').replace('/', '_')
                    unit_name = str(row[name_field]).replace(' ', '_').replace('/', '_')
                    output_file = rasters_dir / f'wsi_{unit_id}_{unit_name}.tif'
                    
                    # Save cropped raster
                    with rasterio.open(
                        output_file,
                        'w',
                        driver='GTiff',
                        height=masked_data.shape[1],
                        width=masked_data.shape[2],
                        count=1,
                        dtype=masked_data.dtype,
                        crs=wsi_raster.crs,
                        transform=transform,
                        nodata=np.nan
                    ) as dst:
                        dst.write(masked_data[0], 1)
                    
                except Exception as e:
                    self.logger.warning(f"Failed to save cropped raster for {row.get(name_field, 'unknown')}: {e}")
                    continue
            
            self.logger.info(f"Saved cropped rasters to {rasters_dir}")
            
        except Exception as e:
            self.logger.error(f"Failed to save cropped rasters: {e}")
            raise

    def _get_raster_crs(self, wsi_raster: rasterio.DatasetReader):
        """Safely get raster CRS object or EPSG string."""
        try:
            return wsi_raster.crs
        except Exception:
            return None
    
    def _save_statistics(self, 
                        stats_df: pd.DataFrame, 
                        output_dir: str, 
                        level: str,
                        is_ocha: bool = False) -> None:
        """Save statistics to JSON file."""
        try:
            # Create output directory
            reports_dir = Path(output_dir) / 'reports'
            reports_dir.mkdir(parents=True, exist_ok=True)
            
            # Save using appropriate adapter
            if is_ocha:
                self.ocha_adapter.save_ocha_stats(stats_df, str(reports_dir), level)
            else:
                self.admin_adapter.save_admin_stats(stats_df, str(reports_dir), level)
            
        except Exception as e:
            self.logger.error(f"Failed to save statistics: {e}")
            raise
