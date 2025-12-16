# PoC GIS Technologies

   🌬️ **A Python-based wind energy suitability analysis tool using multicriteria decision making (MCDM) for spatial analysis.**

## 🎯 Purpose

This proof-of-concept generates a Wind Suitability Index (WSI) by analyzing multiple spatial criteria including wind speed, terrain slope, and grid proximity. The system provides automated spatial analysis and reporting for wind energy site assessment.

## 🏗️ Architecture

The project follows **Clean Architecture** principles:

```
src/
├─ domain/                # Business logic
│  ├─ entities.py         # Core entities
│  └─ policies.py         # Business rules
├─ use_cases/             # Application logic
│  ├─ compute_wsi.py      # WSI computation
│  └─ generate_report.py   # Report generation
├─ infrastructure/        # External adapters
│  ├─ rasterio_adapter.py # Raster I/O operations
│  ├─ folium_map.py       # Interactive maps
│  └─ metrics.py          # Performance monitoring
└─ interface/
   └─ cli.py              # Command-line interface
```

## 🚀 Quick Start

### Prerequisites

- Python 3.10+
- pip

### Installation

1. **Clone and setup:**
   ```bash
   git clone <repository-url>
   cd wind-gis-tech-validation
   pip install -r requirements.txt
   ```

2. **Verify installation:**
   ```bash
   python -m src.interface.cli info
   ```

### Basic Usage

1. **Prepare your data:**
   - Place wind speed data in `data/raw/wind_100m.tif`
   - Place elevation data in `data/raw/dem.tif`
   - Place electrical grid data in `data/raw/grid_lines.gpkg`
   - Create AOI in `data/interim/aoi.geojson`

2. **Configure analysis:**
   - Edit `configs/wsi.yaml` with your parameters
   - Adjust weights, thresholds, and output settings

3. **Run analysis:**
   ```bash
   # Compute Wind Suitability Index
   python -m src.interface.cli compute-wsi configs/wsi.yaml
   
   # Generate viability report
   python -m src.interface.cli generate-report configs/wsi.yaml
   
   # Validate configuration
   python -m src.interface.cli validate-config configs/wsi.yaml
   ```

## 📊 Features

- **Multicriteria Analysis**: Weighted sum of normalized spatial criteria
- **Interactive Visualization**: Folium maps with WSI visualization
- **Report Generation**: HTML reports with analysis statistics
- **Performance Monitoring**: Metrics collection and analysis

### Supported Criteria

- **Wind Speed**: Primary factor for wind energy potential
- **Terrain Slope**: Affects turbine installation feasibility
- **Grid Distance**: Proximity to electrical infrastructure
- **Restrictions**: Protected areas, water bodies, urban zones

### Output Formats

- **Raster**: WSI as GeoTIFF (`wsi.tif`)
- **Interactive Map**: HTML with Folium (`wsi_map.html`)
- **Reports**: HTML with analysis results
- **Metrics**: JSON with performance data

## 🔧 Configuration

The system uses YAML configuration files for easy parameter adjustment:

```yaml
# configs/wsi.yaml
aoi: "data/interim/aoi.geojson"
crs: "EPSG:3857"
resolution_m: 100

layers:
  wind: "data/raw/wind_100m.tif"
  dem: "data/raw/dem.tif"
  grid: "data/raw/grid_lines.gpkg"

weights:
  wind: 0.6
  slope: 0.2
  grid_distance: 0.2
```

## 📈 Outputs

- `outputs/rasters/wsi.tif` - Wind Suitability Index raster
- `outputs/maps/wsi_map.html` - Interactive visualization
- `outputs/reports/viability_report.html` - Analysis report
- `outputs/reports/metrics.json` - Performance metrics

## 🧪 Testing

```bash
pytest tests/
```

## 📚 Documentation

- **[Methodology](docs/methodology.md)**: Analysis methodology

---
