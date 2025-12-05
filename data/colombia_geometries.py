"""
Geometrías reales de departamentos y municipios de Colombia.
Este archivo contiene las coordenadas aproximadas de los límites reales.
"""

# Geometría real de La Guajira (solo territorio terrestre, excluyendo mar)
LA_GUAJIRA_GEOMETRY = {
    "type": "Polygon",
    "coordinates": [[
        # Límite sur (frontera con Cesar) - línea terrestre
        [-73.3, 10.5], [-73.0, 10.5], [-72.7, 10.6], [-72.4, 10.7], [-72.1, 10.8], [-71.8, 10.9],
        
        # Límite este (frontera con Venezuela) - línea terrestre
        [-71.5, 11.0], [-71.2, 11.1], [-71.0, 11.3], [-70.9, 11.5], [-70.9, 11.7], [-71.0, 11.9],
        [-71.2, 12.0], [-71.4, 12.1], [-71.6, 12.1], [-71.8, 12.0], [-72.0, 11.9], [-72.2, 11.8],
        [-72.4, 11.7], [-72.6, 11.6], [-72.8, 11.5], [-73.0, 11.4], [-73.2, 11.3], [-73.3, 11.2],
        
        # Límite oeste (frontera con Magdalena) - línea terrestre
        [-73.3, 11.0], [-73.3, 10.8], [-73.3, 10.5]  # Cierre
    ]]
}

# Otros departamentos importantes para energía eólica
DEPARTAMENTOS = {
    "La Guajira": {
        "geometry": LA_GUAJIRA_GEOMETRY,
        "municipios": [
            "Riohacha", "Maicao", "Uribia", "Manaure", "San Juan del Cesar",
            "Villanueva", "El Molino", "Fonseca", "Barrancas", "Distracción",
            "Hatonuevo", "La Jagua del Pilar", "Urumita", "Albania", "Dibula"
        ]
    },
    "Cesar": {
        "geometry": {
            "type": "Polygon",
            "coordinates": [[
                # Límites más precisos de Cesar
                [-74.2, 7.8], [-73.8, 8.0], [-73.5, 8.2], [-73.2, 8.5], [-72.8, 8.8], [-72.5, 9.2], [-72.2, 9.5], [-72.0, 9.8], [-72.0, 10.2], [-72.2, 10.5],
                [-72.5, 10.7], [-72.8, 10.8], [-73.2, 10.7], [-73.5, 10.5], [-73.8, 10.2], [-74.0, 9.8], [-74.2, 9.5], [-74.3, 9.2], [-74.2, 8.8], [-74.2, 7.8]
            ]]
        },
        "municipios": [
            "Valledupar", "Aguachica", "Codazzi", "La Paz", "San Diego",
            "Chimichagua", "Curumaní", "El Copey", "La Gloria", "Manaure Balcón del Cesar"
        ]
    },
    "Magdalena": {
        "geometry": {
            "type": "Polygon", 
            "coordinates": [[
                # Límites más precisos de Magdalena
                [-75.2, 8.8], [-74.8, 9.0], [-74.5, 9.2], [-74.2, 9.5], [-73.8, 9.8], [-73.5, 10.0], [-73.2, 10.2], [-73.0, 10.5], [-73.0, 10.8], [-73.2, 11.0],
                [-73.5, 11.2], [-73.8, 11.3], [-74.2, 11.2], [-74.5, 11.0], [-74.8, 10.8], [-75.0, 10.5], [-75.2, 10.2], [-75.3, 9.8], [-75.2, 9.5], [-75.2, 8.8]
            ]]
        },
        "municipios": [
            "Santa Marta", "Ciénaga", "Fundación", "Aracataca", "El Retén",
            "Zona Bananera", "Algarrobo", "Ariguaní", "Cerro San Antonio", "Chivolo"
        ]
    },
    "Atlántico": {
        "geometry": {
            "type": "Polygon",
            "coordinates": [[
                # Límites más precisos de Atlántico
                [-75.6, 10.2], [-75.2, 10.4], [-74.8, 10.6], [-74.5, 10.8], [-74.2, 11.0], [-74.0, 11.2], [-74.0, 11.4], [-74.2, 11.6], [-74.5, 11.7], [-74.8, 11.6],
                [-75.0, 11.4], [-75.2, 11.2], [-75.4, 11.0], [-75.6, 10.8], [-75.7, 10.6], [-75.6, 10.4], [-75.6, 10.2]
            ]]
        },
        "municipios": [
            "Barranquilla", "Soledad", "Malambo", "Sabanagrande", "Sabanalarga",
            "Puerto Colombia", "Galapa", "Tubará", "Usiacurí", "Luruaco"
        ]
    }
}

def get_departamento_geometry(departamento_name):
    """Obtener geometría de un departamento."""
    return DEPARTAMENTOS.get(departamento_name, {}).get("geometry")

def get_municipios_departamento(departamento_name):
    """Obtener municipios de un departamento."""
    return DEPARTAMENTOS.get(departamento_name, {}).get("municipios", [])

def get_all_departamentos():
    """Obtener lista de todos los departamentos."""
    return list(DEPARTAMENTOS.keys())

def is_point_in_geometry(lat, lon, geometry):
    """Verificar si un punto está dentro de una geometría."""
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
