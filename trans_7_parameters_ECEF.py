import os
from osgeo import ogr, osr
import math
import numpy as np


"""
# +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489

WGS72 to WGS84 example: https://epsg.io/9606-method
WGS72_ECEF = np.matrix((3657660.66121,255768.54921,5201382.10819))

"""

def transform_ECEF_coordinates(Towgs84, input_ECEF):
    Tx = Towgs84[0]
    Ty = Towgs84[1]
    Tz = Towgs84[2]
    R1 = Towgs84[3]
    R2 = Towgs84[4]
    R3 = Towgs84[5]
    Sc = Towgs84[6]
    
    input_ECEF = np.array(input_ECEF)

    # Unit conversion:
    Tx_c = Tx  # mm to meter
    Ty_c = Ty
    Tz_c = Tz
    Sc_c = Sc * math.pow(10, -6)  # ppb to ppm/parts
    mastorad = math.pi / (180 * 3600)  # second to arc
    R1_c = mastorad * R1
    R2_c = mastorad * R2
    R3_c = mastorad * R3

    T = np.array([Tx_c, Ty_c, Tz_c])
    R = np.matrix(((1, -R3_c, R2_c), (R3_c, 1, -R1_c), (-R2_c, R1_c, 1)))

    result = T.T + (1 + Sc_c) * np.dot(R, input_ECEF.T)
    
    return tuple(result.tolist()[0]) 


def transform_coordinates(input_coords, input_proj, output_proj):
    # Create input and output spatial reference objects
    input_srs = osr.SpatialReference()
    input_srs.ImportFromProj4(input_proj)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromProj4(output_proj)

    # Create a coordinate transformation
    coord_transform = osr.CoordinateTransformation(input_srs, output_srs)

    # Transform the input coordinates to the output coordinates
    transformed_coords = coord_transform.TransformPoint(*input_coords)

    return transformed_coords

if __name__ == "__main__":
    # Define input coordinates and Proj4 strings
    # OSTN15 BNG -> OSGB36
    input_coords = (450102.724,321429.387,52.495)  # Replace with your own input coordinates
    print("input coordinates",input_coords)
    BNG_proj = "+proj=tmerc +k_0=0.9996012717 +lon_0=-2 +lat_0=49 +x_0=400000 +y_0=-100000 +ellps=airy +units=m"
    OSGB36_proj = "+proj=latlong +ellps=airy"
    # Perform coordinate transformation
    OSGB36_coords = transform_coordinates(input_coords, BNG_proj, OSGB36_proj)
    # OSGB36 to ECEF
    OSGB36_ECEF_proj = "+proj=geocent +datum=OSGB36 +units=m +no_defs"
    OSGB36_ECEF_coords = transform_coordinates(OSGB36_coords,OSGB36_proj,OSGB36_ECEF_proj)

    # OSGB36_ECEF 7 parameters to ETRS89 ECEF
    Towgs84 = [446.448, -125.157, 542.06, 0.15, 0.247, 0.842, -20.489]
    ETRS89_ECEF_coords = transform_ECEF_coordinates(Towgs84, OSGB36_ECEF_coords)

    # ETRS89 ECEF to ETRS89 geographic
    ETRS89_geo_proj = "+proj=longlat +ellps=GRS80 +no_defs +type=crs"
    ETRS89_ECEF_proj = "+proj=geocent +ellps=GRS80 +no_defs +type=crs"

    ETRS89_geo_coords = transform_coordinates(ETRS89_ECEF_coords, ETRS89_ECEF_proj,ETRS89_geo_proj)
    print("ETRS89 transformed",ETRS89_geo_coords)
