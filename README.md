# Fix_proj_gdal_7_parameters

Since gdaltranslate and proj4js do not convert ellipsoid heights between different ellipsoid / datum directly, it is still be able to do.

This is a python script to apply the 7 parameters transformation from one datum to the WGS84 equivalent datum, in any way, say ETRS89.

I was using the example from WGS72 to WGS84 example: https://epsg.io/9606-method and a futhre trial from OSTN15 to OSGB36 / BNG and ETRS89.