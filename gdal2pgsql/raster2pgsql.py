#! /usr/bin/env python
#
# This is a simple utility used to dump GDAL dataset into HEX WKB stream.
# It's considered as a prototype of raster2pgsql tool planned to develop
# in future.
# For more details about raster2pgsql tool, see Specification page:
# http://trac.osgeo.org/postgis/wiki/WKTRaster
#
# The script requires Python bindings for GDAL.
# Available at http://trac.osgeo.org/gdal/wiki/GdalOgrInPython
#
################################################################################
# Copyright (C) 2009-2010 Mateusz Loskot <mateusz@loskot.net>
# Copyright (C) 2009-2011 Pierre Racine <pierre.racine@sbf.ulaval.ca>
# Copyright (C) 2009-2010 Jorge Arevalo <jorge.arevalo@deimos-space.com>
# Copyright (C) 2021      Diego Moreira Carvalho <diego@khartes.com.br>
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
################################################################################
import binascii
import math
import osgeo.gdalconst as gdalc
import sys
from osgeo import gdal
from osgeo import osr

# Default version of WKTRaster protocol.
# WARNING: Currently, this is the only valid value
# and option -w, --raster-version is ignored, if specified.
g_rt_version = 0

# Default format of binary output is little-endian (NDR)
# WARNING: Currently, big-endian (XDR) output is not supported
# and option -e, --endian is ignored, if specified.
g_rt_endian = 1

################################################################################
# UTILITIES

def is_nan(x):
    if sys.hexversion < 0x02060000:
        return str(float(x)).lower() == 'nan'
    else:
        return math.isnan(x)  

def gdt2pt(gdt):
    """Translate GDAL data type to WKT Raster pixel type."""
    pixtypes = {
        gdalc.GDT_Byte    : { 'name': 'PT_8BUI',  'id':  4 },
        gdalc.GDT_Int16   : { 'name': 'PT_16BSI', 'id':  5 },
        gdalc.GDT_UInt16  : { 'name': 'PT_16BUI', 'id':  6 },
        gdalc.GDT_Int32   : { 'name': 'PT_32BSI', 'id':  7 },
        gdalc.GDT_UInt32  : { 'name': 'PT_32BUI', 'id':  8 },
        gdalc.GDT_Float32 : { 'name': 'PT_32BF',  'id': 10 },
        gdalc.GDT_Float64 : { 'name': 'PT_64BF',  'id': 11 }
        }    
    # XXX: Uncomment these logs to debug types translation
    return pixtypes.get(gdt, 13)
    
def pt2fmt(pt):
    """Returns binary data type specifier for given pixel type."""
    fmttypes = {
        4: 'B', # PT_8BUI
        5: 'h', # PT_16BSI
        6: 'H', # PT_16BUI
        7: 'i', # PT_32BSI
        8: 'I', # PT_32BUI
        10: 'f', # PT_32BF
        11: 'd'  # PT_64BF
        }
    return fmttypes.get(pt, 'x')

################################################################################
# SQL OPERATIONS

def quote_sql_value(value):
    assert value is not None, "None value given"

    if len(value) > 0 and value[0] != "'" and value[:-1] != "'":
        sql = "'" + str(value) + "'"
    else:
        sql = value
    return sql

def quote_sql_name(name):
    assert name is not None, "None name given"

    if name[0] != "\"" and name[:-1] != "\"":
        sql = "\"" + str(name) + "\""
    else:
        sql = name
    return sql

def make_sql_value_array(values):
    sql = "ARRAY["
    for v in values:
        if isinstance(v, str):
            sql += quote_sql_value(v) + ","
        else:
            sql += str(v) + ','
    sql = sql[:-1] # Trim comma
    sql += "]"
    return sql

def make_sql_schema_table_names(schema_table):
    st = schema_table.split('.')
    if len(st) == 1:
        st.insert(0, "public")
    assert len(st) == 2, "Invalid format of table name, expected [<schema>.]table"
    return (st[0], st[1])

def make_sql_full_table_name(schema_table):
    st = make_sql_schema_table_names(schema_table)
    table = "\"%s\".\"%s\"" % (st[0], st[1])
    return table

def make_sql_drop_table(schema, table):
    return f"DROP TABLE IF EXISTS {quote_sql_name(schema)}.{quote_sql_name(table)} CASCADE;\n"

def make_sql_create_table(schema, table, column):
    return f"CREATE TABLE {quote_sql_name(schema)}.{quote_sql_name(table)} (rid serial PRIMARY KEY, {quote_sql_name(column)} raster, \"filename\" text);\n"

def make_sql_create_gist(schema, table, column):
    sql = f"CREATE INDEX ON {quote_sql_name(schema)}.{quote_sql_name(table)} USING gist (st_convexhull({quote_sql_name(column)}));\n"
    sql += f"ANALYZE {quote_sql_name(schema)}.{quote_sql_name(table)};\n"
    return  sql

def make_sql_insert_raster(schema, table, column, hexwkb, file_name):
    sql = f"INSERT INTO {quote_sql_name(schema)}.{quote_sql_name(table)} "
    sql += f"({quote_sql_name(column)},{quote_sql_name('filename')}) "
    sql += f"VALUES (('{hexwkb.decode('utf-8')}')::raster, ('{file_name}')::text);\n"
    return sql

def make_sql_addrasterconstraints(schema, table, column):
    srid=True 
    scale_x=True 
    scale_y=True 
    blocksize_x=True 
    blocksize_y=True 
    same_alignment=True 
    regular_blocking=False 
    num_bands=True  
    pixel_types=True  
    nodata_values=True  
    out_db=True  
    extent=True 
    list_ = [srid, scale_x, scale_y, blocksize_x, blocksize_y, same_alignment, regular_blocking, num_bands, pixel_types, nodata_values, out_db, extent]
    sql = f"SELECT AddRasterConstraints('{schema}','{table}','{column}',{','.join([str(i).upper() for i in list_])});\n"

    return sql   

################################################################################
# RASTER OPERATIONS

def collect_pixel_types(ds, band_from, band_to):
    """Collect pixel types of bands in requested range.
       Use names of pixel types in format as returned
       by rt_core function rt_pixtype_name()"""

    pt =[]
    for i in range(band_from, band_to):
        band = ds.GetRasterBand(i)
        pixel_type = gdt2pt(band.DataType)['name'][3:]
        pt.append(pixel_type)
    
    return pt

def collect_nodata_values(ds, band_from, band_to):
    """Collect nodata values of bands in requested range"""

    nd = []
    for i in range(band_from, band_to):
        band = ds.GetRasterBand(i)
        nodata = band.GetNoDataValue()
        if nodata is not None and not is_nan(nodata):
            nd.append(nodata)

    return nd    

def calculate_grid_size(raster_size, block_size):
     """Dimensions of grid made up with blocks of requested size"""
 
     # Exact number of grid dimensions
     nx = float(raster_size[0]) / float(block_size[0])
     ny = float(raster_size[1]) / float(block_size[1])
 
     return ( int(math.ceil(nx)), int(math.ceil(ny)))

def calculate_block_pad_size(raster_size, xoff, yoff, block_size):
    """Calculates number of columns [0] and rows [1] of padding""" 
    assert raster_size is not None

    xpad = 0
    ypad= 0
    block_bound = ( xoff + block_size[0], yoff + block_size[1] )

    if block_bound[0] > raster_size[0]:
        xpad = block_bound[0] - raster_size[0]
    if block_bound[1] > raster_size[1]:
        ypad = block_bound[1] - raster_size[1]

    return (xpad, ypad)

def calculate_geoxy(gt, xy):
    """Calculate georeferenced coordinate from given x and y"""
    assert gt is not None
    assert xy is not None
    assert len(xy) == 2

    xgeo = gt[0] + gt[1] * xy[0] + gt[2] * xy[1];
    ygeo = gt[3] + gt[4] * xy[0] + gt[5] * xy[1];

    return (xgeo, ygeo)

def calculate_bounding_box(ds, gt):
    """Calculate georeferenced coordinates of spatial extent of raster dataset"""
    assert ds is not None

    # UL, LL, UR, LR
    dim = ( (0,0),(0,ds.RasterYSize),(ds.RasterXSize,0),(ds.RasterXSize,ds.RasterYSize) )

    ext = (calculate_geoxy(gt, dim[0]), calculate_geoxy(gt, dim[1]),
           calculate_geoxy(gt, dim[2]), calculate_geoxy(gt, dim[3]))

    return ext

def check_hex(hex, bytes_size = None):
    assert hex is not None, "Error: Missing hex string"
    size = len(hex)
    assert size > 0, "Error: hex string is empty"
    assert size % 2 == 0, "Error: invalid size of hex string"
    if bytes_size is not None:
        n = int(size / 2)
        assert n == bytes_size, "Error: invalid number of bytes %d, expected %d" % (n, bytes_size)

def fetch_band_nodata(band, default = 0):
    assert band is not None

    nodata = default
    if band.GetNoDataValue() is not None:
        nodata = band.GetNoDataValue()
    # else:
    #     print("WARNING: No nodata flagged in raster_columns metadata. "
    #           "In serialized raster, nodata bytes will have value of 0.\n")
    return nodata

def wkblify(fmt, data):
    """Writes raw binary data into HEX-encoded string using binascii module."""
    import struct
    # Binary to HEX
    fmt_little = '<' +fmt
    hexstr = binascii.hexlify(struct.pack(fmt_little, data)).upper()
    return hexstr

def wkblify_raster_header(ds, level, ulp, xsize = None, ysize = None):
    """Writes WKT Raster header based on given GDAL into HEX-encoded WKB."""
    assert ds is not None, "Error: Missing GDAL dataset"
    assert level >= 1
    assert len(ulp) == 2, "Error: invalid upper-left corner"

    if xsize is None or ysize is None:
        assert xsize is None and ysize is None
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize

    # Collect GeoReference information
    gt = tuple(list(ds.GetGeoTransform()))
    ul = calculate_geoxy(gt, (ulp[0], ulp[1]))
    rt_ip = ( ul[0], ul[1] )
    rt_skew = ( gt[2], gt[4] )
    rt_scale = ( gt[1] * level, gt[5] * level )
    
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    srid = int(srs.GetAttrValue("AUTHORITY", 1))

    # Burn input raster as WKTRaster WKB format
    hexwkb = b''
    ### Endiannes
    hexwkb += wkblify('B', g_rt_endian)
    ### Version
    hexwkb += wkblify('H', g_rt_version)
    ### Number of bands
    hexwkb += wkblify('H', ds.RasterCount)
    check_hex(hexwkb, 5)
    ### Georeference
    hexwkb += wkblify('d', rt_scale[0])
    hexwkb += wkblify('d', rt_scale[1])
    hexwkb += wkblify('d', rt_ip[0])
    hexwkb += wkblify('d', rt_ip[1])
    hexwkb += wkblify('d', rt_skew[0])
    hexwkb += wkblify('d', rt_skew[1])
    hexwkb += wkblify('i', srid)
    check_hex(hexwkb, 57)
    ### Number of columns and rows
    hexwkb += wkblify('H', xsize)
    hexwkb += wkblify('H', ysize)
    check_hex(hexwkb, 61)
    return hexwkb

def wkblify_band_header(band):
    """Writes band header into HEX-encoded WKB"""
    assert band is not None, "Error: Missing GDAL raster band"

    hexwkb = b""

    first4bits = 0   
       
    nodata = band.GetNoDataValue()
    # If there is no nodata value, set it to 0. Otherwise set the HasNodata bit to 1
    if nodata is not None:
        first4bits += 64
    else:
        nodata = 0
    
    # Encode pixel type
    pixtype = gdt2pt(band.DataType)['id']
    hexwkb += wkblify('B', pixtype + first4bits)
    
    # Encode nodata value (or Zero, if nodata unavailable) 
    hexwkb += wkblify(pt2fmt(pixtype), nodata)

    check_hex(hexwkb)
    return hexwkb

def wkblify_band(band, level, xoff, yoff, block_size):
    """Writes band of given GDAL dataset into HEX-encoded WKB for WKT Raster output."""
    assert band is not None, "Error: Missing GDAL raster band"

    hexwkb = ''
 
    pixels = band.ReadAsArray(xoff, yoff, block_size[0], block_size[1],
                                block_size[0], block_size[1])

    hexwkb = binascii.hexlify(pixels)

    check_hex(hexwkb)
    return hexwkb

def get_raster_as_hexwkb(execute_sql_insert_raster, file_name):
    """Writes given raster dataset using GDAL features into HEX-encoded of
    WKB for WKT Raster output."""
    
    assert file_name is not None, "Input file is none, expected file name"

    ds = gdal.Open(file_name, gdalc.GA_ReadOnly);
    if ds is None:
        raise Exception(f"Error: Cannot open input file: {file_name}")

    band_from, band_to = 1, ds.RasterCount + 1
    raster_size = ( ds.RasterXSize, ds.RasterYSize )
    block_size  =  (2048,2048)
    grid_size  =   calculate_grid_size(raster_size, block_size)
    
    for ycell in range(0, grid_size[1]):
        for xcell in range(0, grid_size[0]):
            hexwkb = b''
            xoff = xcell * block_size[0]
            yoff = ycell * block_size[1]
            read_padding_size = calculate_block_pad_size(( ds.RasterXSize, ds.RasterYSize ), xoff, yoff, block_size)
            valid_read_block_size = ( block_size[0] - read_padding_size[0],
                                        block_size[1] - read_padding_size[1] )


            hexwkb += wkblify_raster_header(ds, 1, (xoff, yoff),valid_read_block_size[0], valid_read_block_size[1] )

            for b in range(band_from, band_to):
                band = ds.GetRasterBand(b)
                assert band is not None, "Missing GDAL raster band %d" % b

                hexwkb += wkblify_band_header(band)
                hexwkb += wkblify_band(band, 1, xoff, yoff, valid_read_block_size)
            check_hex(hexwkb) 
            # TODO: Remove to not to sriddecrease performance
            execute_sql_insert_raster(hexwkb)

################################################################################
# INTERFACE

def raster2pgsql(raster_file_name, sql_file_name, table, schema='public', column="rast"):
    output = open(sql_file_name, "w")
    output.write(raster2pgsql(raster_file_name, table, schema, column))
    output.close()

def raster2pgsql(execute_sql, raster_file_name, table, schema='public', column="rast"):
    def execute_sql_insert_raster(hexwkb):
        execute_sql(make_sql_insert_raster(schema, table, column, hexwkb, raster_file_name))
    execute_sql('BEGIN;\n')
    execute_sql(make_sql_drop_table(schema, table))
    execute_sql(make_sql_create_table(schema, table, column))

    get_raster_as_hexwkb(execute_sql_insert_raster, raster_file_name)

    execute_sql(make_sql_create_gist(schema, table, column))
    execute_sql(make_sql_addrasterconstraints(schema, table, column))

    execute_sql('END;\n')