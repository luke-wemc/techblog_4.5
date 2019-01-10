# contact luke.sanger@wemcouncil.org

# import packages
import iris
import geopandas as gpd
import numpy as np
import cartopy.crs as ccrs
import shapely
import iris.pandas
import iris.analysis.cartography
import os
import time

# place this file in directory where ERA5 files are located and define path below:
directory_in_str = '/data/private/wemc/100WS'

# read the NUTS shapefile and extract the polygons for a individual countries
nuts=gpd.read_file('/data/private/wemc/ref-nuts-2016-20m.shp/NUTS_RG_20M_2016_4258_LEVL_2.shp')

# get the latitude-longitude grid from LSM file
cubelist2=iris.load('/data/private/wemc/ERA5_Europe_LSM.nc')
cube2=cubelist2[0]
cube2=cube2.intersection(longitude=(-180, 180))

# function to check whether data is within shapefile - source: http://bit.ly/2pKXnWa

def geom_to_masked_cube(cube, geometry, x_coord, y_coord,
                        mask_excludes=False):
    """
    Convert a shapefile geometry into a mask for a cube's data.

    Args:

    * cube:
        The cube to mask.
    * geometry:
        A geometry from a shapefile to define a mask.
    * x_coord: (str or coord)
        A reference to a coord describing the cube's x-axis.
    * y_coord: (str or coord)
        A reference to a coord describing the cube's y-axis.

    Kwargs:

    * mask_excludes: (bool, default False)
        If False, the mask will exclude the area of the geometry from the
        cube's data. If True, the mask will include *only* the area of the
        geometry in the cube's data.

    .. note::
        This function does *not* preserve lazy cube data.

    """
    # Get horizontal coords for masking purposes.
    lats = cube.coord(y_coord).points
    lons = cube.coord(x_coord).points
    lon2d, lat2d = np.meshgrid(lons,lats)

    # Reshape to 1D for easier iteration.
    lon2 = lon2d.reshape(-1)
    lat2 = lat2d.reshape(-1)

    mask = []
    # Iterate through all horizontal points in cube, and
    # check for containment within the specified geometry.
    for lat, lon in zip(lat2, lon2):
        this_point = gpd.geoseries.Point(lon, lat)
        res = geometry.contains(this_point)
        mask.append(res.values[0])

    mask = np.array(mask).reshape(lon2d.shape)
    if mask_excludes:
        # Invert the mask if we want to include the geometry's area.
        mask = ~mask
    # Make sure the mask is the same shape as the cube.
    dim_map = (cube.coord_dims(y_coord)[0],
               cube.coord_dims(x_coord)[0])
    cube_mask = iris.util.broadcast_to_shape(mask, cube.shape, dim_map)

    # Apply the mask to the cube's data.
    data = cube.data
    masked_data = np.ma.masked_array(data, cube_mask)
    cube.data = masked_data
    return cube

directory = os.fsencode(directory_in_str)

# iterate through NUTS file
for i, row in nuts.iterrows():
    
    # select NUTS_ID as a string
    cntry = str(row.NUTS_ID)
    
    # assign row to variable based on NUTS_ID string
    geom = nuts[nuts['NUTS_ID'] == cntry]
    
    # loop through directory of NetCDF files
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".nc"):

            # get the latitude-longitude grid from netcdf file
            cubelist=iris.load(filename)
            cube=cubelist[0]
            cube=cube.intersection(longitude=(-180, 180))

            # apply the geom_to_masked_cube function
            geometry = geom
            masked_cube = geom_to_masked_cube(cube, geometry, 'longitude', 'latitude', mask_excludes=True)
            
            # apply the geom_to_masked_cube function to the LSM
            masked_cube2 = geom_to_masked_cube(cube2, geometry, 'longitude', 'latitude', mask_excludes=True)

            # use iris function to get area weights
            grid_areas = iris.analysis.cartography.cosine_latitude_weights(masked_cube)
            grid_areaslsm = iris.analysis.cartography.cosine_latitude_weights(masked_cube2)

            grid_lsm = grid_areas * grid_areaslsm

            # create new time series cube with area averages
            new_cube = masked_cube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_lsm)

            # convert to series and retain the time series - use this!
            dfs = iris.pandas.as_series(new_cube, copy=True)

            #generate a new filename
            edit = str(filename)
            edit2 = edit.rstrip(".nc")
            csv_name = cntry + "_" + edit2 + ".csv"

            # save as csv
            dfs.to_csv(csv_name)

            # pause loop to save CPU drain
            time.sleep(0.2)
