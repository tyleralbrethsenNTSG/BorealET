'''
Script to assemble the inputs necessary for Boreal ET

----------------------------------------------------------
Requirements:
ROI: region of interest (ee.FeatureCollection)
params: from calibration results
Access to HISTARFM data ('projects/borealet/assets/histarfm_lai',
                         'projects/borealet/assets/histarfm_fpar')
    Email tyler.albrethsen@mso.umt.edu for permission.
----------------------------------------------------------

Relies on data from:
    MERRA-2 (NASA)
    MODIS/VIIRS (NASA)
    3D Elevation Program (USGS)
    Canadian Digital Elevation Model (Natural Resources Canada)
    ABoVE: Landsat-derived Annual Dominant Land Cover (NASA)
    HISTARFM (ISP Valencia)
    

Returns an ee.ImageCollection with one image per day and bands for each of the following inputs:
    'LAI',
    'FPAR',
    'Albedo',
    'Elevation',
    'Landcover',
    'Annual_Mean_Temp',
    'SWGDN_day',
    'LWGNT_day',
    'SWGNT_day',
    'PS_day',
    'QV10M_day',
    'T10M_day',
    'SWGDN_night',
    'LWGNT_night',
    'SWGNT_night',
    'PS_night',
    'QV10M_night',
    'T10M_night',
    'Tmin',
    'tmin_close',
    'tmin_open',
    'vpd_open',
    'vpd_close',
    'gl_sh',
    'gl_wv',
    'g_cuticular',
    'csl',
    'rbl_min',
    'rbl_max',
    'beta'

'''

import numpy as np
import datetime
from datetime import timedelta
import ee
ee.Authenticate()
ee.Initialize(project = 'borealet')

ROI = ee.FeatureCollection('projects/borealet/assets/BorealET_ROI').geometry()
DAY_NIGHT_THRESHOLD = 1.0   # Threshold for incoming short-wave radiation to determine daytime (w/m2)

# Params for BorealET_1.0 are based on Calibration_Results_20240507
params = {
    'tmin_close_list': [-8,-6,-8,-7,-8,-8,-999,-8,-999],
    'tmin_open_list':  [5,21.4,3.6,17.8,5,19.4,-999,21.5,-999],
    'vpd_open_list': [650,650,650,650,650,650,-999,650,-999],
    'vpd_close_list': [4156,3434,4636,3923,4156,4673,-999,4963,-999],
    'gl_sh_list': [0.04018,0.014263667,0.090083,0.030417667,0.04018,0.020077,-999,0.034125333,-999],
    'gl_wv_list': [0.020606667,0.032551,0.017276,0.035062667,0.020606667,0.027628667,-999,0.034698667,-999],
    'g_cuticular_list': [0.000151,0.00006,0.000061,0.000053,0.000151,0.00004,-999,0.000051,-999],
    'csl_list': [0.002994667,0.005932333,0.004523333,0.007365333,0.002994667,0.00462,-999,0.006002,-999],
    'rbl_min_list': [734,620,1461.333333,914.3333333,734,1318.666667,-999,490.6666667,-999],
    'rbl_max_list': [2000,2000,2000,2000,2000,2000,-999,2000,-999],
    'beta_list': [1017,3951.333333,741.3333333,1500.666667,1017,770.3333333,-999,5746,-999]
}

def roi_clip(image):
    return image.clip(ROI)


#-----------------------------------
# Define Input Collection
#-----------------------------------


# Weather Reannalysis (MERRA-2)
#-----------------------------------

# Create Masks for Daytime and Nighttime based on downwelling shortwave radiation
def day_mask(image):
    day_night_threshold = 1.0           # Threshold for incoming short-wave radiation to determine daytime (w/m2)
    return image.gte(day_night_threshold)
def night_mask(image):
    day_night_threshold = 1.0           # Threshold for incoming short-wave radiation to determine daytime (w/m2)
    return image.lt(day_night_threshold)


# Define Aggregation Functions for Daily Averages
def aggregate_merra2_daytime(image_collection, mask_list):
    image_list = image_collection.toList(image_collection.size())
    masked_image_list = ee.List.sequence(0, mask_list.size().subtract(1)).map(lambda i: 
        ee.ImageCollection.fromImages(image_list.slice(ee.Number(i), ee.Number(i).add(1))).first().updateMask(
            ee.ImageCollection.fromImages(mask_list.slice(ee.Number(i), ee.Number(i).add(1))).first()))
    daily_daytime_means = ee.ImageCollection.fromImages(
        ee.List.sequence(0, masked_image_list.size().subtract(1), 24).map(lambda i: 
            ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).mean().set({
                'year': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('year'),
                'month': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('month'),
                'day': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('day'),
                'system:time_start': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('system:time_start')
            })))
    return daily_daytime_means

def aggregate_merra2_nighttime(image_collection, mask_list):
    image_list = image_collection.toList(image_collection.size())
    masked_image_list = ee.List.sequence(0, mask_list.size().subtract(1)).map(lambda i: 
        ee.ImageCollection.fromImages(image_list.slice(ee.Number(i), ee.Number(i).add(1))).first().updateMask(
            ee.ImageCollection.fromImages(mask_list.slice(ee.Number(i), ee.Number(i).add(1))).first()))
    daily_nighttime_means = ee.ImageCollection.fromImages(
        ee.List.sequence(0, masked_image_list.size().subtract(1), 24).map(lambda i: 
            ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).mean().set({
                'year': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('year'),
                'month': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('month'),
                'day': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('day'),
                'system:time_start': ee.ImageCollection.fromImages(masked_image_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('system:time_start')
            })))
    return daily_nighttime_means

def get_daily_tmin(image_collection):
    temp_list = image_collection.select('T10M').toList(image_collection.size())
    daily_tmin = ee.ImageCollection.fromImages(
        ee.List.sequence(0, temp_list.size().subtract(1), 24).map(lambda i: 
            ee.ImageCollection.fromImages(temp_list.slice(ee.Number(i), ee.Number(i).add(24))).min().set({
                'year': ee.ImageCollection.fromImages(temp_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('year'),
                'month': ee.ImageCollection.fromImages(temp_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('month'),
                'day': ee.ImageCollection.fromImages(temp_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('day'),
                'system:time_start': ee.ImageCollection.fromImages(temp_list.slice(ee.Number(i), ee.Number(i).add(24))).first().get('system:time_start')
            })))
    return daily_tmin

def annual_mean_collection(first_year, last_year, image_collection):
    def get_annual_mean(year):
        date_start = ee.Date.fromYMD(year, 1, 1)
        date_end = date_start.advance(1, 'year')
        annual_mean = (ee.ImageCollection(image_collection)
                       .select('T10M')
                       .mean()
                       .set({'year': year, 'system:time_start': date_start}))
        return annual_mean
    years = ee.List.sequence(first_year, last_year)
    mat = ee.ImageCollection(years.map(get_annual_mean))
    return mat

# Aggregate MERRA-2 Data into Daily values
def aggregate_all_merra2(start_date, end_date, day_night_threshold):
    radiation_diagnostics = (ee.ImageCollection('NASA/GSFC/MERRA/rad/2')
        .select(['SWGDN', 'LWGNT', 'SWGNT'])
        .filterDate(str(start_date), str(end_date)))

    single_level_diagnostics = (ee.ImageCollection('NASA/GSFC/MERRA/slv/2')
        .select(['PS', 'QV10M', 'T10M'])
        .filterDate(str(start_date), str(end_date)))

    day_mask_collection = radiation_diagnostics.select('SWGDN').map(day_mask)
    day_mask_list = day_mask_collection.toList(day_mask_collection.size())
    night_mask_collection = radiation_diagnostics.select('SWGDN').map(night_mask)
    night_mask_list = night_mask_collection.toList(night_mask_collection.size())

    rad_daytime = aggregate_merra2_daytime(radiation_diagnostics, day_mask_list).map(lambda img: img.rename(['SWGDN_day', 'LWGNT_day', 'SWGNT_day']))
    rad_nighttime = aggregate_merra2_nighttime(radiation_diagnostics, night_mask_list).map(lambda img: img.rename(['SWGDN_night', 'LWGNT_night', 'SWGNT_night']))
    slv_daytime = aggregate_merra2_daytime(single_level_diagnostics, day_mask_list).map(lambda img: img.rename(['PS_day', 'QV10M_day', 'T10M_day']))
    slv_nighttime = aggregate_merra2_nighttime(single_level_diagnostics, night_mask_list).map(lambda img: img.rename(['PS_night', 'QV10M_night', 'T10M_night']))
    tmin = get_daily_tmin(single_level_diagnostics).map(lambda img: img.rename(['Tmin']))
    annual_temp = annual_mean_collection(start_date.year, end_date.year, single_level_diagnostics).map(lambda img: img.rename(['Annual_Mean_Temp']))
    
    return {'rad_daytime': rad_daytime, 'rad_nighttime': rad_nighttime, 'slv_daytime': slv_daytime, 'slv_nighttime': slv_nighttime, 'tmin': tmin, 'annual_temp': annual_temp}

# MODIS Albedo
#-----------------------------------
def aggregate_albedo(start_date, end_date):
    # Mosaic over a one-week rolling window
    albedo_list = []
    for i in range((end_date - start_date).days):
        image_date = start_date + timedelta(days = i)
        mosaic_start = image_date - timedelta(days = 3)
        mosaic_end = image_date + timedelta(days = 4)
        albedo = (ee.ImageCollection('MODIS/061/MCD43A3')
            .filterDate(str(mosaic_start), str(mosaic_end))
            .select('Albedo_BSA_shortwave')
            .map(roi_clip))
        albedo_mosaic = albedo.mosaic()
        albedo_resampled = albedo_mosaic.resample('bicubic')
        blended_image = albedo_resampled.blend(albedo_mosaic)   # Replace any masked pixels in the mosaic with interpolated data
        albedo_list.append(ee.Image(blended_image).multiply(0.001))
    return ee.ImageCollection(albedo_list)

# Elevation
#-----------------------------------
def aggregate_elevation():
    CA_dem = (ee.ImageCollection('NRCan/CDEM')
        .mosaic()
        .toFloat()
        .reproject(
            crs = ee.ImageCollection('NRCan/CDEM').first().projection(),
            scale = 30))
    AK_dem = (ee.Image('USGS/3DEP/10m')
        .reproject(
            crs = CA_dem.projection(),
            scale = 30))
    dem = roi_clip(ee.ImageCollection([CA_dem, AK_dem]).mosaic().rename(['Elevation']))
    return dem

# Landcover
#-----------------------------------
def aggregate_lc(start_date):
    lc = roi_clip(ee.Image('projects/borealet/assets/Wang_LC'))
    bnames = lc.bandNames().getInfo()
    new_names = []
    for name in bnames:
        new_names.append(f'Landcover_{name[15:19]}')
    lc = lc.select(bnames, new_names).select(f'Landcover_{str(start_date.year)}')
    return lc

# HISTARFM LAI and FPAR
#-----------------------------------
def aggregate_lai(start_date, end_date):
    lai_monthly = ee.ImageCollection('projects/borealet/assets/histarfm_lai')
    lai_img_list = []
    for i in range((end_date - start_date).days):
        image_date = start_date + timedelta(days = i)
        year = image_date.year
        month = image_date.month
        day = image_date.day
        if month < 5:
            lai_img = lai_monthly.filter(ee.Filter.eq('year', year)).filter(ee.Filter.eq('month', 4)).first().set('day', day)
        elif month < 11:
            lai_img = lai_monthly.filter(ee.Filter.eq('year', year)).filter(ee.Filter.eq('month', month)).first().set('day', day)
        elif month > 10:
            lai_img = lai_monthly.filter(ee.Filter.eq('year', year)).filter(ee.Filter.eq('month', 10)).first().set('day', day)
        lai_img_list.append(lai_img)
    lai_filled = ee.ImageCollection(lai_img_list)
    return lai_filled

def aggregate_fpar(start_date, end_date):
    fpar_monthly = ee.ImageCollection('projects/borealet/assets/histarfm_fpar')
    fpar_img_list = []
    for i in range((end_date - start_date).days):
        image_date = start_date + timedelta(days = i)
        year = image_date.year
        month = image_date.month
        day = image_date.day
        if month < 5:
            fpar_img = fpar_monthly.filter(ee.Filter.eq('year', year)).filter(ee.Filter.eq('month', 4)).first().set('day', day)
        elif month < 11:
            fpar_img = fpar_monthly.filter(ee.Filter.eq('year', year)).filter(ee.Filter.eq('month', month)).first().set('day', day)
        elif month > 10:
            fpar_img = fpar_monthly.filter(ee.Filter.eq('year', year)).filter(ee.Filter.eq('month', 10)).first().set('day', day)
        fpar_img_list.append(fpar_img)
    fpar_filled = ee.ImageCollection(fpar_img_list)
    return fpar_filled

# Combine all inputs
#-----------------------------------
def combine_images(annual_temp_image, rad_daytime_image, slv_daytime_image, rad_nighttime_image, 
        slv_nighttime_image, tmin_image, albedo_image, lai_image, fpar_image, dem, lc_image):
    lai = lai_image
    fpar = fpar_image
    elevation = dem.reproject(crs = lai.projection())
    lc = lc_image.rename('Landcover').reproject(crs = lai.projection())
    annual_temp = annual_temp_image.resample('bilinear').reproject(crs = lai.projection(), scale = lai.projection().nominalScale())
    rad_daytime = rad_daytime_image.resample('bilinear').reproject(crs = lai.projection(), scale = lai.projection().nominalScale())
    slv_daytime = slv_daytime_image.resample('bilinear').reproject(crs = lai.projection(), scale = lai.projection().nominalScale())
    rad_nighttime = rad_nighttime_image.resample('bilinear').reproject(crs = lai.projection(), scale = lai.projection().nominalScale())
    slv_nighttime = slv_nighttime_image.resample('bilinear').reproject(crs = lai.projection(), scale = lai.projection().nominalScale())
    tmin = tmin_image.resample('bilinear').reproject(crs = lai.projection(), scale = lai.projection().nominalScale())
    albedo = albedo_image.reproject(crs = lai.projection(), scale = lai.projection().nominalScale()).rename('Albedo')
    
    # Create Parameter Images
    tmin_close = ee.Image(-999).rename('tmin_close')
    tmin_open = ee.Image(-999).rename('tmin_open')
    vpd_open = ee.Image(-999).rename('vpd_open')
    vpd_close = ee.Image(-999).rename('vpd_close')
    gl_sh = ee.Image(-999).rename('gl_sh')
    gl_wv = ee.Image(-999).rename('gl_wv')
    g_cuticular = ee.Image(-999).rename('g_cuticular')
    csl = ee.Image(-999).rename('csl')
    rbl_min = ee.Image(-999).rename('rbl_min')
    rbl_max = ee.Image(-999).rename('rbl_max')
    beta = ee.Image(-999).rename('beta')
    
    # Populate Parameter Images with values based on Landcover
    for i in range(9):
        tmin_close = tmin_close.where(lc.eq(i + 1), ee.Image(params['tmin_close_list'][i]))
        tmin_open = tmin_open.where(lc.eq(i + 1), ee.Image(params['tmin_open_list'][i]))
        vpd_open = vpd_open.where(lc.eq(i + 1), ee.Image(params['vpd_open_list'][i]))
        vpd_close = vpd_close.where(lc.eq(i + 1), ee.Image(params['vpd_close_list'][i]))
        gl_sh = gl_sh.where(lc.eq(i + 1), ee.Image(params['gl_sh_list'][i]))
        gl_wv = gl_wv.where(lc.eq(i + 1), ee.Image(params['gl_wv_list'][i]))
        g_cuticular = g_cuticular.where(lc.eq(i + 1), ee.Image(params['g_cuticular_list'][i]))
        csl = csl.where(lc.eq(i + 1), ee.Image(params['csl_list'][i]))
        rbl_min = rbl_min.where(lc.eq(i + 1), ee.Image(params['rbl_min_list'][i]))
        rbl_max = rbl_max.where(lc.eq(i + 1), ee.Image(params['rbl_max_list'][i]))
        beta = beta.where(lc.eq(i + 1), ee.Image(params['beta_list'][i]))
    
    combined = ee.Image(lai.addBands([fpar, albedo, elevation, lc, annual_temp, rad_daytime, slv_daytime, rad_nighttime, slv_nighttime, tmin,
                tmin_close, tmin_open, vpd_open, vpd_close, gl_sh, gl_wv, g_cuticular, csl, rbl_min, rbl_max, beta])
                # Mask all pixels with missing model parameters
                .updateMask(tmin_close.neq(-999))
                .set({'year': fpar.get('year'), 'month': fpar.get('month'), 'day': fpar.get('day')}))
    return combined.set('system:time_start', rad_daytime.get('system:time_start'))

# Map over the collections using the combine_images function
def aggregate_all_inputs(start_date, end_date):
    annual_temp = aggregate_all_merra2(start_date, end_date, 1.0)['annual_temp'].first()
    rad_daytime = aggregate_all_merra2(start_date, end_date, 1.0)['rad_daytime']
    slv_daytime = aggregate_all_merra2(start_date, end_date, 1.0)['slv_daytime']
    rad_nighttime = aggregate_all_merra2(start_date, end_date, 1.0)['rad_nighttime']
    slv_nighttime = aggregate_all_merra2(start_date, end_date, 1.0)['slv_nighttime']
    tmin = aggregate_all_merra2(start_date, end_date, 1.0)['tmin']
    albedo_filled = aggregate_albedo(start_date, end_date)
    lai_filled = aggregate_lai(start_date, end_date)
    fpar_filled = aggregate_fpar(start_date, end_date)
    dem = aggregate_elevation()
    lc = aggregate_lc(start_date)
    
    all_inputs = ee.ImageCollection(ee.List.sequence(0, lai_filled.size().subtract(1)).map(
        lambda image: combine_images(
            annual_temp,
            ee.Image(rad_daytime.toList(lai_filled.size()).get(image)),
            ee.Image(slv_daytime.toList(lai_filled.size()).get(image)),
            ee.Image(rad_nighttime.toList(lai_filled.size()).get(image)),
            ee.Image(slv_nighttime.toList(lai_filled.size()).get(image)),
            ee.Image(tmin.toList(lai_filled.size()).get(image)),
            ee.Image(albedo_filled.toList(lai_filled.size()).get(image)),
            ee.Image(lai_filled.toList(lai_filled.size()).get(image)),
            ee.Image(fpar_filled.toList(lai_filled.size()).get(image)),
            dem,
            lc)))
    return all_inputs
