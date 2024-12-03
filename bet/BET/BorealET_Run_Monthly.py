'''
Script to run the BorealET model
'''
import datetime
from dateutil.relativedelta import relativedelta
import BorealET_Input_Collection_v1 as bet_input
import BorealET_v1 as bet
import time
import ee
ee.Authenticate()
ee.Initialize(project = 'borealet')

ROI = ee.FeatureCollection('projects/borealet/assets/BorealET_ROI').geometry()

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





# Function to Calculate ET with BorealET
def et_fun(image):
    return bet.evapotranspiration(image).set('system:time_start', image.get('system:time_start'))
    
def scale_and_add_totals(image):
    # Get the system:time_start
    time_start = image.get('system:time_start')
    
    # Divide by the scale factor
    image = image.divide(1e7)
    
    # Add bands for daytime ET, nighttime ET, and total daily ET
    et_day = image.select('e_canopy_day').add(image.select('e_soil_day')).add(image.select('t_day')).rename('et_day')
    et_night = image.select('e_canopy_night').add(image.select('e_soil_night')).add(image.select('t_night')).rename('et_night')
    et_total = et_day.add(et_night).rename('et_total')
    
    return image.addBands([et_day, et_night, et_total]).set('system:time_start', time_start)

# Function to check if any tasks are running
def any_task_running():
    tasks = ee.batch.Task.list()
    return any(task.state in ['READY', 'RUNNING'] for task in tasks)


for year in range(2003, 2009):
    for month in range(4, 11):
        # Define eport period:
        start = datetime.date(year, month, 1)
        end = start + relativedelta(months = 1)
        
        # Skip dates that are already done
        if start < datetime.date(2006, 8, 1):
            pass
        else:
            print('Export start: ', start, '  Export end: ', end)
            
            # Make daily ET images for the period
            all_inputs = bet_input.aggregate_all_inputs(start, end)
            monthly_et_image = (all_inputs
                    .map(et_fun)
                    .map(scale_and_add_totals)
                    .mean()
                    .select(['e_canopy_day', 'e_canopy_night', 'e_soil_day', 'e_soil_night', 't_day', 't_night'])
                    .multiply(1e7).toUint32()
                    .set({'month': start.month, 'year': start.year}))
            
            # Define the export task
            task = ee.batch.Export.image.toAsset(
                image=monthly_et_image,
                description=f'et_{start.year}_{start.month}',
                assetId=f'projects/borealet/assets/BorealET_Monthly_20241121/BorealET_{start.year}_{start.month}',
                region=ROI,
                scale=30,
                maxPixels=1e13
            )
            
            # Start the export task
            task.start()
            print(f'Started export task for Boreal_ET_{start.year}_{start.month}')
            
            # Wait for the current batch of tasks to finish before starting the next one
            minutes = 0
            while any_task_running():
                print(f'{year}-{month} task was submitted {minutes} minutes ago and is still running')
                minutes += 5
                time.sleep(60 * 5)  # Wait for 5 minutes before checking again
