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

# ROI must be the same or within the ROI defined in the bet_input script.
ROI = ee.FeatureCollection('projects/borealet/assets/taylor_bounds').geometry()

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

# Function to check if any tasks are running
def any_task_running():
    tasks = ee.batch.Task.list()
    return any(task.state in ['READY', 'RUNNING'] for task in tasks)


for year in range(2001, 2009):
    for month in range(1, 13):
        # Define eport period:
        start = datetime.date(year, month, 1)
        end = start + relativedelta(months = 1)
        
        # Skip dates that are already done
        if start < datetime.date(2001, 1, 1):
            pass
        else:
            print('Export start: ', start, '  Export end: ', end)
            
            # Make daily ET images for the period
            all_inputs = bet_input.aggregate_all_inputs(start, end)
            et = all_inputs.map(et_fun)
            
            # Export the ImageCollection
            et_list = et.toList(et.size().getInfo())
            for day in range(et_list.size().getInfo()):
                et_img = ee.Image(et_list.get(day))
                date = ee.Date(et_img.get('system:time_start').getInfo()).format('YYYYMMdd').getInfo()
                
                # Define the export task
                task = ee.batch.Export.image.toAsset(
                    image=et_img,
                    description=f'et_{date}',
                    assetId=f'projects/borealet/assets/BorealET_Taylor_v2/BorealET_{date}',
                    region=ROI,
                    scale=30,
                    maxPixels=1e13
                )
                
                # Start the export task
                task.start()
                print(f'Started export task for Boreal_ET_{date}')
            
            # Wait for the current batch of tasks to finish before starting the next one
            minutes = 0
            while any_task_running():
                print(f'{year}-{month} tasks were submitted {minutes} minutes ago and are still running')
                minutes += 5
                time.sleep(60 * 5)  # Wait for 5 minutes before checking again





'''
start = datetime.date(2001, 11, 25)
end = datetime.date(2002, 1, 1)
print('start: ', start, '  end: ', end)


all_inputs = bet_input.aggregate_all_inputs(start, end)
print('Number of input images: ', all_inputs.size().getInfo())

et = all_inputs.map(et_fun)
print('Number of ET images: ', et.size().getInfo())


# Export the ImageCollection
et_list = et.toList(et.size().getInfo())
for i in range(et_list.size().getInfo()):
    et_img = ee.Image(et_list.get(i))
    date = ee.Date(et_img.get('system:time_start').getInfo()).format('YYYYMMdd').getInfo()
    
    # Define the export task
    task = ee.batch.Export.image.toAsset(
        image=et_img,
        description=f'et_{date}',
        assetId=f'projects/borealet/assets/BorealET_20240813/BorealET_{date}',
        region=ROI,
        scale=30,
        maxPixels=1e13
    )
    
    # Start the export task
    task.start()
    print(f'Started export task for Boreal_ET_{date}')
'''










'''
# Test Inputs
#-----------------------------------
start = datetime.date(2003, 1, 1)
end = start + relativedelta(months = 1)

merra2 = bet_input.aggregate_all_merra2(start, end, 1.0)
print(merra2.keys())
print(merra2['rad_daytime'].size().getInfo())

albedo = bet_input.aggregate_albedo(start, end)
print(albedo.size().getInfo())

lc = bet_input.aggregate_lc(start)
print(lc.bandNames().getInfo())

lai = bet_input.aggregate_lai(start, end)
print(lai.size().getInfo())



# For testing:
test_img = all_inputs.first()
# Define the export task
input_img_date = ee.Date(test_img.get('system:time_start').getInfo()).format('YYYYMMdd').getInfo()
task = ee.batch.Export.image.toAsset(
    image=test_img,
    description=f'TEST_INPUT_{input_img_date}',
    assetId=f'projects/borealet/assets/TEST_INPUT_{input_img_date}',
    region=ROI,
    scale=30,
    maxPixels=1e13)
# Start the export task
task.start()
print(f'Started export task for BorealET_TEST_INPUT_{input_img_date}')

# Calculate ET
#-----------------------------------
def et_fun(image):
    return bet.evapotranspiration(image).set('system:time_start', image.get('system:time_start'))
et = all_inputs.map(et_fun)
print(et.size().getInfo())
print(et.first().bandNames().getInfo())
task2 = ee.batch.Export.image.toAsset(
    image=et.first(),
    description=f'TEST_ET_IMG_{input_img_date}',
    assetId=f'projects/borealet/assets/TEST_ET_IMG_{input_img_date}',
    region=ROI,
    scale=30,
    maxPixels=1e13)
# Start the export task
task2.start()
print(f'Started export task for TEST_ET_IMG_{input_img_date}')
'''


