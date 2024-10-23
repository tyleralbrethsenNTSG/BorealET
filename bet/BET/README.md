BorealET Evapotranspiration Model
=================================

The BorealET terrestrial evapotranspiration model is a Google Earth Engine implementation of MOD16 that utilizes gap-filled Lansat imagery (Moreno-Martínez, 2020) and annual landcover classification (Wang, 2019) and is calibrated on ground data from North American boreal and arctic flux tower sites.
The model is designed to address questions related to the effects of landcover change, disturbance, and vegetation shifts on evapotranspiration.

[Details on the MOD16 algorithm can be found here.](https://arthur-e.github.io/MOD16/)


BorealET Documentation
----------------------

BorealET is run in the Python API for Google Earth Engine. A forward run is conducted by first assembling an input collection, then running the model. Sample scripts are available to run BorealET at daily and monthly time-steps.

BorealET returns an ee.Image at daily time-steps and 30-meter spatial resolution with the following quantities:

- Daytime Evaporation from bare soil surfaces
- Daytime Evaporation from wet canopy surfaces
- Daytime Transpiration from terrestrial vegetation
- Nighttime Evaporation from bare soil surfaces
- Nighttime Evaporation from wet canopy surfaces
- Nighttime Transpiration from terrestrial vegetation

The sum of these six quatities is daily actual evapotranspiration. Values are scaled to convert to int32. Google Earth Engine helper functions (Javascript) are provided to rescale and add a band for daily total evapotranspiration.



References
----------

- Monteith, J. L., and M. Unsworth. 2001. Principles of Environmental Physics. Second Ed.
- Moreno-Martínez, Á., Izquierdo-Verdiguier, E., Maneta, M. P., Camps-Valls, G., Robinson, N., Muñoz-Marí, J., Sedano, F., Clinton, N., & Running, S. W. 2020. Multispectral high resolution sensor fusion for smoothing and gap-filling in the cloud. *Remote Sensing of Environment* 247, 111901. https://doi.org/10.1016/j.rse.2020.111901
- Mu, Q., Heinsch, F. A., Zhao, M., & Running, S. W. 2007. Development of a global evapotranspiration algorithm based on MODIS and global meteorology data. *Remote Sensing of Environment,* 111(4), 519–536.
- Mu, Q., M. Zhao, and S. W. Running. 2011. Improvements to a MODIS global terrestrial evapotranspiration algorithm. *Remote Sensing of Environment* 115 (8):1781–1800.
- Wang, J.A., D. Sulla-Menashe, C.E. Woodcock, O. Sonnentag, R.F. Keeling, and M.A. Friedl. 2019. ABoVE: Landsat-derived Annual Dominant Land Cover Across ABoVE Core Domain, 1984-2014. *ORNL DAAC* Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/1691
