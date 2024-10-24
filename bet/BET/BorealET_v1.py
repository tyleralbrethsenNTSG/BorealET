'''
Boreal ET
The algorithm for terrestrial boreal evapotranspiration.

This model is adapted from the the MOD16 architecture to accept alternate driver data for 
Landcover, FPAR, and LAI to produce a 30-meter evapotranspiration product that is driven
by a variable landcover classification. It was calibrated for Boreal and Tundra ecosystems.
'''

import ee
ee.Authenticate()
ee.Initialize(project = 'ee-tyleralbrethsen')


# Define constants:
STEFAN_BOLTZMANN = 5.67e-8 # Stefan-Boltzmann constant, W m-2 K-4
SPECIFIC_HEAT_CAPACITY_AIR = 1013 # J kg-1 K-1, Monteith & Unsworth (2001)
MOL_WEIGHT_WET_DRY_RATIO_AIR = 0.622 # Ratio of molecular weight of water vapor to that of dry air (ibid.)
TEMP_LAPSE_RATE = 0.0065 # Standard temperature lapse rate [-(deg K) m-1]
GRAV_ACCEL = 9.80665 # Gravitational acceleration [m s-2]
GAS_LAW_CONST = 8.3143 # Ideal gas law constant [m3 Pa (mol)-1 K-1]
AIR_MOL_WEIGHT = 28.9644e-3 # Molecular weight of air [kg (mol)-1]
STD_TEMP_K = 288.15 # Standard temperature at sea level [deg K]
STD_PRESSURE_PASCALS = 101325.0 # Standard pressure at sea level [Pa]
# A pre-determined quantity, not physically meaningful, used in air_pressure()
AIR_PRESSURE_RATE = GRAV_ACCEL / (TEMP_LAPSE_RATE * (GAS_LAW_CONST / AIR_MOL_WEIGHT))


# Latent Heat of Vaporization (J kg-1)
def latent_heat_vaporization(temp_image):
    lhv = ee.Image.constant(2.501).subtract(temp_image.subtract(273.15).multiply(0.002361)).multiply(1e6)
    return lhv.rename('lhv')

# Air Density
def air_density(image, suffix):
    
    temp = image.select(f'T10M{suffix}')
    rh = rhumidity(image, suffix)     # This won't work until you update rhumidity
    pressure = image.select(f'PS{suffix}')
    
    pressure_expression = pressure.divide(100).multiply(0.348444)
    humidity_expression = rh.multiply(100)
    temp_expression = temp.subtract(273.15).multiply(0.00252).subtract(0.020582)
    density = ee.Image(pressure_expression.subtract(humidity_expression.multiply(temp_expression)).divide(temp))
    return density.rename(f'density{suffix}')


# Air Pressure
def air_pressure(image):
    
    elevation = image.select('Elevation')
    AIR_PRESSURE_RATE = GRAV_ACCEL / (TEMP_LAPSE_RATE * (GAS_LAW_CONST / AIR_MOL_WEIGHT))
    
    temp_ratio = ee.Image.constant(1).subtract(elevation.multiply(TEMP_LAPSE_RATE).divide(STD_TEMP_K))
    air_pressure = temp_ratio.pow(AIR_PRESSURE_RATE).multiply(STD_PRESSURE_PASCALS)
    
    return air_pressure.rename('air_pressure')

# Vapor Pressure Deficit
def vapor_pressure_deficit(image, suffix):
    
    temp = image.select(f'T10M{suffix}')
    qv10m = image.select(f'QV10M{suffix}')
    pressure = image.select(f'PS{suffix}')
    
    temp_c = temp.subtract(273.15)
    avp_numerator = qv10m.multiply(pressure)
    avp_denom = ee.Image.constant(0.622).add(ee.Image.constant(0.379).multiply(qv10m))
    avp = avp_numerator.divide(avp_denom)
    svp = ee.Image.constant(610.7).multiply(ee.Image.constant(2.71828).pow(temp_c.multiply(17.38).divide(temp_c.add(239))))
    vpd = svp.subtract(avp)
    return vpd.rename(f'vpd{suffix}')

# Saturation Vapor Pressure
def saturation_vapor_pressure(image, suffix):
    temp = image.select(f'T10M{suffix}')
    svp_constant = 1e3 * 0.6108
    
    temp_c = temp.subtract(273.15)
    svp = ee.Image.constant(svp_constant).multiply(ee.Image.constant(2.71828).pow(temp_c.multiply(17.27).divide(temp_c.add(237.3))))
    return svp.rename(f'svp{suffix}')

# Relative Humidity
def rhumidity(image, suffix):
    
    esat = saturation_vapor_pressure(image, suffix)
    vpd = vapor_pressure_deficit(image, suffix)
    avp = esat.subtract(vpd)

    # Where avp < 0, set avp to zero
    avp = avp.max(0)
    
    rh = avp.divide(esat)
    # Where rh > 1, set rh to 1.
    rh = rh.min(1)
    return rh.rename(f'rh{suffix}')

# Psychrometric Constant
def psychrometric_constant(image, suffix):
    temp = image.select(f'T10M{suffix}')
    lhv = latent_heat_vaporization(temp)
    pressure = air_pressure(image)

    psy_constant = pressure.multiply(SPECIFIC_HEAT_CAPACITY_AIR).divide(lhv.multiply(MOL_WEIGHT_WET_DRY_RATIO_AIR))
    return psy_constant.rename(f'psy{suffix}')

# Slope of the Saturation Vapor Pressure Curve
def svp_slope(image, suffix):
    temp = image.select(f'T10M{suffix}')
    svp = saturation_vapor_pressure(image, suffix)

    numerator = svp.multiply(17.38).multiply(239)
    denom = temp.subtract(273.15).add(239).pow(2)
    svp_slope = numerator.divide(denom)

    return svp_slope.rename(f'svp_slope{suffix}')

# Linear Constraint (from MOD17)
def linear_constraint(min_image, max_image, form: str = None):
    assert form is None or form in ('reversed', 'binary'),\
        'Argument "form" must be None or one of: "reversed", "binary"'
    if form == 'reversed':
        return lambda img: img.max(min_image).min(max_image).subtract(max_image).abs().divide(max_image.subtract(min_image))
    if form == 'binary':
        return lambda img: img.where(img.neq(1), xmin).where(img.eq(1), xmax)
    return lambda img: img.max(min_image).min(max_image).subtract(min_image).divide(max_image.subtract(min_image))

# Surface Conductance
def surface_conductance(image, suffix):
    m_tmin = linear_constraint(image.select('tmin_close'), image.select('tmin_open'))
    m_vpd = linear_constraint(image.select('vpd_open'), image.select('vpd_close'), 'reversed')
    
    m_tmin_image = m_tmin(image.select('Tmin').subtract(273.15))
    m_vpd_image = m_vpd(vapor_pressure_deficit(image, suffix))
    s_conductance = ee.Image(m_tmin_image.multiply(m_vpd_image).multiply(image.select('csl')))
    return s_conductance.rename(f'surface_conductance{suffix}')

# Soil Heat Flux
def soil_heat_flux(image, suffix):

    r'''
    Soil heat flux [MJ m-2], based on surface energy balance. If the mean
    annual temperature is between `Tmin_close` and 25 deg C and the
    contrast between daytime and nighttime air temperatures is greater
    than 5 deg C:

    $$
    G_{\text{soil}} = 4.73 T - 20.87
        \quad\text{iff}\quad T_{\text{min,close}}
        \le T_{\text{annual}} < 25 ,\,
            T_{\text{day}} - T_{\text{night}} \ge 5
    $$

    Otherwise, soil heat flux is zero.

    Finally, soil heat flux should be no greater than 39% of the net
    radiation to the land surface.

    $$
    G_{\text{soil}} = 0.39 A \quad\text{iff}\quad G_{soil} > 0.39 |A|
    $$

    See ca. Line 4355 in MODPR16_main.c and ca. User Guide Eq. 9
    '''
    albedo = image.select('Albedo')
    temp_annual = image.select('Annual_Mean_Temp')
    sw_rad_day = image.select('SWGDN_day')
    temp_day = image.select('T10M_day')
    temp_night = image.select('T10M_night')
    
    temp_annual_mask = temp_annual.lt(273.15 + 25).And(temp_annual.gte(image.select('tmin_close').add(273.15)))
    temp_diff_mask = temp_day.subtract(temp_night).gte(5.0)
    condition = temp_annual_mask.And(temp_diff_mask)
    
    if suffix == '_day':
        temp = image.select('T10M_day')
        lw_net = image.select('LWGNT_day')
        rad_net = sw_rad_day.multiply(ee.Image(1).subtract(albedo)).add(lw_net)   # TALK TO ARTHUR
    
    if suffix == '_night':
        temp = image.select('T10M_night')
        lw_net = image.select('LWGNT_night')
        rad_net = lw_net                                                        # TALK TO ARTHUR
    
    rad_net_masked = rad_net.updateMask(condition)
    temp_masked = temp.updateMask(condition)
    zero_image = temp.multiply(0)
    
    g_masked = temp_masked.subtract(273.15).multiply(4.73).subtract(20.87)
    g = g_masked.unmask(0).where(condition.Not(), zero_image)
    condition2 = g.abs().gt(rad_net.abs().multiply(0.39))
    g = g.where(condition2, rad_net.multiply(0.39))
    
    return g.rename(f'g{suffix}')

# Soil Radiation
def radiation_soil(image, suffix):

    albedo = image.select('Albedo')
    temp_annual = image.select('Annual_Mean_Temp')
    fpar = image.select('FPAR')
    sw_rad_day = image.select('SWGDN_day')
    lw_net_day = image.select('LWGNT_day')
    temp_day = image.select('T10M_day')
    temp_night = image.select('T10M_night')
    
    g_soil = soil_heat_flux(image, suffix)
    rad_net_day = sw_rad_day.multiply(ee.Image(1).subtract(albedo)).add(lw_net_day)      # TALK TO ARTHUR
    rad_net_night = image.select('LWGNT_night')
    
    if suffix == '_day':
        temp = image.select('T10M_day')
        condition = rad_net_day.subtract(g_soil).lt(0).And(rad_net_day.gt(0))
        g_soil = g_soil.where(condition, rad_net_day)
        rad_soil = ee.Image(1).subtract(fpar).multiply(rad_net_day.subtract(g_soil))
    
    if suffix == '_night':
        temp = image.select('T10M_night')
        condition = rad_net_day.gt(0).And(rad_net_night.subtract(g_soil).lt(rad_net_day.multiply(-0.5)))
        g_soil = g_soil.where(condition, rad_net_day.multiply(0.5).add(rad_net_night))
        rad_soil = ee.Image(1).subtract(fpar).multiply(rad_net_night.subtract(g_soil))

    return rad_soil.rename(f'rad_soil{suffix}')


# Evapotranspiration Components

# Evaporation from Wet Canopy
def evaporation_wet_canopy(image, suffix, lhv = ee.Image(None), rh = ee.Image(None), f_wet = ee.Image(None), tiny: float = 1e-7):

    albedo = image.select('Albedo')
    lai = image.select('LAI')
    fpar = image.select('FPAR')
    temp = image.select(f'T10M{suffix}')
    sw_rad = image.select(f'SWGDN{suffix}')
    lw_net = image.select(f'LWGNT{suffix}')
    pressure = air_pressure(image) # Only depends on elevation
    vpd = vapor_pressure_deficit(image, suffix)
    lc = image.select('Landcover')
    
    # From lines 588-592 in MOD16.evapotranspiration()
    rad_net = sw_rad.multiply(ee.Image(1).subtract(albedo)).add(lw_net)
    rad_canopy = rad_net.multiply(fpar)
    
    if lhv == ee.Image(None):
        temp = image.select(f'T10M{suffix}')
        lhv = latent_heat_vaporization(temp)
    if rh == ee.Image(None):
        rh = rhumidity(image, suffix)
    if f_wet == ee.Image(None):
        rh_condition = rh.gte(0.7)
        f_wet = ee.Image(0).where(rh_condition, rh.pow(4))
    
    f_wet = f_wet.where(f_wet.eq(0), f_wet.add(tiny))
    lai = lai.where(lai.eq(0), lai.add(tiny))
    s = svp_slope(image, suffix)
    rho = air_density(image, suffix)

     # Wet canopy resistance to sensible heat ("rhc")
    r_h = ee.Image(1).divide(image.select('gl_sh').multiply(lai).multiply(f_wet))
    # Wet canopy resistance to evaporated water on the surface ("rvc")
    r_e = ee.Image(1).divide(image.select('gl_wv').multiply(lai).multiply(f_wet))
    # Resistance to radiative heat transfer through air ("rrc")
    r_r = rho.multiply(SPECIFIC_HEAT_CAPACITY_AIR).divide(temp.pow(3).multiply(STEFAN_BOLTZMANN).multiply(4))

    # Aerodynamic resistance to evaporated water on the wet canopy surface ("rhrc")
    r_a_wet = r_h.multiply(r_r).divide(r_h.add(r_r)) # (s m-1)
    numer = (f_wet.multiply(s.multiply(rad_canopy)
                                 .add(rho.multiply(SPECIFIC_HEAT_CAPACITY_AIR)
                                      .multiply(fpar)
                                      .multiply(vpd)
                                      .multiply(ee.Image(1).divide(r_a_wet)))))
    denom = (s.add(pressure
                   .multiply(SPECIFIC_HEAT_CAPACITY_AIR)
                   .multiply(r_e)
                   .multiply(ee.Image(1)
                   .divide(lhv
                           .multiply(MOL_WEIGHT_WET_DRY_RATIO_AIR)
                           .multiply(r_a_wet)))))

    # Mu et al. (2011), Equation 17; PET (J m-2 s-1) is divided by the
    #   latent heat of vaporization (J kg-1) to obtain mass flux (kg m-2 s-1)
    numer = numer.max(0)
    evap = numer.divide(denom).divide(lhv)
    # If f_wet or LAI are ~zero, then wet canopy evaporation is zero
    tiny_condition = f_wet.lte(tiny).Or(lai.lte(tiny))
    evap = evap.where(tiny_condition, 0)
    # For bare soil pixels, canopy evaporation is zero
    evap = evap.where(lc.eq(1), 0)
    # None of the components of total ET can be negative:
    evap = evap.max(0)
    return evap.rename(f'e_canopy{suffix}')

# Evaporation from Bare Soil
def evaporation_soil(image, suffix, lhv = ee.Image(None), rh = ee.Image(None), f_wet = ee.Image(None)):

    fpar = image.select('FPAR')
    temp = image.select(f'T10M{suffix}')
    pressure = air_pressure(image)          # Only depends on elevation
    vpd = vapor_pressure_deficit(image, suffix)
    rad_soil = radiation_soil(image, suffix)
    
    # Correction for atmospheric temperature and pressure
    r_corr = ee.Image(101300).divide(pressure).multiply(temp.divide(293.15).pow(1.75))

    if lhv == ee.Image(None):
        temp = image.select(f'T10M{suffix}')
        lhv = latent_heat_vaporization(temp)
    if rh == ee.Image(None):
        rh = rhumidity(image, suffix)
    if f_wet == ee.Image(None):
        rh_condition = rh.gte(0.7)
        f_wet = ee.Image(0).where(rh_condition, rh.pow(4))

    s = svp_slope(image, suffix)
    rho = air_density(image, suffix)
    gamma = psychrometric_constant(image, suffix)

    # Resistance to radiative heat transfer through air ("rrc")
    r_r = rho.multiply(SPECIFIC_HEAT_CAPACITY_AIR).divide(temp.pow(3).multiply(STEFAN_BOLTZMANN).multiply(4))

    r_tot = (image.select('rbl_max')
             .subtract(image.select('rbl_max')
                       .subtract(image.select('rbl_min'))
                       .multiply(image.select('vpd_close')
                                 .subtract(vpd))
                       .divide(image.select('vpd_close')
                               .subtract(image.select('vpd_open')))))
    r_tot = r_tot.where(vpd.lte(image.select('vpd_open')), image.select('rbl_min')).where(vpd.gte(image.select('vpd_close')), image.select('rbl_max'))
    r_tot = r_tot.divide(r_corr)
    r_as = r_tot.multiply(r_r).divide(r_tot.add(r_r))

    numer = (s.multiply(rad_soil)
             .add(rho
                  .multiply(SPECIFIC_HEAT_CAPACITY_AIR)
                  .multiply(ee.Image(1)
                            .subtract(fpar))
                  .multiply(vpd
                            .divide(r_as))))
    denom = (s.add(gamma
                       .multiply(r_tot
                       .divide(r_as))))
    evap_sat = (numer.multiply(f_wet).divide(denom))
    evap_unsat = (numer.multiply(ee.Image(1).subtract(f_wet)).divide(denom))
    e = evap_sat.max(0)
    e = (e.where(evap_unsat.gte(0), e.add(evap_unsat
                                          .multiply(rh
                                                   .pow(vpd.divide(image.select('beta')))))))
    e_soil = e.divide(lhv)
    # None of the components of total ET can be negative:
    e_soil = e_soil.max(0)
    return e_soil.rename(f'e_soil{suffix}')

# Transpiration
def transpiration(image, suffix, lhv = ee.Image(None), rh = ee.Image(None), f_wet = ee.Image(None), tiny: float = 1e-7):
    
    albedo = image.select('Albedo')
    lai = image.select('LAI')
    fpar = image.select('FPAR')
    tmin = image.select('Tmin')
    temp = image.select(f'T10M{suffix}')
    sw_rad = image.select(f'SWGDN{suffix}')
    lw_net = image.select(f'LWGNT{suffix}')
    pressure = air_pressure(image)          # Only depends on elevation
    vpd = vapor_pressure_deficit(image, suffix)
    lc = image.select('Landcover')
    
    # From lines 588-592 in MOD16.evapotranspiration()
    rad_net = sw_rad.multiply(ee.Image(1).subtract(albedo)).add(lw_net)
    rad_canopy = rad_net.multiply(fpar)
    
    # Correction for atmospheric temperature and pressure
    r_corr = ee.Image(101300).divide(pressure).multiply(temp.divide(293.15).pow(1.75))
    if lhv == ee.Image(None):
        lhv = latent_heat_vaporization(temp)
    if rh == ee.Image(None):
        rh = rhumidity(image, suffix)
    if f_wet == ee.Image(None):
        rh_condition = rh.gte(0.7)
        f_wet = ee.Image(0).where(rh_condition, rh.pow(4))

    s = svp_slope(image, suffix)
    rho = air_density(image, suffix)
    gamma = psychrometric_constant(image, suffix)
    # Resistance to radiative heat transfer through air ("rrc")
    r_r = rho.multiply(SPECIFIC_HEAT_CAPACITY_AIR).divide(temp.pow(3).multiply(STEFAN_BOLTZMANN).multiply(4))

    # Surface conductance (assumed to be zero at nighttime)
    g_surf = ee.Image(0) # (Equation 13, MOD16 C6.1 User's Guide)
    if suffix == '_day':
        g_surf = surface_conductance(image, suffix).divide(r_corr)

    g_cuticular = image.select('g_cuticular').divide(r_corr)
    # -- Canopy conductance, should be zero when LAI or f_wet are zero;
    #   updated calculation for conductance to sensible heat, see User
    #   Guide's Equation 15
    gl_sh = image.select('gl_sh').multiply(lai).multiply(ee.Image(1).subtract(f_wet))
    g = (gl_sh.multiply(g_surf.add(g_cuticular))
         .divide(gl_sh.add(g_surf).add(g_cuticular)))
    g_canopy = (ee.Image(tiny)
                .where(lai.gt(0).And(ee.Image(1).subtract(f_wet).gt(0)), g))
    # -- Aerodynamic resistance to heat, water vapor from dry canopy
    #   surface into the air (Equation 16, MOD16 C6.1 User's Guide)
    r_a_dry = (ee.Image(1).divide(image.select('gl_sh')).multiply(r_r)
                .divide(ee.Image(1).divide(image.select('gl_sh')).add(r_r))) # (s m-1)

    # If canopy radiation is negative, drop (s * A_c) term in the
    #   transpriration calculation (L4046 in MODPR16_main.c)
    
    rad_canopy = rad_canopy.max(0)
    # PLANT TRANSPIRATION
    t = ((ee.Image(1).subtract(f_wet)).multiply(s.multiply(rad_canopy)
         .add(rho.multiply(SPECIFIC_HEAT_CAPACITY_AIR).multiply(fpar).multiply(vpd.divide(r_a_dry)))))
    t = (t.divide(s.add(gamma.multiply(ee.Image(1).add(ee.Image(1).divide(g_canopy).divide(r_a_dry))))))

    # Divide by the latent heat of vaporization (J kg-1) to obtain mass
    #   flux (kg m-2 s-1)
    t = t.divide(lhv).where(g_canopy.lte(tiny), 0)
    
    # For bare soil pixels, t is zero
    t = t.where(lc.eq(1), 0)
    # None of the components of total ET can be negative:
    t = t.max(0)
    return t.rename(f't{suffix}')


# Instantaneous Evapotranspiration
def evapotranspiration(image, lhv = ee.Image(None), rh = ee.Image(None), f_wet = ee.Image(None)):
    r'''
    Instaneous evapotranspiration (ET) [kg m-2 day-1], the sum of bare soil
    evaporation, wet canopy evaporation, and canopy transpiration.
    
    Parameter 'image' must be an ee.Image with the following bands:
    For details on input data prep, see Boreal_ET_Input_Collection_20240710.py
    ---------------------------------------------------------------
    The following bands represent model drivers:
        'Landcover': Annual lancover classification for the ABoVE study region
            From Wang et al (in prep)
        'Annual_Mean_Temp': Annual mean temperature (deg K)
            From MERRA-2 single-level diagnostics
        'SWGDN_day': Daytime surface incoming shortwave flux (W m-2)
            From MERRA-2 radiation diagnostics
        'LWGNT_day': Daytime surface net downward longwave flux (W m-2)
            From MERRA-2 radiation diagnostics
        'SWGNT_day': Daytime surface net downward shortwave flux (W m-2)
            From MERRA-2 radiation diagnostics
        'PS_day': Daytime surface pressure (Pa)
            From MERRA-2 single-level diagnostics
        'QV10M_day': Daytime specific humidity at 10m (mass fraction)
            From MERRA-2 single-level diagnostics
        'T10M_day': Daytime mean temperature (deg K)
            From MERRA-2 single-level diagnostics
        'SWGDN_night': Nighttime surface incoming shortwave flux (W m-2)
            From MERRA-2 radiation diagnostics
        'LWGNT_night': Nighttime surface net downward longwave flux (W m-2)
            From MERRA-2 radiation diagnostics
        'SWGNT_night': Nighttime surface net downward shortwave flux (W m-2)
            From MERRA-2 radiation diagnostics
        'PS_night': Nighttime surface pressure (Pa)
            From MERRA-2 single-level diagnostics
        'QV10M_night': Nighttime specific humidity at 10m (mass fraction)
            From MERRA-2 single-level diagnostics
        'T10M_night': Nighttime mean temperature (deg K)
            From MERRA-2 single-level diagnostics
        'Tmin': Daily minimum temperature (deg K)
            From MERRA-2 single-level diagnostics
        'Albedo': Black-sky albedo for shortwave broadband
            From MCD43A3 (Albedo_BSA_shortwave)
        'LAI': Leaf Area Index
            From HISTARFM
        'FPAR': Fraction of Photosynthetically Active Radiation
            From HISTARFM
        'Elevation': Elevation (m)
            From Natural Resources Canada and USGS
    
    The following bands represent the required model parameters.
    Values are determined in calibration and depend on landcover classification:
        'tmin_close': Temperature at which stomata are almost completely
            closed due to (minimum) temperature stress (deg C)
        'tmin_open': Temperature at which stomata are completely open, i.e.,
            there is no effect of temperature on transpiration (deg C)
        'vpd_open': The VPD at which stomata are completely open, i.e.,
            there is no effect of water stress on transpiration (Pa)
        'vpd_close': The VPD at which stomata are almost completely closed
            due to water stress (Pa)
        'gl_sh': Leaf conductance to sensible heat per unit LAI
            (m s-1 LAI-1)
        'gl_wv'Leaf conductance to evaporated water per unit LAI
            (m s-1 LAI-1)
        'g_cuticular': Leaf cuticular conductance (m s-1)
        'csl': Mean potential stomatal conductance per unit leaf area (m s-1)
        'rbl_min': Minimum atmospheric boundary layer resistance (s m-1)
        'rbl_max': Maximum atmospheric boundary layer resistance (s m-1)
        'beta': Factor in soil moisture constraint on potential soil
            evaporation, i.e., (VPD / beta); from Bouchet (1963)
    
    All bands are reprojected to EPSG:4326 (WGS 84) with a 30-meter resolution.
    '''
    
    et_img = ee.Image(0)
    for suffix in ['_day', '_night']:
        e_canopy = evaporation_wet_canopy(image, suffix, lhv, rh, f_wet)
        e_soil = evaporation_soil(image, suffix, lhv, rh, f_wet)
        t = transpiration(image, suffix, lhv, rh, f_wet)
        et = e_canopy.add(e_soil).add(t).rename(f'et{suffix}')
        et_img = et_img.addBands([et, e_canopy, e_soil, t])
    # Export all bands
    #et_img = et_img.select(['e_canopy_day', 'e_canopy_night', 'e_soil_day', 'e_soil_night', 't_day', 't_night', 'et_day', 'et_night'])
    # Export only components
    et_img = et_img.select(['e_canopy_day', 'e_canopy_night', 'e_soil_day', 'e_soil_night', 't_day', 't_night'])
    # Export only totals
    #et_img = et_img.select(['et_day', 'et_night'])
    # Convert from kg m-2 s-1 to kg m-2 day-1
    et_img = et_img.multiply(60).multiply(60).multiply(24)
    # Apply scale factor to store as an int32
    et_img = et_img.multiply(1e7).toUint32()
    return et_img