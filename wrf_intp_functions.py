### wrf_intp_functions
## Megan Schiede
from wrf import (slp)
def preprocessed(ds):
    ##Save Memory
    droplst=['C1H','C2H','C1F','C2F','C3H','C4H','C3F','C4F','PCB','PC','LU_INDEX','VAR_SSO','HFX_FORCE',
            'LH_FORCE','TSK_FORCE', 'HFX_FORCE_TEND','LH_FORCE_TEND','TSK_FORCE_TEND','NEST_POS','FNM','FNP','RDNW','RDN'] #Variables to drop
    ds.drop(droplst)
    ## Add Height
    ph=ds['PH']
    phb=ds['PHB']
    height=(ph+phb)/9.81
    ds['HEIGHT']=height
    ## Add Pressure
    p=ds['P']
    pb=ds['PB']
    P=p+pb
    ds['PRESSURE']=P
    return ds

#Temperature on isobaric field
def temp_isobaric(ds,pressure_level,label):
    '''
    ds = wrf dataset
    pressure_level = desired pressure level
    '''
    t = ds['T'].values
    g = 9.81 

    ph = ds['PH'] # Perturbation geopotential
    phb = ds['PHB']  # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt #).values
    ds['AGL']=height_agl


    theta_p = ds['T']  # Perturbation potential temperature
    p = ds['P']        # Perturbation pressure
    pb = ds['PB']      # Base-state pressure

    p0 = 1000.0 * 100  # Reference pressure in Pa (hPa to Pa)
    R = 287.0          # Specific gas constant for dry air
    Cp = 1004.0        # Specific heat at constant pressure
    kappa = R / Cp
    P = p + pb

    theta = theta_p + 300.0  # Assume base-state theta is 300 K

    T = (theta * ((P / p0) ** kappa)).values
    t = T
    p_at_time=ds.PRESSURE

    intp_sel=interplevel(t, p_at_time,pressure_level*100)
    intp_sel=intp_sel.rename({'dim_1':'south_north','dim_2':'west_east'})
    if label == True:
        return intp_sel,pressure_level
    else:
        return intp_sel

def u_isobaric(ds,pressure_level,label):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE 
    us=ds.U.isel(west_east_stag=slice(0, 264))
    intp = interplevel(us, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    if label == True:
        return intp,pressure_level
    else:
        return intp
def v_isobaric(ds,pressure_level, label):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE
    vs=ds.V.isel(south_north_stag=slice(0, 222))
    intp = interplevel(vs, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    if label == True:
        return intp,pressure_level
    else:
        return intp

def ptemp_isobaric(ds,pressure_level,label):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE
    T=ds.T+300 # Assuming theta_0 = 300 K
    intp = interplevel(T, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    if label == True:
        return intp,pressure_level
    else:
        return intp

## to get to mslp
def mslp(ds):
    p0 = 1000.0 * 100  # Reference pressure in Pa (hPa to Pa)
    R = 287.0          # Specific gas constant for dry air
    Cp = 1004.0        # Specific heat at constant pressure
    kappa = R / Cp
    theta_p = ds.variables['T'][:]
    p = ds.variables['P'][:]        # Perturbation pressure
    pb = ds.variables['PB'][:]
    P = p + pb
    theta = theta_p + 300.0  # Assume base-state theta is 300 K
    
    T = (theta * ((P / p0) ** kappa)) 
    g = 9.81 

    ph = ds.variables['PH'][:]  # Perturbation geopotential
    phb = ds.variables['PHB'][:]  # Base-state geopotential
    hgt = ds.variables['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters
   
    z=z[1:]
    SLP=slp(z.values,T.values,P.values,ds.QVAPOR.values) ## Use slp fn from WRF package
    SLP=SLP.rename({'dim_0':'lon','dim_1':'lat'})
    return(SLP)

def sh_isobaric(ds,pressure_level,label):    
    g = 9.81 

    ph = ds['PH']  # Perturbation geopotential
    phb = ds['PHB']  # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt 
    ds['AGL']=height_agl


    theta_p = ds['T']  # Perturbation potential temperature
    p = ds['P']      # Perturbation pressure
    pb = ds.['PB']      # Base-state pressure

    p0 = 1000.0 * 100  # Reference pressure in Pa (hPa to Pa)
    R = 287.0          # Specific gas constant for dry air
    Cp = 1004.0        # Specific heat at constant pressure
    kappa = R / Cp
    P = p + pb

    theta = theta_p + 300.0  # Assume base-state theta is 300 K

    T = (theta * ((P / p0) ** kappa)).values

    qvapor = ds['QVAPOR'] # Water vapor mixing ratio
    q = qvapor / (1 + qvapor)
    qvapor = q  # Specific humidity (kg/kg), Water vapor mixing ratio the same as SH to a few percent
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE
    intp = interplevel(qvapor, p_at_time, pressure_levels)
    if label == True:
        return intp,pressure_level
    else:
        return intp

def rh_isobaric(ds,pressure_level,label):    
    g = 9.81 

    ph = ds['PH']  # Perturbation geopotential
    phb = ds['PHB'] # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt #).values
    ds['AGL']=height_agl


    theta_p = ds['T']  # Perturbation potential temperature
    p = ds['P']       # Perturbation pressure
    pb = ds['PB']     # Base-state pressure

    p0 = 1000.0 * 100  # Reference pressure in Pa (hPa to Pa)
    R = 287.0          # Specific gas constant for dry air
    Cp = 1004.0        # Specific heat at constant pressure
    kappa = R / Cp
    P = p + pb

    theta = theta_p + 300.0  # Assume base-state theta is 300 K

    T = (theta * ((P / p0) ** kappa)).values

    qvapor = ds['QVAPOR']  # Specific humidity (kg/kg)
    
    # Calculate saturation vapor pressure (E_s) in hPa
    es = 6.112 * np.exp(17.67 * (T-273.15) / ((T-273.15) + 243.5))

    # Calculate actual vapor pressure (E) in hPa
    # Using approximation: E = q * P / (0.622 + 0.378 * q)
    e = qvapor * P / (0.622 + 0.378 * qvapor) / 100  # Convert Pa to hPa

    # Calculate relative humidity
    rh = (e / es) * 100

    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE
    intp = interplevel(rh, p_at_time, pressure_levels)
    
    if label == True:
        return intp,pressure_level
    else:
        return intp

def qrain_isobaric(ds,pressure_level):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE #[0]
    intp = interplevel(ds.QRAIN, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    return intp
def qcloud_isobaric(ds,pressure_level):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE #[0]
    intp = interplevel(ds.QCLOUD, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    return intp
def qsnow_isobaric(ds,pressure_level):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE #[0]
    intp = interplevel(ds.QSNOW, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    return intp
def qice_isobaric(ds,pressure_level):
    pressure_levels=pressure_level*100
    p_at_time=ds.PRESSURE #[0]
    intp = interplevel(ds.QICE, p_at_time, pressure_levels)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    return intp

##########
#Height#
##########

def temp_height(ds,height_level,label):
    g = 9.81 

    ph = ds['PH']  # Perturbation geopotential
    phb = ds['PHB'] # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt
    ds['AGL']=height_agl


    theta_p = ds['T'] # Perturbation potential temperature
    p = ds['P']        # Perturbation pressure
    pb = ds['PB']     # Base-state pressure

    p0 = 1000.0 * 100  # Reference pressure in Pa (hPa to Pa)
    R = 287.0          # Specific gas constant for dry air
    Cp = 1004.0        # Specific heat at constant pressure
    kappa = R / Cp
    P = p + pb

    theta = theta_p + 300.0  # Assume base-state theta is 300 K

    T = (theta * ((P / p0) ** kappa)).values
    
    h_at_time=ds.AGL[1::]
    
    intp_sel=interplevel(T, h_at_time,height_level*1000)
    intp_sel=intp_sel.rename({'dim_1':'south_north','dim_2':'west_east'})
    if label == True:
        return intp_sel,height_level
    else:
        return intp_sel
    
def pressure_height(ds,height_level,label):
    g = 9.81 

    ph = ds['PH'] # Perturbation geopotential
    phb = ds['PHB'] # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt 
    ds['AGL']=height_agl

    p = ds['P']     # Perturbation pressure
    pb = ds['PB']      # Base-state pressure

    P = p + pb
    
    h_at_time=ds.AGL[1::]
    intp_sel=interplevel(P.values, h_at_time,height_level*1000)
    intp=intp_sel.rename({'dim_1':'south_north','dim_2':'west_east'})
    if label == True:
        return intp,height_level
    else:
        return intp

def u_height(ds,height_level,label):
    g = 9.81 

    ph = ds['PH'] # Perturbation geopotential
    phb = ds['PHB']  # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt 
    ds['AGL']=height_agl
    height_level=height_level*1000
    h_at_time=ds.AGL[1::] 
    us=ds.U.isel(west_east_stag=slice(0, 264))
    intp = interplevel(us, h_at_time, height_level) 
    if label == True:
        return intp,height_level
    else:
        return intp
def v_height(ds,height_level,label):
    g = 9.81 

    ph = ds['PH'] # Perturbation geopotential
    phb = ds['PHB']  # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt 
    ds['AGL']=height_agl
    height_level=height_level*1000
    h_at_time=ds.AGL[1::] 
    vs=ds.V.isel(south_north_stag=slice(0, 222))
    intp = interplevel(vs, h_at_time, height_level)
    if label == True:
        return intp,height_level
    else:
        return intp

def sh_height(ds,height_level, label):    
    g = 9.81 

    ph = ds['PH']  # Perturbation geopotential
    phb = ds['PHB']  # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt 
    ds['AGL']=height_agl
    
    qvapor = ds['QVAPOR'][:]  # Specific humidity (kg/kg)
   
    h_at_time=ds.AGL[1::]
    intp = interplevel(qvapor, h_at_time, height_level*1000)
    if label == True:
        return intp,height_level
    else:
        return intp

def rh_height(ds,height_level,label):    
    g = 9.81 

    ph = ds['PH']  # Perturbation geopotential
    phb = ds['PHB'] # Base-state geopotential
    hgt = ds['HGT']  # Terrain height

    z = (ph + phb) / g  # Total geopotential height in meters

    height_agl = z - hgt
    ds['AGL']=height_agl


    theta_p = ds['T']  # Perturbation potential temperature
    p = ds['P']       # Perturbation pressure
    pb = ds['PB']      # Base-state pressure

    p0 = 1000.0 * 100  # Reference pressure in Pa (hPa to Pa)
    R = 287.0          # Specific gas constant for dry air
    Cp = 1004.0        # Specific heat at constant pressure
    kappa = R / Cp
    P = p + pb

    theta = theta_p + 300.0  # Assume base-state theta is 300 K

    T = (theta * ((P / p0) ** kappa)).values

    qvapor = ds['QVAPOR'] # Water vapor mixing ratio

    q = qvapor / (1 + qvapor)


    qvapor = q  # Specific humidity (kg/kg), Water vapor mixing ratio the same as SH to a few percent
    # Calculate saturation vapor pressure (E_s) in hPa
    es = 6.112 * np.exp(17.67 * (T-273.15) / ((T-273.15) + 243.5))

    # Calculate actual vapor pressure (E) in hPa
    # Using approximation: E = q * P / (0.622 + 0.378 * q)
    e = qvapor * P / (0.622 + 0.378 * qvapor) / 100  # Convert Pa to hPa

    # Calculate relative humidity
    rh = (e / es) * 100
    h_at_time=ds.AGL[1::]
    intp = interplevel(rh, h_at_time, height_level*1000)
    if label == True:
        return intp,height_level
    else:
        return intp

def ref_height(ds,height):
    height2=height*1000
    h_at_time=ds.HEIGHT #[0]
    intp = interplevel(ds.REFL_10CM, h_at_time[1::], height2)  #interplevel(variable, pressureindexedatonetime, pressure_desired)
    return intp, height