# -*- coding: utf-8 -*-
"""
Tools for calculating open-water evaporation using the aerodynmaic mass-transfer approach.
"""

import cmath as cm
import numpy as np
import math as m
import pandas as pd
import multiprocessing as mp


class Aero(object):
    """
    Manages meterological time series input/output for aerodynamic
    mass-transfer evaporation calculation and contains methods for batch and
    single calculations.  
    
    An :obj:`Aero` object allows the aerodynamic mass-transfer evaporation
    estimation to be calculated from meterological data that is stored in a
    :obj:`pandas.DataFrame` with a date or datetime-like index. The
    :attr:`Aero.df` can be assigned on initialization or later, it can also be
    reassigned at anytime. 

    The :meth:`Aero.single_calc` static method calculates evaporation for a 
    single measurement set and can be used without creating an :obj:`Aero` 
    object, e.g. in another module. For calculating evaporation for a time 
    series of input meterological data use the :meth:`Aero.run` method which
    uses multiple processors (if they are available).
    """

    def __init__(self, df=None):
        if df is not None and not isinstance(df, pd.DataFrame):
            raise TypeError("Must assign a pandas.DataFrame object")
        self._df = df


    def run(self, sensor_height, timestep, variable_names=None, nproc=None):
        """
        Run aerodynamic mass-transfer evaporation routine on time series data
        that contains necessary input in-place and in parallel. 

        Arguments:
            sensor_height (float): height of sensor in meters.
            timestep (float or int): sensor sampling frequency in seconds.

        Keyword Arguments:
            variable_names (None or dict): default None. Dictionary with user 
                variable names as keys and variable names needed for
                :mod:`aeroevap` as values. If None, the needed input variables
                must be named correctly in the :attr:`Aero.df` dataframe: 'WS',
                'P', 'T_air', 'T_skin', and 'RH' for windspeed, air pressure,
                air temperature, skin temperature, and relative humidity
                resepctively.
            nproc (None or int): default None. If none use half of the available
                cores for parallel calculations. 

        Returns:
            None

        Hint:
            A :obj:`pandas.DataFrame` must be assigned to the :attr:`Aero.df`
            instance property before calling :meth:`Aero.run`. If the names of
            the required meterological variables in the dataframe are not
            named correctly you may pass a dictionary to the ``variable_names``
            argument which maps your names to those used by ``AeroEvap``. For
            example if your surface temperature column is named 'surface_temp'
            then 
            
            >>> variable_names = {'surface_temp' : 'T_skin'}


           
        """

        if not isinstance(self._df, pd.DataFrame):
            print(
                'ERROR: no pandas.DataFrame assigned to Aero.df, please '
               'assign first.'
            )
            return

        if variable_names is not None:
            df = self._df.rename(columns=variable_names)
        else:
            df = self._df

        df['date'] = df.index
        df['SH'] = sensor_height
        df['dt'] = timestep
        input_vars = ['date', 'WS', 'P', 'T_air', 'T_skin', 'RH', 'SH', 'dt']
        if not set(input_vars).issubset(df.columns):
            print(
                'ERROR: missing on or more needed columns for calculation:\n'
                '{}'.format(', '.join(input_vars))
            )
            return
 
        # run each input using n processors
        inputs = df[input_vars].values.tolist()
        if not nproc:
            nproc = mp.cpu_count() // 2 # use half cores
        pool = mp.Pool(processes=nproc)
        results = pool.map(_calc,inputs)
        pool.close()
        pool.join()

        results_df = pd.concat(results)
        output_vars = ['E', 'Ce', 'VPD', 'stability']
        self._df[output_vars] = results_df[output_vars]
        for el in ['date', 'SH', 'dt']:
            if el in self._df:
                self._df.drop(el, axis=1, inplace=True)

    @property
    def df(self):
        """
        :obj:`pandas.DataFrame` containing input time series meterological data
        and calculated variables after running :meth:`Aero.run`.
        """
        if isinstance(self._df, pd.DataFrame):
            return self._df

    @df.setter
    def df(self, df):
        if not isinstance(df, pd.DataFrame):
            raise TypeError("Must assign a pandas.DataFrame object")
        self._df = df

    @staticmethod
    def single_calc(datetime, wind, pressure, T_air, T_skin, RH, 
            sensor_height, timestep):
        """
        Estimates open water evaporation using the aerodynamic mass transfer
        approach for a single time period.
        
        Arguments:
            datetime (datetime.datetime or str): date-time of measurements for
                error logging.
            wind (float): windspeed in m/s
            pressure (float): air pressure in mbar
            T_air (float): air temperature in C
            T_skin (float): skin temperature (water surface) in C
            RH (float): relative humidity (0-100)
            sensor_height (float): sensor height in m
            timestep (int or float): measurement frequency in seconds

        Returns (tuple): 
            evaporation (mm/timestep), bulk transfer coefficient (Ce), vapor pressure deficit (kPa), and MOST stability (z/L)

        
        """
       
        check=np.array(
            [wind, pressure, T_air, T_skin, RH, sensor_height, timestep]
        )
        if np.isnan(check).any():
            #print('One or more variables missing on {}'.format(datetime))
            return np.nan, np.nan, np.nan

        ###########################################################
        #Constants
        K=0.41 #von Karman constant
        g=9.81 #gravity (m/s^2)
        a=0.0123 #Charnock constant
        
        ############################################################
        #Calculate meterological variables
        #Sensort height (m)
        z=sensor_height
        
        #Convert from Celcius to Kelvin
        T_air=T_air+273.15
        T_skin=T_skin+273.15
        
        #Potential temperatures (air and skin) Kelvin
        T_air_pot=(T_air)*(1000./pressure)**(0.286)
        T_skin_pot=(T_skin)*(1000./pressure)**(0.286)
        
        #Atmospheric vapor pressure (kPa) (2m)
        e_air=(RH/100)*(0.6108*m.exp(((17.27*(T_air-273.15))/((T_air-273.15)+237.3))))
        
        #Atmospheric specific humidity (kg/kg) (2m)
        q_air=0.62*e_air/(pressure/10-0.38*e_air) 
        
        #Saturated Water-Surface vapor pressure (kPa) (0m)
        e_sat=0.6108*m.exp(((17.27*(T_skin-273.15))/((T_skin-273.15)+237.3)))
        
        #Saturated specific humidity at water surface (kg/kg)  (0m)
        q_sat=0.62*e_sat/(pressure/10-0.38*e_sat)
        
        #Vapor Pressure Deficit (kPa)
        VPD=e_sat-e_air
            
        #Density of air (kg/m^3)
        density_air=(pressure/10*1000)/((T_air)*286.9*(1+0.61*q_air))
        
        #Kinematic viscocity
        #Estimated using data from Montgomery, 1947 in Verburg, 2010
        v=(4.94*10**-8*(T_air-273.15)+1.7185*10**-5)/density_air        

        #Virtual Temperature
        Tv=T_air*(1+q_air*0.61)
        ###############################################################
        # Bulk Transfer Coefficient Iteration, Ce
        
        #Stable Condition (z/L > 0)
        #Initial Values for Iteration
        #Stability Function (Momentum)
        Sm=0
        #Stability Function (Temperature)
        St=0
        #Stability Function (Vapor)
        Sq=0
        #Friction velocity
        Us=0
        #Roughness Length of Momentem
        zo=.00010
        
        #Friction Velocity of Momentum
        u_f=(K*(wind-Us))/(cm.log(z/zo)-Sm)
        #Roughness Length of Vapor
        zoq=7.4*zo*cm.exp(-2.25*(zo*u_f/v)**.25)
        #Roughness Legnth of Temperature
        zot=7.4*zo*cm.exp(-2.25*(zo*u_f/v)**.25)
        #Scaling Potential Temperature
        t_fv=(K*(T_air_pot-T_skin_pot))/(cm.log(z/zot)-St)
        #Scaling Humidity
        q_f=np.divide(K*(q_air-q_sat), cm.log(z/zoq)-Sq)
        #Monin-Obhukov Length
        L=np.divide(Tv*u_f**2, K*g*t_fv)
        
        try:
            # avoid crash with bad values causing log of 0 or neg values
            for x in range(0, 199):
                #Friction Velocity of Momentum
                u_f=np.divide((K*(wind-u_f)),(cm.log(z/zo)-Sm))
                #Scaling Potential Temperature
                t_fv=np.divide((K*(T_air_pot-T_skin_pot)),(cm.log(z/zot)-St))
                #Scaling Humidity
                q_f=np.divide((K*(q_air-q_sat)),(cm.log(z/zoq)-Sq))
                #Stability Function of Momentum
                Sm=np.float64(-5.2*(z))/L
                #Stability Function of Vapor
                Sq=np.divide(np.float64(-5.2*(z)),L)
                #Roughness Length of Momemtum
                zc=np.divide(a*u_f**2,g)
                zs=np.divide(0.11*v,u_f)
                zo=zc+zs;
                #Roughness Length of Vapor
                zoq=7.4*zo*cm.exp(-2.25*(zo*u_f/v)**.25)
                #Monin-Obhukov Length
                L=np.divide((Tv*u_f**2),(K*g*t_fv))
        except:
            print('Could not converge on {} for stable conditions'.format(
                    datetime
                )
            )
            return np.nan, np.nan, np.nan

        stability_s = z/L

        if ~np.isreal(L):
            Ce_s=np.nan
        else:
            if np.real(stability_s) > 0:
                Ce_s=np.real(K**2/((cm.log(z/zo)-Sm)*(cm.log(z/zoq)-Sq)))
            else:
                Ce_s=np.nan
        ###############################################################
        #Unstable Conditions (z/L< 0)
        #Initial Values for Iteration
        #Stability Function (Momentum)
        Sm=0
        #Stability Function (Temperature)
        St=0
        #Stability Function (Vapor)
        Sq=0
        #Friction velocity
        Us=0
        #Roughness Length of Momentem
        zo=.00010
        
        #Friction Velocity of Momentum
        u_f=(K*(wind-Us))/(m.log(z/zo)-Sm);
        #Roughness Length of Vapor
        zoq=7.4*zo*m.exp(-2.25*(zo*u_f/v)**.25);
        #Roughness Legnth of Temperature
        zot=7.4*zo*m.exp(-2.25*(zo*u_f/v)**.25);
        #Scaling Potential Temperature
        t_fv=(K*(T_air_pot-T_skin_pot))/(m.log(z/zot)-St);
        #Scaling Humidity
        q_f=(K*(q_air-q_sat))/(m.log(z/zoq)-Sq);
        #Monin-Obhukov Length
        L=np.divide((Tv*u_f**2),(K*g*t_fv))
        
        try:
            for x in range(0, 199):
                #Friction Velocity of Momentum
                u_f=np.divide((K*(wind-u_f)),(cm.log(z/zo)-Sm))
                #Scaling Temperature
                t_fv=np.divide((K*(T_air_pot-T_skin_pot)),(cm.log(z/zot)-St))
                #Scaling Humidity
                q_f=np.divide((K*(q_air-q_sat)),(cm.log(z/zoq)-Sq))
                #Input for Stability function calculations
                x=(1-16*(z/L))**.25
                #Stability Function of Momentum
                Sm=2*cm.log((1+x)/2)+cm.log((1+x**2)/2)-2*cm.atan(x)+(m.pi/2)
                #Stability Function of Vapor
                Sq=2*cm.log((1+x**2)/2)
                #Roughness Length of Momemtum
                zc=np.divide(a*u_f**2,g)
                zs=np.divide(0.11*v,u_f)
                zo=zc+zs
                #Roughness Length of Vapor
                zoq=7.4*zo*cm.exp(-2.25*(zo*u_f/v)**.25)
                #Monin-Obhukov Length
                L=np.divide((Tv*u_f**2),(K*g*t_fv))
        except:
            print('Could not converge on {} for unstable conditions'.format(
                    datetime
                )
            )
        
        stability_u = z/L
        if ~np.isreal(L):
            Ce_u=np.nan
        else:
            if np.real(stability_u) < 0:
                Ce_u=np.real(K**2/((cm.log(z/zo)-Sm)*(cm.log(z/zoq)-Sq)))
            else:
                Ce_u=np.nan
        #################################################################
        #Neutral Conditions, z/L=0
        #Initial Conditions
        zo=.00010
        
        try:
            # avoid crash with bad values causing log of 0 or neg values
            for x in range(0,199):
                #Friction Velocity of Momentum
                u_f=np.divide((K*wind),(np.emath.log(z/zo)))
                #Roughness Length of Momemtum
                zc=np.divide((a*u_f**2),g)
                zs=np.divide(0.11*v,u_f)
                zo=zc+zs
                #Roughness Length of Vapor
                zoq=7.4*zo*m.exp(-2.25*(zo*u_f/v)**.25)
        except:
            print('Could not converge on {} for neutral conditions'.format(
                    datetime
                )
            )

        Ce_n=np.divide(
            (K**2),
            ((np.emath.log(np.divide(z,zo)))*(np.emath.log(np.divide(z,zoq))))
        )
        ################################################################
        #Assign correct Ce value (stable, unstable, or neutral)
        
        if cm.isfinite(Ce_s):
            Ce=Ce_s
            stability = stability_s
        else:
            if cm.isfinite(Ce_u):
                Ce=Ce_u
                stability = stability_u
            else:
                Ce=Ce_n
                stability = 0
        ################################################################
        #Calculated evaporation in mm/timestep
        E=density_air*Ce*(q_sat-q_air)*wind*timestep

        return E, Ce, VPD, np.real(stability)
        
    
def _calc(input_list):
    """
    Helper function for parrallel calculations, needs to be at top-level for
    pickling.

    date = input_list[0]
    WS = input_list[1]
    P = input_list[2]
    T_air = input_list[3]
    T_skin = input_list[4]
    RH = input_list[5]
    SH = input_list[6]
    dt = input_list[7]
    """
    date = input_list[0]
    return pd.DataFrame(
        index=[date], columns=['E', 'Ce', 'VPD', 'stability'], data=[
            Aero.single_calc(
                date, input_list[1], input_list[2], input_list[3], 
                input_list[4], input_list[5],input_list[6],input_list[7]
            )
        ]
    )
