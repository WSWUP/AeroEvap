# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:13:07 2015

@author: cpearson
Aerodynmaic mass transfer code to calculate E, Ce, and VPD from buoy input parameters
"""

def aero(datetime,wind,pressure,T_air,T_skin,RH, sensor_height, timestep):
    "Esimates open water evaporation using the aerodynamica mass transfer approach"
#Input Variable Units:
#Wind: m/s
#Pressure: mb
#T_air: C
#T_skin: C
#RH: Percent
#Sensor height: meters
#Timestep: sampling rate in seconds "
#TEST DATA

#Neutral Conditions 
#datetime=7.3573e+05
#wind=0.597000000000000
#pressure=878.400000000000
#T_air=5.72000000000000
#T_skin=5.42425040400000
#RH=91.1000000000000
#sensor_height=2
#timestep=1800

#Unstable (731947,11.7661091674962,920.0627,0.4288,6.3593,82.2214,2,86400)
#datetime=731947
#wind=11.7661091674962
#pressure=920.0627
#T_air=0.4288
#T_skin=6.3593
#RH=82.2214
#RH=float('NaN')
#sensor_height=2
#timestep=86400
   
    import cmath
    import numpy as np
    import math as m
    import sys as sys
    check=np.array([datetime, wind, pressure, T_air, T_skin, RH, sensor_height, timestep])
    if np.isnan(check).any():
       Ce=float('NaN')
       E=float('NaN')
       #print('Incomplete Dataset')
       
       sys.exit()

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
    v=(4.94*10**-8*(T_air-273.15)+1.7185*10**-5)/density_air #Estimated using data from Montgomery, 1947 in Verburg, 2010
    
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
    u_f=(K*(wind-Us))/(cmath.log(z/zo)-Sm)
    #Roughness Length of Vapor
    zoq=7.4*zo*cmath.exp(-2.25*(zo*u_f/v)**.25)
    #Roughness Legnth of Temperature
    zot=7.4*zo*cmath.exp(-2.25*(zo*u_f/v)**.25)
    #Scaling Potential Temperature
    t_fv=(K*(T_air_pot-T_skin_pot))/(cmath.log(z/zot)-St)
    #Scaling Humidity
    q_f=(K*(q_air-q_sat))/(cmath.log(z/zoq)-Sq)
    #Monin-Obhukov Length
    L=(Tv*u_f**2)/(K*g*t_fv)
    
    for x in range(0, 199):
        #Friction Velocity of Momentum
        u_f=(K*(wind-u_f))/(cmath.log(z/zo)-Sm)
        #Scaling Potential Temperature
        t_fv=(K*(T_air_pot-T_skin_pot))/(cmath.log(z/zot)-St)
        #Scaling Humidity
        q_f=(K*(q_air-q_sat))/(cmath.log(z/zoq)-Sq)
        #Stability Function of Momentum
        Sm=np.float64(-5.2*(z))/L
        #Stability Function of Vapor
        Sq=np.float64(-5.2*(z))/L
        #Roughness Length of Momemtum
        zc=a*u_f**2/g;
        zs=0.11*v/u_f;
        zo=zc+zs;
        #Roughness Length of Vapor
        zoq=7.4*zo*cmath.exp(-2.25*(zo*u_f/v)**.25)
        #Monin-Obhukov Length
        L=(Tv*u_f**2)/(K*g*t_fv);
    stability=z/L
    #stability=z/L
    if ~np.isreal(L):
        Ce_s=float('NaN')    
    else:
        if np.real(stability)>0:
            Ce_s=np.real(K**2/((cmath.log(z/zo)-Sm)*(cmath.log(z/zoq)-Sq)))
        else:
            Ce_s=float('NaN')  
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
    L=(Tv*u_f**2)/(K*g*t_fv);
    
    for x in range(0, 199):
        #Friction Velocity of Momentum
        u_f=(K*(wind-u_f))/(cmath.log(z/zo)-Sm)
        #Scaling Temperature
        t_fv=(K*(T_air_pot-T_skin_pot))/(cmath.log(z/zot)-St)
        #Scaling Humidity
        q_f=(K*(q_air-q_sat))/(cmath.log(z/zoq)-Sq)
        #Input for Stability function calculations
        x=(1-16*(z/L))**.25
        #Stability Function of Momentum
        Sm=2*cmath.log((1+x)/2)+cmath.log((1+x**2)/2)-2*cmath.atan(x)+(m.pi/2)
        #Stability Function of Vapor
        Sq=2*cmath.log((1+x**2)/2)
        #Roughness Length of Momemtum
        zc=a*u_f**2/g
        zs=0.11*v/u_f
        zo=zc+zs
        #Roughness Length of Vapor
        zoq=7.4*zo*cmath.exp(-2.25*(zo*u_f/v)**.25)
        #Monin-Obhukov Length
        L=(Tv*u_f**2)/(K*g*t_fv)
    
    stability=z/L
    if ~np.isreal(L):
        Ce_u=float('NaN')
    else:
        if np.real(stability)<0:
            Ce_u=np.real(K**2/((cmath.log(z/zo)-Sm)*(cmath.log(z/zoq)-Sq)))
        else:
            Ce_u=float('NaN')  
    #################################################################
    #Nuetral Conditions, z/L=0
    #Initial Conditions
    zo=.00010
    
    for x in range(0,199):
        #Friction Velocity of Momentum
        u_f=K*wind/m.log(z/zo)
        #Roughness Length of Momemtum
        zc=a*u_f**2/g
        zs=0.11*v/u_f
        zo=zc+zs
        #Roughness Length of Vapor
        zoq=7.4*zo*m.exp(-2.25*(zo*u_f/v)**.25)
    Ce_n=K**2/((m.log(z/zo))*(m.log(z/zoq)));
    ################################################################
    #Assign correct Ce value (stable, unstable, or neutral)
    
    if cmath.isfinite(Ce_s):
        Ce=Ce_s
    else:
        if cmath.isfinite(Ce_u):
            Ce=Ce_u
        else:
            Ce=Ce_n
    ################################################################
    #Calculated evaporation in mm/timestep
    E=np.array(density_air*Ce*(q_sat-q_air)*wind*timestep)
    return E, Ce, VPD
        
#if __name__ == "__main__":
#        sys.exit(main())
    