%%C. Pearson, Desert Research Institute, 10/2015
%This function estimates open water evaporation by the bulk aerodynamic 
%tranfer approach. The bulk transfer coefficient, Ce, is solved for iteratively
%based on Monin-Obukhov Similarity equations presented in Brutsaert, 1982 'Evaporation into the Atmosphere'
%Input Variable Units:
%Wind: m/s
%Pressure: mb
%T_air: C
%T_skin: C
%RH: Percent
%Sensor height: meters
%Timestep: sampling rate in seconds
%Returns E in mm/timestep; Ce; VPD (kPa) 


%testdata=[731947, 11.7661091674962, 920.0627, 0.4288, 6.3593, 82.2214, 2, 86400]

function [E, Ce, VPD]=aero(datetime, wind, pressure, T_air, T_skin, RH, sensor_height, timestep) 
% datetime=7.3573e+05
% wind=0.597000000000000
% pressure=878.400000000000
% T_air=5.72000000000000
% T_skin=5.42425040400000
% RH=91.1000000000000
% sensor_height=2
% timestep=1800

%Returns NaNs if input data contains NaNs
check=[datetime, wind, pressure, T_air, T_skin, RH, sensor_height, timestep];
if sum(isnan(check))>0
    Ce=NaN;
    E=NaN;
    VPD=NaN;
    Ce_all=[NaN,NaN,NaN];
return
end

%Constants
K=0.41; %von Karman constant
g=9.81; %gravity (m/s^2)
a=0.0123; %Charnock constant

% Sensort height
z=sensor_height;

%Convert from C to Kelvin
T_air=T_air+273.15;
T_skin=T_skin+273.15;

%Potential temperatures (air and skin) Kelvin
T_air_pot=(T_air).*(1000./pressure).^(0.286);
T_skin_pot=(T_skin).*(1000./pressure).^(0.286);

%Atmospheric vapor pressure (kPa) (2m)
e_air=(RH/100).*(0.6108*exp(((17.27*(T_air-273.15))./((T_air-273.15)+237.3))));

%Atmospheric specific humidity (kg/kg) (2m)
q_air=0.62*e_air./(pressure/10-0.38*e_air); 

%Saturated Water-Surface vapor pressure (kPa) (0m)
e_sat=0.6108*exp(((17.27*(T_skin-273.15))./((T_skin-273.15)+237.3)));

%Saturated specific humidity at water surface (kg/kg)  (0m)
q_sat=0.62*e_sat./(pressure/10-0.38*e_sat);

%Vapor Pressure Deficit (e_sat - e_air); (kPa)
VPD=e_sat-e_air;

%Density of air (kg/m^3)
density_air=(pressure/10*1000)./((T_air).*286.9.*(1+0.61*q_air));

%Kinematic viscocity
v=(4.94*10^-8*(T_air-273.15)+1.7185*10^-5)./density_air; % Estimated using data from Montgomery, 1947 in Verburg, 2010

%Virtual Temperature
Tv=T_air*(1+q_air*0.61);

%% Bulk Transfer Coefficient Iteration, Ce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stable Conditions (z/L > 0)

%Initial Values for Iteration
%Stability Function (Momentum)
Sm=0;
%Stability Function (Temperature)
St=0;
%Stability Function (Vapor)
Sq=0;
%Friction velocity
Us=0;
%Roughness Length of Momentem
zo=.00010;

%Friction Velocity of Momentum
u_f=(K*(wind-Us))/(log(z/zo)-Sm);
%Roughness Length of Vapor
zoq=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
%Roughness Legnth of Temperature
zot=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
%Scaling Potential Temperature
t_fv=(K*(T_air_pot-T_skin_pot))/(log(z/zot)-St);
%Scaling Humidity
q_f=(K*(q_air-q_sat))/(log(z/zoq)-Sq);
%Monin-Obhukov Length
L=(Tv*u_f^2)/(K*g*t_fv);

for i=1:200
    %Friction Velocity of Momentum
    u_f=(K*(wind-u_f))/(log(z/zo)-Sm); %(7)
    %Scaling Potential Temperature
    t_fv=(K*(T_air_pot-T_skin_pot))/(log(z/zot)-St); %(8)
    %Scaling Humidity
    q_f=(K*(q_air-q_sat))/(log(z/zoq)-Sq); %(9)
    %Stability Function of Momentum
    Sm=-5.2*(z/L);
    %Stability Function of Vapor
    Sq=-5.2*(z/L);
    %Roughness Length of Momemtum
    zc=a*u_f^2/g;
    zs=0.11*v/u_f;
    zo=zc+zs;
    %Roughness Length of Vapor
    zoq=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
    %Monin-Obhukov Length
    L=(Tv*u_f^2)/(K*g*t_fv);
end
%Monin-Obukov Stability Parameter 
stability=z/L;
%Check for Stable Conditions (z/L>0)
if stability>0
%     Ce_all(1,1)=q_f*u_f/((wind-u_f)*(q_air-q_sat)); % Eq. Format Test (SAME)
    Ce_all(1,1)=K^2/((log(z/zo)-Sm)*(log(z/zoq)-Sq));
else
    Ce_all(1,1)=NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unstable Conditions (z/L < 0)

%Initial Values for Iteration
%Stability Function (Momentum)
Sm=0;
%Stability Function (Temperature)
St=0;
%Stability Function (Vapor)
Sq=0;
%Friction velocity
Us=0;
%Roughness Length of Momentem
zo=.00010;

%Friction Velocity of Momentum
u_f=(K*(wind-Us))/(log(z/zo)-Sm);
%Roughness Length of Vapor
zoq=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
%Roughness Legnth of Temperature
zot=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
%Scaling Potential Temperature
t_fv=(K*(T_air_pot-T_skin_pot))/(log(z/zot)-St);
%Scaling Humidity
q_f=(K*(q_air-q_sat))/(log(z/zoq)-Sq);
%Monin-Obhukov Length
L=(Tv*u_f^2)/(K*g*t_fv);

for i=1:200
    %Friction Velocity of Momentum
    u_f=(K*(wind-u_f))/(log(z/zo)-Sm); %(7)
    %Scaling Temperature
    t_fv=(K*(T_air_pot-T_skin_pot))/(log(z/zot)-St); %(8)
    %Scaling Humidity
    q_f=(K*(q_air-q_sat))/(log(z/zoq)-Sq); %(9)
    %Input for Stability function calculations
    x=(1-16*(z/L))^.25;
    %Stability Function of Momentum
    Sm=2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+(pi/2);
    %Stability Function of Vapor
    Sq=2*log((1+x^2)/2);
    %Roughness Length of Momemtum
    zc=a*u_f^2/g;
    zs=0.11*v/u_f;
    zo=zc+zs;
    %Roughness Length of Vapor
    zoq=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
    %Monin-Obhukov Length
    L=(Tv*u_f^2)/(K*g*t_fv);
end
%Monin-Obukov Stability Parameter 
stability=z/L;
%Check for Stable Conditions (z/L<0)
if stability<0
%     Ce_all(1,2)=q_f*u_f/((wind-u_f)*(q_air-q_sat)); % Eq. Format Test (SAME)
    Ce_all(1,2)=K^2/((log(z/zo)-Sm)*(log(z/zoq)-Sq));
else
    Ce_all(1,2)=NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Neutral Conditions (z/L=0)

%Initial Conditions
zo=.00010;

for i=1:200
    %Friction Velocity of Momentum
    u_f=K*wind/log(z/zo);
    %Roughness Length of Momemtum
    zc=a*u_f^2/g;
    zs=0.11*v/u_f;
    zo=zc+zs;
%Roughness Length of Vapor
zoq=7.4*zo*exp(-2.25*(zo*u_f/v)^.25);
end

% Ce_all(1,3)=K^2/(log(z/zo)*log(z/zoq));
Ce_all(1,3)=K^2/((log(z/zo))*(log(z/zoq)));

%Selects and creates final array of Ce values based on stability (stable or unstable)
a=Ce_all(1,~isnan(Ce_all(1,:)));
b=isempty(a);
%No Ce value 
if b==1
    Ce=NaN;
else
    %Chooses unstable solution if non-unique solutions are found
    if length(a)>1
        Ce=a(1);
    else
        Ce=a;
    end
end
%Removes imaginary part of solutions (0i)
Ce=real(Ce);
%Remove Negative Values
Ce(find(Ce<0))=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates evaporation in mm/timestep
E=density_air.*Ce.*(q_sat-q_air).*wind*timestep;

disp(E);
disp(Ce);
disp(VPD);
end



