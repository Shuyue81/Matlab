function [h]=heatcoefficient(D,E,F,z)

Vmet=load('wind speed 24.txt');%wind speed (m/s) measured at the meteorological station%

amet=0.14;%wind speed profile exponent at the meteorological station%
sigmamet=270;%wind speed profile boundary layer thickness at the meteorological station (m/s)%
zmet=10;%height above ground of the wind speed sensor at the meteorological station (m)%

a=0.33;%wind speed profile exponent at the site%
sigma=460;%wind speed profile boundary layer thickness at the site (m/s)%

Vz=Vmet.*(sigmamet/zmet)^amet*(z/sigma)^a;%wind speed (m/s) measured at the site%

h=D+Vz.*E+F*Vz.^2;%h is convection heat transfer coefficiency of exterior wall surface%
%D,E,F are material roughness coefficient%

end

