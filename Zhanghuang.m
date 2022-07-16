clear

I0=1355;
c0=0.5598;
c1=0.4982;
c2=-0.6762;
c3=0.02842;
c4=-0.00317;
c5=0.014;
d=-17.853;
k=0.843;
CC=load('cloud cover 27.txt');
h=load('solar altitude 27.txt');
fai=load('relative humidity 27.txt');
Vw=load('wind speed 27.txt');

T=load('dry bulb 27.txt');
azimuth=load('solar azimuth 27.txt');


for n=4:27
   
I(n)=(I0*sind(h(n))*(c0+c1*CC(n)+c2*CC(n)^2+c3*(T(n)-T(n-3))+c4*fai(n)+c5*Vw(n))+d)/k;
KT(n)=I(n)/(I0*sind(h(n)));
KTC(n)=0.4268+0.1934*sind(h(n));
KDS(n)=KT(n)-(1.107+0.03569*sind(h(n))+1.681*(sind(h(n)))^2)*(1-KT(n))^2;
DHr(n)=I0*sind(h(n))*KDS(n)*(1-KT(n))/(1-KDS(n));
DHw1(n)=DHr(n)*abs(cosd(0-azimuth(n))/tand(h(n)));
DHw2(n)=DHr(n)*abs(cosd(90-azimuth(n))/tand(h(n)));
SH(n)=I0*sind(h(n))*(KT(n)-KDS(n))/(1-KDS(n));

end


