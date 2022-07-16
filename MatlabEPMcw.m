clear


%%%%%%%%%%%%%%%%%%%%%%%%% solar radiation %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I0=1355;%I0 is global solar constant%
c0=0.5598;
c1=0.4982;
c2=-0.6762;
c3=0.02842;
c4=-0.00317;
c5=0.014;
d=-17.853;
k=0.843;%c0,c1,c2,c3,c4,c5,d,k are regression coefficients%
CC=load('cloud cover 24.txt');%CC is the cloud cover in tenth on July 21st%
h=load('altitude 24.txt');%h is solar altitude on July 21st(degree)%
fai=load('relative humidity 24.txt');% fai is relative humidity on July 21st(%)%
Vw=load('wind speed 24.txt');%Vw is wind speed on July 21st(m/s)% 

T=load('dry bulb 24.txt');%T is dry bulb temperature of air on July 21st(degree C)%
Tbefore=load('dry bulb before 24.txt');%Tbefore is dry bulb temperature of air from July 20th 21pm. to July 21st 21pm.(degree C)%
azimuth=load('solar azimuth 24.txt');%azimuth is solar azimuth on July 21st(degree)%


for n=1:24%24 hours during on day%
   
I(n)=(I0*sind(h(n))*(c0+c1*CC(n)/8+c2*(CC(n)/8)^2+c3*(T(n)-Tbefore(n))+c4*fai(n)+c5*Vw(n))+d)/k;
%equation to calculate estimated hourly solar radiation(W/m^2)%
if I(n)<0
    I(n)=0;
end
%solar radiation should be bigger than 0, so if calculated solar radiation
%is smaller than 0, assume solar radiation as 0(W/m^2)%

if h(n)==0
    KT(n)=0;
else
    KT(n)=I(n)/(I0*sind(h(n)));
end

KTC(n)=0.4268+0.1934*sind(h(n));

if KT(n)>=KTC(n)
    KDS(n)=KT(n)-(1.107+0.03569*sind(h(n))+1.681*(sind(h(n)))^2)*(1-KT(n))^2;
else 
    KDS(n)=(3.996-3.862*sind(h(n))+1.540*(sind(h(n)))^2)*KT(n)^3;
end
%This step is to determine the formular of KDS%


if h(n)==0
    DHw1(n)=0;
    DHw2(n)=0;
     DHr(n)=I0*sind(h(n))*KDS(n)*(1-KT(n))/(1-KDS(n));
     SH(n)=I0*sind(h(n))*(KT(n)-KDS(n))/(1-KDS(n));
else
 DHr(n)=I0*sind(h(n))*KDS(n)*(1-KT(n))/(1-KDS(n));%direct solar radiation on the roof (W/m^2)%
DHw1(n)=DHr(n)*abs(cosd(0-azimuth(n))/tand(h(n)));%direct solar radiation (W/m^2) on the northern vertical wall (wall 1)%
DHw2(n)=DHr(n)*abs(cosd(90-azimuth(n))/tand(h(n)));%direct solar radiation (W/m^2) on the eastern vertical wall (wall 2)%
SH(n)=I0*sind(h(n))*(KT(n)-KDS(n))/(1-KDS(n));%diffuse radiation (W/m^2)%

end
end





%%%%%%%%% convection transfer heat coefficient of exterior wall surface %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hw_ex=heatcoefficient(12.49,4.065,0.028,1.5);
%material of wall is brick, so it is rough, thus, D=12.49,E=4.065,F=0.028%
hr_ex=heatcoefficient(10.79,4.192,0,3);
%material of roof is concrete, so it is very rough, thus, D=11.58,E=5.894,F=0%

time=load('time.txt');%24 hours%

 plot(time,hw_ex)
 hold on
 plot(time,hr_ex)
 xlabel('Time (h)')
 ylabel('heat transfer coefficient h')
 plot the figure of convection heat coefficient of wall and roof agianst time%




%%%%%%this step is to expand 1*24 matrics into 144*1 matrics %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x=1:24
    for y=6*x-5:6*x
        DHr144(y,1)=DHr(1,x);
        DHw1144(y,1)=DHw1(1,x);
        DHw2144(y,1)=DHw2(1,x);
        SH144(y,1)=SH(1,x);
        hw_ex144(y,1)=hw_ex(x,1);
        hr_ex144(y,1)=hr_ex(x,1);
        T144(y,1)=T(x,1);
    end
end
%change the time interval from one hour to ten minutes%




%%%%%%%%%%%%%%long wave radiation for interior wall%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(1)=15;A(2)=15;A(3)=15;A(4)=15;A(5)=25;A(6)=25;F=zeros(10,6);F(1,:)=1;

for c=1:10
    for e=1:6
        F(c+1,e)=1/(1-A(e)*F(c,e)/(15*(F(c,1)+F(c,2)+F(c,3)+F(c,4))+25*(F(c,5)+F(c,6))));
    end
end




[z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15]=Matlabcw20(25,25,25);



%%%%%%%%%%%%%%%%%%%%%% northern wall (wall 1) %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lw=0.2;%the thickness of vertical wall (m)%
kw=0.84;%thermal conductiviy of wall(W/m)%
pw=1700;%dencity of wall (kg/m^3)%
cw=800;%specific heat capacity of wall (J/(kg*K))%
pcw=pw*cw;%the product of dencity and specific heat capacity of wall (J/(K*m^3))%
aw=1/(pcw/kw);%to simplify the calculation%
t0=0;%when thime is 0 am. on July 21st%
deltat=600;%the time interval is 600 seconds%
deltaxw=0.04;%the distance interval of wall is 0.04m%
F0w=aw*deltat/(deltaxw^2);%to calculate fourior number%
t=86400;%number of seconds in the whole day%
m=(t-t0)/deltat;%number of time nodes%
nw=Lw/deltaxw;%number of distance nodes of wall%
Tw1=zeros(m,nw);%matrics of temperature of vertical wall 1 at different time nodes and distance nodes%
Tw1(1,1)=z1;Tw1(1,2)=z2;Tw1(1,3)=z3;Tw1(1,4)=z4;Tw1(1,5)=z5;%initial tmperature of vertical wall 1 is 25 degree C%



%%%%%%%%%%%%%%%%%%%%%eastern wall (wall 2) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tw2=zeros(m,nw);%matrics of temperature of vertical wall 2 at different time nodes and distance nodes%
Tw2(1,1)=z6;Tw2(1,2)=z7;Tw2(1,3)=z8;Tw2(1,4)=z9;Tw2(1,5)=z10;%initial tmperature of vertical wall 1 is 25 degree C%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% roof %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lr=0.2;%the thickness of roof (m)%
kr=1.13;%thermal conductiviy of roof(W/m)%
pr=2000;%dencity of roof (kg/m^3)%
cr=1000;%specific heat capacity of roof (J/(kg*K))%
pcr=pr*cr;%the product of dencity and specific heat capacity of roof (J/(K*m^3))%
ar=1/(pcr/kr);%to simplify the calculation%
deltaxr=0.04;%the distance interval of roof is 0.04m%
F0r=ar*deltat/(deltaxr*deltaxr);%to calculate fourior number%
nr=Lr/deltaxr;%number of distance nodes%
Tr=zeros(m,nr);%matrics of temperature of roof at different time nodes and distance nodes%
Tr(1,1)=z11;Tr(1,2)=z12;Tr(1,3)=z13;Tr(1,4)=z14;Tr(1,5)=z15;%initial tmperature of roof is 25 degree C%
Tw3=zeros(144,1);
Tw3(:,1)=25;

for i=1:m-1
    for j=1:nw
       
        if j==1
            
            %calcutaion of interior long wave radiation for each side of wall%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            epsilon=0.9;
            for e=1:6
               href(e)=4*5.67*10^(-8)*300^3/(1/F(11,e)+(1-epsilon)/epsilon);
            end
            
            
hradw1_1(i)=(0.865+Tw1(i,j)/200)*href(1);%first estimates of radiation of each wall%
hradw2_1(i)=(0.865+Tw2(i,j)/200)*href(2);
hradw3_1(i)=(0.865+Tw3(i,j)/200)*href(3);
hradw4_1(i)=(0.865+Tw3(i,j)/200)*href(4);
hradr_1(i)=(0.865+Tr(i,j)/200)*href(5);
hradf_1(i)=(0.865+Tw3(i,j)/200)*href(6);

TMR_1(i)=(15*hradw1_1(i)*Tw1(i,j)+15*hradw2_1(i)*Tw2(i,j)+15*hradw3_1(i)*25+15*hradw4_1(i)*25+25*hradr_1(i)*Tr(i,j)+25*hradf_1(i)*25)/(15*hradw1_1(i)+15*hradw2_1(i)+15*hradw3_1(i)+15*hradw4_1(i)+25*hradr_1(i)+25*hradf_1(i));
%the first estimate of MRT%

hradw1_2(i)=(0.865+TMR_1(i)/200)*hradw1_1(i);%the second estimate of each wall%
hradw2_2(i)=(0.865+TMR_1(i)/200)*hradw2_1(i);
hradw3_2(i)=(0.865+TMR_1(i)/200)*hradw3_1(i);
hradw4_2(i)=(0.865+TMR_1(i)/200)*hradw4_1(i);
hradr_2(i)=(0.865+TMR_1(i)/200)*hradr_1(i);
hradf_2(i)=(0.865+TMR_1(i)/200)*hradf_1(i);

TMR_2(i)=(15*hradw1_2(i)*Tw1(i,j)+15*hradw2_2(i)*Tw2(i,j)+15*hradw3_2(i)*25+15*hradw4_2(i)*25+25*hradr_2(i)*Tr(i,j)+25*hradf_2(i)*25)/(15*hradw1_2(i)+15*hradw2_2(i)+15*hradw3_2(i)+15*hradw4_2(i)+25*hradr_2(i)+25*hradf_2(i));
%the second estimate of MRT%

qradflux_w1(i)=hradw1_2(i)*(Tw1(i,j)-TMR_2(i));%the long-wave radiation heat transfer flux of internal surface, W/m^2%
qradflux_w2(i)=hradw2_2(i)*(Tw2(i,j)-TMR_2(i));
qradflux_r(i)=hradr_2(i)*(Tr(i,j)-TMR_2(i));


hinconw1(i)=1.31*abs(25-Tw1(i,j))^(1/3);%the heat tranfer coefficient of convection of internal surface%
hinconw2(i)=1.31*abs(25-Tw2(i,j))^(1/3);
hinconr(i)=1.81*abs(25-Tr(i,j))^(1/3)/(1.328+1);


qinconw1(i)=hinconw1(i)*(Tw1(i,j)-25);%the heat tranfer energy of convection of internal surface, W/m^2%
qinconw2(i)=hinconw2(i)*(Tw2(i,j)-25);
qinconr(i)=hinconr(i)*(Tr(i,j)-25);

qconvection(i)=qinconr(i)*25+qinconw1(i)*15+qinconw2(i)*15;%total energy of convection of internal surface of three surfaces%


       Tw1(i+1,j)=Tw1(i,j)*(1-F0w)+F0w*Tw1(i,j+1)-F0w*(qinconw1(i)+qradflux_w1(i))/kw*deltaxw;
       Tw2(i+1,j)=Tw2(i,j)*(1-F0w)+F0w*Tw2(i,j+1)-F0w*(qinconw2(i)+qradflux_w2(i))/kw*deltaxw;
       Tr(i+1,j)=Tr(i,j)*(1-F0r)+F0r*Tr(i,j+1)-F0r*(qinconr(i)+qradflux_r(i))/kr*deltaxr;
       %to calculate the interior surface's temperature at each time node for wall 1, wall 2 and roof%
       
        elseif j>1&&j<nw
            Tw1(i+1,j)=F0w*(Tw1(i,j+1)+Tw1(i,j-1))+Tw1(i,j)*(1-2*F0w);
            Tw2(i+1,j)=F0w*(Tw2(i,j+1)+Tw2(i,j-1))+Tw2(i,j)*(1-2*F0w);
             Tr(i+1,j)=F0r*(Tr(i,j+1)+Tr(i,j-1))+Tr(i,j)*(1-2*F0r);
            %to calculate the temperature of distance node 2, 3,4 at each time node for wall1, wall 2 and roof%
            
        elseif j==nw
           
         
            qradw1_ex(i)=DHw1144(i+1)+SH144(i+1);
            qradw2_ex(i)=DHw2144(i+1)+SH144(i+1);
            qradr_ex(i)=DHr144(i+1)+SH144(i+1);
            
            
            qconw1_ex(i)=hw_ex144(i)*(T144(i)-Tw1(i,j));
            qconw2_ex(i)=hw_ex144(i)*(T144(i)-Tw1(i,j));
            qconr_ex(i)=hr_ex144(i)*(T144(i)-Tw1(i,j));
            
            Tw1(i+1,j)=Tw1(i,j)*(1-F0w)+F0w*Tw1(i,j-1)+F0w*deltaxw*(qradw1_ex(i)*0.7+qconw1_ex(i))/kw;
             Tw2(i+1,j)=Tw2(i,j)*(1-F0w)+F0w*Tw2(i,j-1)+F0w*deltaxw*(qradw2_ex(i)*0.7+qconw2_ex(i))/kw;
             Tr(i+1,j)=Tr(i,j)*(1-F0r)+F0r*Tr(i,j-1)+F0r*deltaxr*(qradr_ex(i)*0.7+qconr_ex(i))/kr;
            %to calculate the xterior surface's temperature at each time node for wall 1, wall 2 and roof%
        end
        
    end
end



%plot(1:1:143,qconvection)


