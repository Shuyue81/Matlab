b = 200;
t = 0:1:b;
V = 180; % room volume (m3)
x_initial = 400; % initial pollutant concentration at T =0 (ppm) = 400ppm
x0 = 400; %background pollutant concentration (ppm)
gt = 6; % emission rate (cm3/s)
air_flow_rate_1200 = 0.0075; % air change rate plus deposition rate (m3/s/person)
air_flow_rate_800 = 0.015;
fai_1200 = air_flow_rate_1200*30*3600/V;
xt_1200 = x0 + (gt/air_flow_rate_1200) + (x_initial - x0 - (gt/air_flow_rate_1200))*exp(-fai_1200.*t/60);
fai_800 = air_flow_rate_800*30*3600/V;
xt_800 = x0 + (gt/air_flow_rate_800) + (x_initial - x0 - (gt/air_flow_rate_800))*exp(-fai_800.*t/60);


figure
x = t;
y1 = xt_1200;
y2 = xt_800;
plot(x,[y1;y2]);
axis([0,200,0,1300]);
title("Carbon Dioxide Concentration Change","FontSize",14);
xlabel("Time(min)","FontSize",14);
ylabel("Carbon Dioxide Concentration (ppm)","FontSize",14);
legend("carbon dioxide concentration(air change rate = 0.0075m3/s/person)","carbon dioxide concentration(air change rate = 0.015m3/s/person)","FontSize",14);


xt = 700:1:1200;
x0_co2 = 400;
gt = 6;
air_flow_rate = gt./(xt-x0_co2);

figure
x = xt;
y = air_flow_rate;
plot(x,y);
title("Ventilation Rate Change With Carbon Dioxide Concentration","FontSize",14);
xlabel("Carbon Dioxide Concentration (ppm)","FontSize",14);
ylabel("Ventilation Rate(m^{3}/s/person)","FontSize",14);
legend("Ventilation Rate","FontSize",14);

virus_emission_gt = 10;%virus emission copies/s/person
xt_virus = virus_emission_gt./air_flow_rate; % concentration of virus in air 

figure
x = air_flow_rate;
y = xt_virus;
plot(x,y);
title("Concentration Of Virus Change With Ventilation","FontSize",14);
xlabel("Ventilation Rate(m^{3}/s/person)","FontSize",14);
ylabel("Concentration Of Virus In Air(copies/m^{3})","FontSize",14);
legend("Concentration Of Virus In Air","FontSize",14);

