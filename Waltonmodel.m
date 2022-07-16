clear

Lw=0.2;
kw=0.84;
pw=1700;
cw=800;
pcw=pw*cw;
aw=1/(pcw/kw);
t0=0;
deltat=600;
deltaxw=0.05;
F0w=aw*deltat/(deltaxw^2);
t=86400;
m=(t-t0)/deltat;
nw=Lw/deltaxw+1;
Tw=zeros(m,nw);
Tw(1,1:nw)=25;

for i=1:m-1
    for j=1:nw
        if j==1
        Tw(i+1,j)=Tw(i,j)*(1-F0w)+F0w*Tw(i,j+1);
        elseif j>1&&j<nw
            Tw(i+1,j)=F0w*(Tw(i,j+1)+Tw(i,j-1))+Tw(i,j)*(1-2*F0w);
        elseif j==nw
            Tw(i+1,j)=Tw(i,j)*(1-3*F0w)+F0w*Tw(i,j-1);
        end
    end
end


Aw=15;
Pw=16;
Dhw=4*Aw/Pw;
Twall=Tw(1:144,1);
Tout=load('dry bulb 144.txt');
dTw1(i)=25-Tw1(i,j);
hinw1=1.31*abs(25-Tw1(i,j))^(1/3);
qw1in=1.31*abs(25-Tw1(i,j))^(1/3)*15*(25-Tw1(i,j));%dTw1(i)=25-Tw1(i,j)    hinw1=1.31*abs(25-Tw1(i,j))^(1/3)%


Lr=0.2;
kr=1.13;
pr=2000;
cr=1000;
pcr=pr*cr;
ar=1/(pcr/kr);
deltaxr=0.05;
F0r=ar*deltat/(deltaxr*deltaxr);
nr=Lr/deltaxr+1;
Tr=zeros(m,nr);
Tr(1,1:nr)=25;
for i=1:m-1
    for j=1:nr
        if j==1
        Tr(i+1,j)=Tr(i,j)*(1-F0r)+F0r*Tr(i,j+1);
        elseif j>1&&j<nr
            Tr(i+1,j)=F0r*(Tr(i,j+1)+Tr(i,j-1))+Tr(i,j)*(1-2*F0r);
        elseif j==nr
            Tr(i+1,j)=Tr(i,j)*(1-3*F0r)+F0r*Tr(i,j-1);
        end
    end
end



Ar=25;
Pr=20;
Dh_r=4*Ar/Pr;
Troof=Tr(1:144,1);
dTr=Troof-Tout;
hinr=1.81*abs(25-Tr(i,j))^(1/3)/(1.328+1);

timeminute=load('time minute.txt');

plot(timeminute, hinw)
hold on
plot(timeminute,hinr)



