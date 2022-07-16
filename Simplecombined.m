hw=heatcoefficient(12.49,4.065,0.028);
hr=heatcoefficient(10.79,4.792,0);
time=load('time.txt');

plot(time,hw)
hold on
plot(time,hr)
xlabel('Time (h)')
ylabel('heat transfer coefficient h')