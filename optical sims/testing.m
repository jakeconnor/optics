%% testing: for testing fdfiber class for simulation
function [Voltage] = testing()
clear all;
timesteps=4000;
fiberlength=1;
beta1=6;
N=6;
loss=0;
deltat=0.001;


fiberOne=fdfiber(fiberlength,beta1,N,loss);
fiberOne.plotting=0; %so we get output on the fly, for testing purposes
lasttime=0;
for n=1:timesteps
    Vin=sin(n*deltat);
    fiberOne=fiberOne.simulateFiber(Vin, deltat*n - lasttime);
    lasttime = deltat*n;
    Voltage(n)=fiberOne.Vout;
end
figure;
plot(deltat*(1:length(Voltage)),real(Voltage),'r-',deltat*(1:length(Voltage)),imag(Voltage),'g-');
title('Vout');
xlabel('Time');
ylabel('Amplitude');

legend('Real','Imaginary');
end