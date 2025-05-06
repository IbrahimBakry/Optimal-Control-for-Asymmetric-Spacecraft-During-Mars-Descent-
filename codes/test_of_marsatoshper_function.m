% Test of marsatmospher function
clc;
clear;
clf;

i = 0;
for h = 0:1000:100000
    i = i +1; ii(i) = i;
    [rho,t,p] = marsatmoshper(h);
    Rho(i) = rho;
    T(i) = t;
    P(i) = p;
end

figure(1); plot(ii,Rho,'linewidth',2); title('Density')
figure(2); plot(ii,T,'linewidth',2); title('Temperature')
figure(3); plot(ii,P,'linewidth',2); title('Pressure')