clc;
close all;
clear;

dzdfk = linspace(0.001,0.0037,10);
dydfk = linspace(0.001,0.0037,10);
my0fk = (linspace(0,0.0133,10));
mz0fk = (linspace(0,0.0133,10));
Ixydk = 0.0439;
Ixzdk = 0.0439; 

dzdf = linspace(0.001,0.0183,10);
dydf = linspace(0.001,0.0183,10);
my0f = (linspace(0,0.078,10));
mz0f = (linspace(0,0.078,10));


plot(dydfk,my0fk,dydf,my0f)


