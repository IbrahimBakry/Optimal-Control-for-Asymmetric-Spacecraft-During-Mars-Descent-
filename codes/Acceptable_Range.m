% This code: For Aero and mass asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
clc;
clear all;
close all;
 
mzn0 = -0.01; my0f = 0.3; mz0f = 0.3; 
dzdash = 0.05; dydash = 0.05;
Ixzd = 0.1; Ixyd = 0.1;


%% Acceptaple range
leng = 50;
mz1 = -sin(linspace(0.1,0.3,leng));
Cx1 = cos(linspace(0.1,1.55,leng));
Cy1 = sin(linspace(0.1,0.6,leng));

Ixd = 270 / 443;
p1 = 1.2; p2 = 1.2; p3 = 1.1; p4 = 1.1; p5 = 0.7; p6 = 0.7;
u=3;

% % My
% C = (Cx1/p2 + 1/p4 + 1/p5).^2 + (Cx1/p1 + 1/p3 + 1/p6).^2;
% omega = u*sqrt(C);
% 
% myf01 = -sqrt(omega.^2./C).*mz1./p5;
% mzf01 = -sqrt(omega.^2./C).*mz1./p6;
% Ixyd1 = sqrt(omega.^2./C).*(1-Ixd)./p3;
% Ixzd1 = sqrt(omega.^2./C).*(1-Ixd)./p4;
% dyd1 = -sqrt(omega.^2./C).*mz1./(p1);
% dzd1 = -sqrt(omega.^2./C).*mz1./(p2);
% myf01v = (-mz1./(1-Ixd)).*(p3/p1).*Ixyd; % Мои
% mzf01v = (-mz1./(1-Ixd)).*(p4/p2).*Ixzd; % Мои
% figure; plot(Ixyd1,myf01,'LineWidth',3)
% 
% 1

% Kurkina
p1k = 1.5; p2k = 1.5; p3k = 1; p5k = 0.5; p4k = 1; p6k = 0.5;
A = (Cx1+Cy1)./(p4k*Cy1)+(Cx1+Cy1)./(p3k*Cy1)+1/(p1k*p5k)+1/(p2k*p6k)+1/(p4k*p6k);
Ck = 1/p1k^2 + 1/p2k^2 - 2/(p1k*p3k) + 2/(p2k*p4k) + 1/(4*p3k^2)+ 1/(4*p4k^2);
uk = 25;%5.7
B = abs(A)/(Ck)^0.5;
omegak = B/sqrt(uk^3);

myf02 = -sqrt(omegak.^2./B).^(2/3).*mz1./p5k;
mzf02 = -sqrt(omegak.^2./B).^(2/3).*mz1./p6k;
Ixyd2 = sqrt(omegak.^2./B).^(2/3).*(1-Ixd)./p3k;
Ixzd2 = sqrt(omegak.^2./B).^(2/3).*(1-Ixd)./p4k;
dyd2 = -sqrt(omegak.^2./B).^(2/3).*mz1./(Cy1.*p1k);
dzd2 = -sqrt(omegak.^2./B).^(2/3).*mz1./(Cy1.*p2k);
% myf02 = (-mz1./(1-Ixd)).*(p3k/p5k).*Ixy; % Куркина
% mzf02 = (-mz1./(1-Ixd)).*(p4k/p6k).*Ixz; % Куркина

figure(1);
    plot(Ixyd2,myf02/2.3,'LineWidth',3); hold on
    plot(Ixyd2,myf02/3.15,'LineWidth',3);
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
    hold all; grid on; box on; xlabel('Ixyd2');ylabel('myf0') 
    
%     myf01 = [0.0086 0.0087 0.0089 0.009 0.0091 0.0093 0.0094 0.0096 0.0097 myf02/2.3];
%     Ixyd1 = [Ixyd2 0.0278 0.0277 0.0275 0.0274 0.0273 0.0271 0.0269 0.0278 0.0287];
%     plot(Ixyd1,myf01,'LineWidth',3);
    
%     tt = linspace(0.02,0.078,100); aa1 = polyfit((Ixyd2),(myf02/2.3),1); aa1v = polyval(aa1,tt);
%     plot(tt,aa1v,'LineWidth',3);
    axis([0.025 0.045 0.005 0.02])

    
% figure(2);
%     plot(Ixzd2,dzd2,'LineWidth',3); hold on
%     ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%     ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     hold all; grid on; box on; xlabel('Ixyd2');ylabel('myf0') 

    
