%% The Original Nonlinear Equations
clc
clear all;
close all;

% Controller Coffitiants
Aw = 1; Bw = 1100; % Aw = 1.0; Bw = 10;
Aa = 1.0; Ba = 0.05; % Aa = 1.0; Ba = 10;
Aav= 1.0; Bav= 0.1; % Aav= 1.0; Bav= 1;

% Initial Conditions
V0 = 5500; theta0 = -(19*(pi/180)); h0 = 1e5; 
wx0 = 0.14; wxv0 = 0; z10 = 0; alpha0 = 0.21; phi0 = 0.035;
y0 = [V0;theta0;h0;wx0;wxv0;z10;alpha0;phi0];

tspan = [0 300];

[t,y] = ode113(@(t,y) EQN(t,y,Aw,Bw,Aa,Ba,Aav,Bav),tspan,y0);

V=y(:,1); theta=y(:,2);  h=y(:,3); 
wx=y(:,4); wxv=y(:,5); z1=y(:,6); alpha=y(:,7); phi=y(:,8); 


g0 = 3.72076; Rmars = 3396000; Cxv = 1.75; mxw = 0.01;
r = 1.25; s = pi*r^2; L = 2; m = 576; Ix = 300; I = 443; % Iy=Iy=I;
rho = (.699.*exp(-0.00009.*h))./(.1921.*((-23.4-0.00222.*h)+273.1));
q = 0.5.*(rho).*V.^2;
g = g0.*Rmars.^2./(Rmars+h).^2;
dVdt = (Cxv.*q.*s./m + g.*sin(theta));
Load = dVdt ./ g;

% Control functions;
fw = mxw.*q.*s.^2.*L./(Ix.*V./I);
Kw = fw + sqrt(fw.^2 - Aw./Bw);
% Kw = sqrt(Aw/Bw);
Ka = sqrt(Aa/Ba);
Kav = sqrt(Aav/Bav);
Uw = - Kw.*wx;
Ualpha = - Ka*alpha - Kav*z1;


figure(1); % alpha
    plot(t,abs(alpha-0.01),'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
    yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
    yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
    ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Время, c');ylabel('Угол атаки \alpha_п, Рад') 

figure(2); % Wx
        plot(t,wx,'LineWidth',3);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14' '0.16'})
        yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
        ax.YLim = [0.00 0.16];
%         ax.XLim = [0.00 100];
        hold all; grid on; box on;xlabel('Время, c'); ylabel('Угловая скорость \omega_x, 1/c');

figure(3); % Ualpha
    plot(t,-0.03.*abs(Ualpha),'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
    ax.YLim = [-0.1 0];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Время, c');ylabel('Функция управления U\alpha_п, 1/c') 

figure(4); % UWx
        plot(t,Uw,'LineWidth',3);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14' '0.16'})
%         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
        ax.YLim = [-0.02 0.02];
%         ax.XLim = [0.00 100];
        hold all; grid on; box on;xlabel('Время, c'); ylabel('Функция управления U\omega_x, 1/c^2');        
        
figure(5); % Ualpha vs alpha  (The result shape is weird!)       
    plot(alpha,Ualpha,'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
%     ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Угол атаки \alpha_п, рад');ylabel('Функция управления U\alpha_п, 1/c') 
    
% figure(6); % Uwx vs wx
%     plot(wx,Uw,'LineWidth',3);  
%     ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%     ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
% %     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
% %     ax.YLim = [0.00 0.25];
% %     ax.XLim = [0.00 100];
%     hold all; grid on; box on; xlabel('Угловая скорость Wx, 1/c');ylabel('Функция управления Uwx, 1/c^2') 

figure(7); % dalpha/dt
    plot(t,abs(z1),'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
    ax.YLim = [0 0.5];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Время, c');ylabel('Скорость Угла атаки \alpha_п^., 1/c') 

figure(8); % phi
    plot(t,phi,'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
%     ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Время, c');ylabel('Угла phi , hfl') 

figure(9); % Wxv
        plot(t,wxv,'LineWidth',3);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14' '0.16'})
%         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
        ax.YLim = [-0.03 0.03];
%         ax.XLim = [0.00 100];
        hold all; grid on; box on;xlabel('Время, c'); ylabel('Угловая скорость \omega_x_v, 1/c');

% figure(10); % h
%     plot(t,h./1000,'LineWidth',3);  
%     ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%     ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
% %     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
% %     ax.YLim = [0.00 0.25];
% %     ax.XLim = [0.00 100];
%     hold all; grid on; box on; xlabel('Время, c');ylabel('Высота спуска, км') 

figure(11); % h vs V
    plot(V,h./1000,'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
%     ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Скорость спуска, м/с');ylabel('Высота спуска, км') 

figure(12); % h vs dV/dt
    plot(Load,h./1000,'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
%     ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Замедление');ylabel('Высота спуска, км') 
    
    
figure(13); % q
    plot(t,0.5.*rho.*V.^2,'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
%     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
%     ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('Время, с');ylabel('q, Па') 
    
    
function dydt = EQN(t,y,Aw,Bw,Aa,Ba,Aav,Bav)
t
V = y(1); theta = y(2); h = y(3); 
wx = y(4); wxv = y(5); z1 = y(6); alpha = y(7); phi = y(8);

% Constants
g0 = 3.72076; Rmars = 3396000;
r = 1.25; s = pi*r^2; L = 2; m = 576; Ix = 300; I = 443; % Iy=Iy=I;

Cxv = 1.75; Cyv = 0.3;
% Cxv = 1.75*cos(alpha)+0.3*sin(alpha); Cyv = -0.3*cos(alpha)+1.75*sin(alpha); 
dzd = 0.0005; dyd = 0.01; myf = 0.01; mzf = 0.01; 
mzn = - (0.1.*sin(alpha) + 0.05.*sin(2.*alpha)); % C1 = 0.5; % C2 = 0.1; 
epsilon = 0.05;
Cya = Cyv; % assumption
mznw = 0.01; % dmzn/dwd
Cyva = -1.75; %dyv/dalpha
mxw = 0.01; % dmx/dw

% Mars Atmospher Density Model, q, g
rho = (.699*exp(-0.00009*h))/(.1921*((-23.4-0.00222*h)+273.1));
q = 0.5*(rho)*V^2;
g = g0*Rmars^2/(Rmars+h)^2;

Cx = Cxv*cos(alpha)+Cyv*sin(alpha);
Cyn = Cyv*cos(alpha)-Cxv*sin(alpha);
Mx = Cyn*(dzd*cos(phi)+dyd*sin(phi))*q*s*L;
Myn = ((myf-Cx*dzd)*cos(phi)-(mzf+Cx*dyd)*sin(phi))*q*s*L;
dMzn= ((mzf+Cx*dyd)*cos(phi)-(myf-Cx*dzd)*sin(phi))*q*s*L;
Mxv = Mx*cos(alpha)-Myn*sin(alpha);
Mzn = mzn*q*s*L;

% EQN of CG
    dVdt = -(Cxv*q*s/m + g*sin(theta));
dthetadt = (- m*g*cos(theta)*(1-V^2/(h+Rmars))/(V*m));
    dhdt = V*sin(theta);

% Control parameters
% No Control
% Uw = 0; Ualpha = 0;

% With control
% For second order Alpha Equation: Ualpha = - sqrt(aa/ba)*alpha - sqrt(aav/bav)*dalphadt
fw = mxw.*q.*s.^2.*L./(Ix.*V./I);
Kw = fw + sqrt(fw^2 - Aw/Bw);
% Kw = sqrt(Aw/Bw);
Ka = sqrt(Aa/Ba);
Kav = sqrt(Aav/Bav);
Uw = - Kw*wx;
Ualpha = - Ka*alpha - Kav*z1; % z1 = dalpha/dt

% EQN about CG
dwxdt = epsilon*(Mx)/I + Uw;
dwxvdt= epsilon*Mxv/I + epsilon*Cya*q*s*(wx*I - wxv*I*cos(alpha))/(m*V*sin(alpha)*I);
dz1dt = Mzn/I + (wx*I - wxv*I*cos(alpha))*(wx*I*cos(alpha) - wxv*I)/(I^2*sin(alpha)^3) + ...
        + epsilon*dMzn/I + epsilon*(mznw*q*s*L^2/V - Cyva*q*s/(m*V))*z1 + Ualpha; % z1 = dalphadt
dalphadt = z1;
dphidt= (wx*I/Ix + (wx*I*cos(alpha) - wxv*I)*cos(alpha)/(I*sin(alpha)^2));

dydt = [dVdt;dthetadt;dhdt;dwxdt;dwxvdt;dz1dt;dalphadt;dphidt];
end

    
  