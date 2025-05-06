% This code: For Aero and mass asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
% Notes: more hv, less ht: correct he q curve
%        mzn0 = -0.008 makes the alpha, w curves either closer or far from
%        others.
clc;
clear;
close all;

% Initial Conditions
V0 = 5500; theta0 = -(20*(pi/180)); h0 = 1e5; 
alpha10 = 8*(pi/180); wx10 = 4; thetal0= 10; % thetal (thetaL) NOT theta1 is another equation rather than theta (inclination of irectory)

% Stydy effect of different a & b on wx and alpha
mv = [366]; rv = [1.3]; Lv = [1.8]; Ixv= [135]; Izv= [186];

for j = 1:1 % for KA
y0 = [V0;theta0;h0;alpha10;wx10;thetal0];
tspan = [0 300];

[t,y] = ode45(@(t,y) eqn(t,y,j),tspan,y0);

V=y(:,1); theta=y(:,2); h=y(:,3); alpha1=y(:,4); wx1=y(:,5); thetal=y(:,6);

%-------------------------------------------------
r = rv(j); S = pi*r^2; L = Lv(j); m = mv(j); Ix = Ixv(j);Iz = Izv(j);Ixd= Ix/Iz; % Inertia I=Iy=Iz,
I = Ixv(j)/Ixv(j);
rho = (.699.*exp(-0.00009.*h))./(.1921.*((-23.4-0.00222.*h)+273.1));
q = 0.5*(rho).*V.^2;
 
mzn0 = -0.01;
omegak = sqrt(-mzn0.*q.*S.*L.*cot(alpha1)./I);
omega_xr = omegak ./ sqrt(1 - Ixd);

%-------------------------------------------------
figure(1); % alpha
    plot(t,abs(alpha1),'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
    hold all; grid on; box on; xlabel('t, c');ylabel('Угол атаки \alpha_п, Рад') 

figure(2); % Wx
        plot(t,wx1,t,omega_xr,'LineWidth',3);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        hold all; grid on; box on;xlabel('t, c'); ylabel('Угловая скорость \omega_x, 1/c');
end


function [out] = eqn(t,y,j)
V=y(1); theta=y(2); h=y(3); alpha1=y(4); wx1=y(5); thetat1=y(6);

% «Spirit» «Insight» «Schiaparelli» «Mars Polar Lander» «Mars3»
mv = [366, 576, 800]; rv = [1.3, 1.25, 1.6];Lv = [1.8, 2, 1.8];
Ixv= [135, 443, 506];
Izv= [186, 300, 768];

g0 = 3.72076; Rmars = 3396000;
r = rv(j); S = pi*r^2; L = Lv(j); m = mv(j); Ix = Ixv(j);Iz = Izv(j);Iy = Iz;Ixd= Ix/Iz; I = Iz;

% Aerodyamic
Cxv = -2;
epsilon = 0.01; 
mzn0 = -0.01;

% Mars Atmospher Density Model
rho = (.699*exp(-0.00009*h))/(.1921*((-23.4-0.00222*h)+273.1));
q = 0.5*(rho)*V^2;
w = sqrt(-mzn0*q*S*L/Iz);
g = g0*Rmars^2/(Rmars+h)^2;

  dVdt = -(Cxv*q*S/m + g*sin(theta));
dthetadt = (- m*g*cos(theta)*(1-V^2/(h+Rmars))/(V*m));
  dhdt = V*sin(theta);     
drhodh = (6291*exp(-(9*h)/100000)*dhdt)/(100000000*((213231*h)/500000000 - 4796737/100000)) + (149048469*exp(-(9*h)/100000)*dhdt)/(500000000000*((213231*h)/500000000 - 47.9674)^2);

%---------------------
theta1 = pi;
theta2 = 0;
mxf = 0.01;
Ixzd = 0.01;

mzn = -0.1;
mad = 0.09;
mxad = 0.06;

omega_k = sqrt(-mzn*q*S*L*cot(alpha1)/I);
omega_a = sqrt(Ixd^2*wx1^2/4 + omega_k^2);
omega12 = Ixd*wx1/2 + omega_a;

Mzna = mzn*q*S*L*cot(alpha1);
Fa = -Mzna/I + omega12^2/cos(alpha1)^2 + (Ixd*wx1 - omega12)*(Ixd*wx1 - 2*omega12);

dI = 0.5*(Iz - Iy);
dId = dI/I;
md = sqrt(Ixzd^2 + dId^2);
% theta3 = 0.5*asin(dId/md);
theta3 = 0.1;
%---------------------
   
% NO control
% dalpha1dt = - epsilon*(0.5*madash*w*cos(thetat1+theta1));
dalpha1dt = 4*epsilon*omega_a^2*(-omega_k*tan(alpha1)/(4*q*omega_a^2)*((0.5*V^2*drhodh) + rho*V*dVdt) - mxf*q*S*L*tan(alpha1)*omega12/(4*omega_a^2*I) + mad*omega_k^2/(2*omega_a)*cos(thetat1 + theta1) - omega12*tan(alpha1)/(4*omega_a^2)*((10 + Ixd)*wx1*omega12 - 2*(2 + Ixd)*wx1^2)*md*cos(2*thetat1 + theta3) - omega12*tan(alpha1)/(4*omega_a^2)*(tan(alpha1)^2 - 4)*omega12^2*md*cos(2*thetat1 + theta3))/Fa;
% dwx1dt = - epsilon*(mxadash*w^2*sin(thetat1+theta2)/Ixd);  
dwx1dt = -epsilon*mxad*omega_k^2*sin(thetat1 + theta2)/Ixd + epsilon*md*omega12^2*tan(alpha1)^2*cos(2*thetat1 + theta3)/Ixd;
dthetat1dt = wx1-w; %(1-0.5*Ixd)*wx2-w

out = [dVdt; dthetadt; dhdt; dalpha1dt; dwx1dt; dthetat1dt];
end
