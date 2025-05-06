% This code: For Aero and mass asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
% Notes: more hv, less ht: correct he q curve
%        mzn0 = -0.008 makes the alpha, w curves either closer or far from
%        others.
clc;
clear;
close all;

% «Spirit» «Insight» «Schiaparelli» «Mars Polar Lander» «Mars3»
mv = [174, 366, 577, 576, 800];
rv = [1.15, 1.3, 1.2, 1.25, 1.6];
Lv = [1.5, 1.8, 1.8, 2, 1.8];
Ixv= [90, 186, 250, 270, 768];
Izv= [80, 135, 195, 443, 506];
BA = [0.001]; BW = [0.001];

% mv = [174]; rv = [1.15]; Lv = [1.5]; Ixv= [90]; Izv= [80]; % activate it for engine study

for j = 1:length(mv)
% Atmospher Characteristics
g0 = 3.9; rho = 0.019;  % Mars Atmospher
r = rv(j); S = pi*r^2; L = Lv(j);            % Sizing r = 1.25; L = 2
m = mv(j); Ix = Ixv(j); Iz = Izv(j); Ixd= Ix/Iz; % Inertia I=Iy=Iz, m = 576
Rmars = 3396000;

% Assumptions
 Cxv = 0.04; epsilon = 0.01; % epsilon 0.02 changes the shape of convergence
mzn0 = -0.01;  mz1 = -0.3;   my0f = 0.3;    mz0f = 0.3; 
 Cx1 = 0.9;     Cy1 = 0.9; dzdash = 0.05; dydash = 0.05;
 Cya = 0.9;
 
% Controller
aa = 1; ba = 0.005; aw = 1; bw = 0.005;  % quadratic condition parameters
hv = 3; ht = 0.01; hh = 0.4; hg = 0.01; % hv = 15; ht = 0.00001; hh = 12; hg = 0.1;
% high hv, small ht: forms the right q curve 

% Initial Conditions
alpha1(1) = 0.32; wx1(1) = 0.14; thetal(1) = 0; % thetal (thetaL) NOT theta1 is another equation rather than theta (inclination of irectory)
alpha2(1) = 0.32; wx2(1) = 0.14; thetat2(1) = 0;
alphaz(1) = 0.32; wxz(1) = 0.14;
V(1) = 5000; theta(1) = -(12*(pi/180)); gamma(1)= 0.1; h(1) = 1e5; % V(1) = 3500; theta(1) = -0.017

t = 0;
% n = 1:2:300;        % discretization counter
for k = 1:1:300 % integration counter
    t = t + 1;
    tt(t) = t-1;
%     nn(t) = n(t); 
    
% Using Euler           
     [Rho,~] = marsatmoshper(h(t)); % Mars Atmospher Density Model
     rho = Rho;
           q = 0.5*(rho)*V(t)^2; RHO(t) = rho; qv(t)=q;
           w = sqrt(-mzn0*q*S*L/Iz);
           g = g0*(Rmars/(Rmars+h(t)))^2;
           
      V(t+1) = V(t) - hv*(Cxv*q*S/m + g*sin(theta(t)));
  theta(t+1) = theta(t) + ht*(- m*g*cos(theta(t))*(1-V(t)^2/(h(t)+Rmars))/(V(t)*m));
  gamma(t+1) = gamma(t) + hg*(wx1(t)*cos(alpha1(t))+tan(theta(t))*(Cya*q*S*sin(gamma(t))/(m*V(t))));
      h(t+1) = h(t) + hh*V(t)*sin(theta(t));     

     ma1 = (-my0f + Cx1*dzdash)*(w^2/mz1); %-(my0f/mz1)*w^2 + (Cx1/mz1)*w^2*dzdash;
     ma2 = (-mz0f - Cx1*dydash)*(w^2/mz1);
  madash = sqrt(ma1^2 + ma2^2)/w^2;
  
    mxa1 = -(Cy1/mz1)*w^2*dydash;
    mxa2 = -(Cy1/mz1)*w^2*dzdash;
 mxadash = sqrt(mxa1^2 + mxa2^2)/w^2;
 
  theta1 = asin(ma1/(madash*w^2));
  theta2 = asin(-mxa1/(mxadash*w^2));


%   % Control in the given system
% if n(ismember(n,t)) == t
   Ualpha(t) = - 0.5*epsilon*sqrt(aa*ba)*alpha1(t)/ba;          
   Uwx(t) = - 0.5*epsilon*sqrt(aw*bw)*wx1(t)/bw;  % 
% else
%    Ualpha(t) = Ualpha(t-1);
%    Uwx(t) =  Uwx(t-1);    
% end
 
alpha1(t+1) = alpha1(t) - 0.5*epsilon*(0.5*madash*w*cos(thetal(t)+theta1)) +  Ualpha(t);          
   wx1(t+1) = wx1(t) - 0.5*epsilon*(mxadash*w^2*sin(thetal(t)+theta2)/Ixd) + Uwx(t);          
thetal(t+1) = thetal(t) - 0.5*epsilon*(1-0.5*Ixd)*wx1(t)-w;
   
% NO control
alpha2(t+1) = alpha2(t) - (0.5*madash*w*cos(2.5*thetat2(t)+theta1));          
   wx2(t+1) = wx2(t) - (mxadash*w^2*sin(thetat2(t)+theta2)/Ixd);   
thetat2(t+1) = thetat2(t) - 0.5*epsilon*(1-0.5*Ixd)*wx2(t)-w;
% Discontinous Using Inverse Z:  
          Ka = sqrt(aa*ba); 
          Kw = sqrt(aw*bw); 
alphaz(t+1) = alphaz(t)*(-Ka)^t;
   wxz(t+1) = wxz(t)*(-Kw)^t;  

% Control laws in body coordinate system OXYZ   
   Ux(t) = Uwx(t);
   Uy(t) = sin(gamma(t))*Ualpha(t);
   Uz(t) = cos(gamma(t))*Ualpha(t);
   
   
% MAximum deceleration
a1 = 0.699; a2 = 0.00009; a3 = 47.967; a4 = 0.000426; a5 = a2*a3-a4; a6 = a2*a4;
Ve = V(1); hs = h(t);
a_max(t) = -(Ve^2*sin(theta(t))/(2*exp(1)))*((a5-a6*hs)/(a3-a4*hs));
dvdt(t) = Cxv*q*S/m + g*sin(theta(t)); 
n(t) = a_max(t)./g;
nv(t) = dvdt(t)./g;
% % Engines ability
%    mux(k) = Ux(k)*Ix;                        % moment about ox
%    muy(k) = -(Uy(k))*(2*Iz*w); % moment about oy
%    muz(k) = -(Uz(k))*(2*Iz*w); % moment about oz   
end

figure(111);
        plot(h(1:end-1),nv,'LineWidth',2); 
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        hold all; grid on; box on; xlabel('t [C]');ylabel('\Amax [m/s^2]') 
figure(112);
        plot(h(1:end-1),RHO,'LineWidth',2); 
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        hold all; grid on; box on; xlabel('t [C]');ylabel('\RHO ') 


%    Mux = max(abs(mux)); %    Mux = trapz(mux);  % moment [N.m]
%    Muy = max(abs(muy)); %    Muy = trapz(muy);  % moment [N.m]
%    Muz = max(abs(muz)); %    Muz = trapz(muz);  % moment [N.m]
%    
%    Fux = Mux/(4*r) % N
%    Fuy = Muy/(4*r) % N
%    Fuz = Muz/(4*r) % N

% Charting
figure(1); % alpha
        plot(tt,(alpha1(1:end-1)),'LineWidth',2); 
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'})
        yticks([0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.345]);
        ax.YLim = [0.00 0.35];
        ax.XLim = [0.00 tt(end)];
        hold all; grid on; box on; xlabel('t(n) [C]');ylabel('\alpha [Рад]') 
       
figure(2); % Wx
        plot(tt,wx1(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
        yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
%         ax.YLim = [0.00 0.14];
        hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\omega [1/C]');
       
% figure(22); % ThetaL with control
%         plot(tt,thetal(1:end-1),'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
% %         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
% %         ax.YLim = [0.00 0.14];
%         hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\thetaL [1/C]');        
%         
% figure(3); % Ualpha
%         plot(tt,Ualpha,'linewidth',2); hold all
%         grid on; box on; xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         xticks([0 50 100 150 200 250 300]) 
%         yticks([-0.04 -0.03 -0.02 -0.01 0])
%         yticklabels({'-0.04' '-0.03' '-0.02' '-0.01' '0'})
%         %         axis([0 300 -0.1 0])  
%         ax.YLim = [-0.04 0];
% figure(4); % Uwx
%         plot(tt,Uwx,'linewidth',2); hold all
%         grid on; box on; xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         xticks([0 25 50 75 100 125 150]) 
%         yticks([-0.05 -0.04 -0.03 -0.02 -0.01 0])
%         yticklabels({'-0.05' '-0.04' '-0.03' '-0.02' '-0.01' '0'})
% %         axis([0 150 -2 0]) 
%         ax.YLim = [-0.05 0];
% 
% figure(5); %% Dynamic pressure 
%         plot(tt,qv./13,'linewidth',2); ylabel('q')
%         grid on; box on; xlabel('t(n) [C]'); hold all
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         xticks([0 50 100 150 200 250 300]) 
%         yticks([0 100 200 300 400 500 600 700])
%         yticklabels({'0', '100', '200', '300', '400', '500', '600', '700'})
%         axis([0 300 0 700]) 
% 
% % without control
figure(6); % alpha2 no control
        plot(tt,abs(alpha2(1:end-1)),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'-1.0' '-0.75' '-0.50' '-0.25' '0.0' '0.25' '0.50' '0.75' '1.0'})
        yticks([-1.0 -0.75 -0.50 -0.25 0.0 0.25 0.50 0.75 1.0]);
%         ax.YLim = [0.00 0.35];
%         ax.XLim = [0.00 tt(end)];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\alpha [Рад]')     
figure(7); % wx2 no control
        plot(tt,wx2(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'-0.4' '-0.2' '0.0' '0.2' '0.4' '0.6'})
        yticks([-0.4 -0.2 0.0 0.2 0.4 0.6]);
%         ax.YLim = [0.00 0.14];
        hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\omegax [1/C]');
% figure(77); % theta (fast phase) no control
%         plot(tt,thetat2(1:end-1),'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'-0.4' '-0.2' '0.0' '0.2' '0.4' '0.6'})
% %         yticks([-0.4 -0.2 0.0 0.2 0.4 0.6]);
% %         ax.YLim = [0.00 0.14];
%         hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\theta [1/C]');        
% 
figure(8); % h
        plot(tt,h(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
        yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
        ax.YLim = [0.00 100000];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('h');  
%         
figure(9); % discret; alpha
        plot(tt(1:4:end),alpha1(1:4:end-1),'LineWidth',2);hold all
%         stairs(1:16:length(alphaz),abs(alphaz(1:1:10)),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'})
        yticks([0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.345]);
        ax.YLim = [0.00 0.35];
        ax.XLim = [0.00 100];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\alpha [Рад]') 
figure(10); % discret wx euler and z
        plot(tt(1:4:end),wx1(1:4:end-1),'LineWidth',2); hold all
%         stairs(1:16:length(wxz),abs(wxz(1:1:10)),'LineWidth',2); 
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
        yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
        ax.YLim = [0.00 0.14];
        ax.XLim = [0.00 100];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\omega [1/C]');  
figure(11); % Ualphaz
        stairs(tt(1:4:end),Ualpha(1:4:end),'linewidth',2); hold all
        grid on; box on; xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 50 100 150 200 250 300]) 
        yticks([-0.030 -0.025 -0.020 -0.015 -0.010 -0.005 0])
        yticklabels({'-0.030' '-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
        %         axis([0 300 -0.1 0])  
        ax.YLim = [-0.03 0];
        ax.XLim = [0 100];
figure(12); % Uwxz
        stairs(tt(1:4:end),Uwx(1:4:end),'linewidth',2); hold all
        grid on; box on; xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 25 50 75 100 125 150]) 
        yticks([-0.010 -0.008 -0.006 -0.004 -0.002 0])
        yticklabels({'-0.010' '-0.008' '-0.006' '-0.004' '-0.002' '0'})
%         axis([0 150 -2 0]) 
        ax.YLim = [-0.01 0];
        ax.XLim = [0 100];
 
% figure(13); % Gamma
%         plot(tt,gamma(1:end-1),'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
% %         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
% %         ax.YLim = [0.00 100000];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('Gamma');  
% figure(14); % Theta "Tang"
%         plot(tt,theta(1:end-1),'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
% %         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
% %         ax.YLim = [0.00 100000];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('Tang');  
figure(15); % V
        plot(tt,V(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
%         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
%         ax.YLim = [0.00 100000];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('V'); 

% % figure(16); All Control functions
% % subplot(321); % Ualpha
% %         plot(tt,Ualpha,'linewidth',2); hold all
% %         grid on; box on; %xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
% %         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% %         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% % subplot(323); % Uwx
% %         plot(tt,Uwx,'linewidth',2); hold all
% %         grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
% %         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% %         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% % subplot(322); % Ualpha
% %         plot(tt,Ux,'linewidth',2); hold all
% %         grid on; box on; %xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
% %         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% %         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% % subplot(324); % Uwx
% %         plot(tt,Uy,'linewidth',2); hold all
% %         grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
% %         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% %         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% % subplot(326); % Uwx
% %         plot(tt,Uz,'linewidth',2); hold all
% %         grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
% %         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% %         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;

figure(17); % control functions in OXYZ system, make uw function weight coefficient equal to 1
subplot(221); % Ux
        plot(tt,Ux,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticks([-0.01 0])
        yticks([-0.010 -0.0075 -0.0050 -0.0025 0]);
        yticklabels({'-0.01' '-0.0075' '-0.0050' '-0.0025' '0'})
        ax.YLim = [-0.01 0];
        ax.XLim = [0 100];
subplot(222); % Uy
        plot(tt,Uy,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticks([-0.003 -0.002 -0.001 0]); ax.YLim = [-0.003 0];
        yticklabels({'-0.003' '-0.002' '-0.001' '0'})  
        ax.XLim = [0 100];
subplot(223); % Uz
        plot(tt,Uz,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;       
        yticks([-0.03 -0.02 -0.01 0]); ax.YLim = [-0.03 0];
        yticklabels({'-0.03' '-0.02' '-0.01' '0'})  
        ax.XLim = [0 100];

% figure(178); % discrete control functions in OXYZ system
% subplot(221); % Ux
%         stairs(tt(1:2:end),Ux(1:2:end),'linewidth',2); hold all
%         grid on; box on; %xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticks([-0.05 -0.04 -0.03 -0.02 -0.01 0])
%         yticklabels({'-0.05' '-0.04' '-0.03' '-0.02' '-0.01' '0'})
%         ax.YLim = [-0.05 0]; 
% subplot(222); % Uy
%         stairs(tt(1:2:end),Uy(1:2:end),'linewidth',2); hold all
%         grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticks([-0.004 -0.003 -0.002 -0.001 0]); ax.YLim = [-0.004 0];
%         yticklabels({'-0.004' '-0.003' '-0.002' '-0.001' '0'})        
% subplot(223); % Uz
%         stairs(tt(1:2:end),Uz(1:2:end),'linewidth',2); hold all
%         grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;       
%         yticks([-0.04 -0.03 -0.02 -0.01 0]); ax.YLim = [-0.04 0];
%         yticklabels({'-0.04' '-0.03' '-0.02' '-0.01' '0'})         
        
% figure(19); % For Studying Engines Efffect, control functions in OXYZ system 
%         plot(nn,Ux,'linewidth',3); hold all
%         plot(nn,Uy,'linewidth',3); hold all
%         plot(nn,Uz,'linewidth',3); hold all
%         grid on; box on; 
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticks([-0.025 -0.020 -0.015 -0.010 -0.005 0])
% %         yticklabels({'-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
%         ax.YLim = [-0.08 0];
%         
%         st=2;
% figure(20); % For Studying Engines Efffect, control functions in OXYZ system 
%         stairs(nn(1:st:end),Ux(1:st:end),'linewidth',3); hold all
%         stairs(nn(1:st:end),Uy(1:st:end),'linewidth',3); hold all
%         stairs(nn(1:1:end),Uz(1:1:end),'linewidth',3); hold all
%         grid on; box on; 
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticks([-0.025 -0.020 -0.015 -0.010 -0.005 0])
% %         yticklabels({'-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
%         ax.YLim = [-0.08 0];

end