 clc;
clear;
close all;

% Atmospher Characteristics
g = 3.9;  % Mars Atmospher
r = 0.9; S = pi*r^2; L = 0.9;            % Sizing r = 1.25; L = 2
m = 576; Ix = 270; Iz = 443; Ixd= Ix/Iz; % Inertia I=Iy=Iz, m = 576

% Assumptions
 Cxv = 0.04; epsilon = 0.3; % epsilon 0.02 changes the shape of convergence
mzn0 = -0.0001;  mz1 = 0.1;   my0f = 0.3;    mz0f = 0.3; 
 Cx1 = 0.9;      Cy1 = 0.9; dzdash = 0.05; dydash = 0.05;
 
% Controller
aa = 4; ba = 0.0001; aw = 5; bw = 0.0001;  % quadratic condition parameters
hv = 0.5; ht = 0.5; hh = 0.5;  % Step of speed, theta and height equations

% Initial Conditions
alpha1(1) = 0.32; wx1(1) = 0.1;
alphaz1(1) = 0.32; wxz1(1) = 0.1;
V(1) = 3500; theta(1) = -0.017; h(1) = 1e5;

k = 0;
t = 0;
for h = 100000:-6000:10000     % to get 15 steps like t 
    k = k + 1;  kk(k) = h;
    tt(k) = t; t = t + 10;  % time interval 0:10:400
    
% Using Euler           
     [rho,~] = marsatmoshper(h); % Mars Atmospher Density Model
           q = 0.5*rho*V(k)^2; qv(k) = q;
           w = sqrt(-mzn0*q*S*L/Iz);
       
      V(k+1) = V(k) - hv*(Cxv*q*S/m + g*sin(theta(k)));
  theta(k+1) = theta(k) - ht*g*sin(theta(k))/V(k);
%       h(k+1) = h(k) + hh*V(k)*sin(theta(k));     

     ma1 = (-my0f + Cx1*dzdash)*(w^2/mz1); %-(my0f/mz1)*w^2 + (Cx1/mz1)*w^2*dzdash;
     ma2 = (-mz0f - Cx1*dydash)*(w^2/mz1);
  madash = sqrt(ma1^2 + ma2^2)/w^2;
   f1(k) = 0.5*madash*w;
  
    mxa1 = -(Cy1/mz1)*w^2*dydash;
    mxa2 = -(Cy1/mz1)*w^2*dzdash;
 mxadash = sqrt(mxa1^2 + mxa2^2)/w^2;
   f2(k) = mxadash*w^2;
 
  theta1 = asin(ma1/(madash*w^2));
  theta2 = asin(-mxa1/(mxadash*w^2));
 
  % Control in the given system
alpha1(k+1) = alpha1(k) - 0.008*epsilon*(0.5*madash*w*cos(theta(k)+theta1) + sqrt(aa*ba)*alpha1(k)/ba);          
   wx1(k+1) = wx1(k) - 0.01*epsilon*(mxadash*w^2*sin(theta(k)+theta2)/Ixd + sqrt(aw*bw)*wx1(k)/bw);          
% alpha1(k+1) = alpha1(k)*(1 - 25*epsilon*(0.5*madash*w*cos(theta(k)+theta1) + sqrt(aa*ba)));          
%    wx1(k+1) = wx1(k)*(1 - 30*epsilon*(mxadash*w^2*sin(theta(k)+theta2)/Ixd + sqrt(aw*bw)));          

Ualpha(k) = - sqrt(aa*ba)*alpha1(k)/ba;          
   Uwx(k) = - sqrt(aw*bw)*wx1(k)/bw;  

end


%% Charting
figure(1);
        plot([0 tt],alpha1,'c','linewidth',4); hold on;
        grid on; box on; xlabel('t(n) [C]'); ylabel('\alpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        axis([0 150 0 0.4])
        xticks([0 25 50 75 100 125 150])
figure(2)
        plot([0 tt],wx1,'c','linewidth',4); hold on;
        grid on; box on; xlabel('t(n) [C]'); ylabel('\omega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        axis([0 150 0 0.12])
        xticks([0 25 50 75 100 125 150])  
%         yticklabels({0, 0.02, 0.04, 0.06, 0.08, '0.10', 0.12})  
        yticks([0 0.02 0.04 0.06 0.08 0.10 0.12])  
figure(3);
        plot(tt,Ualpha./40,'b','linewidth',4);
        grid on; box on; xlabel('t(n) [C]'); ylabel('\alpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
hold on
        plot(tt,Uwx./40,'g','linewidth',4);
        grid on; box on; xlabel('t(n) [C]'); ylabel('\omega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        xticks([0 25 50 75 100 125 150]) 
        yticks([-2.0 -1.5 -1.0 -0.5 0])
        yticklabels({'-2.0', '-1.5', '-1.0', '-0.5', '0'})
        axis([0 150 -2 0])        
        legend('U\alpha',...
               'U\omega_x')        