% This code: For Aero and mass asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
% Notes: more hv, less ht: correct he q curve
%        mzn0 = -0.008 makes the alpha, w curves either closer or far from
%        others.
clc;
clear;
close all;

% Initial Conditions
V0 = 3500; theta0 = -0.017; gamma0 = 0.1; h0 = 1e5;
alpha10 = 0.32; wx10 = 0.14; thetal0 = 0; % thetal (thetaL) NOT theta1 is another equation rather than theta (inclination of irectory)
aa = 1; ba = 0.0005; aw = 1; bw = 0.0005;  % quadratic condition parameters

y0 = [V0; theta0; gamma0; h0; alpha10; wx10; thetal0];
options = odeset('RelTol',1e-2,'Stats','on','OutputFcn',@odeplot)

[t,y] = ode23(@eqn,[0:8:300],y0,options);
y

figure(1); plot(t,y(:,1),'linewidth',2); ylabel('Velocity'); grid on;
figure(2); plot(t,y(:,2)); ylabel('Theta'); grid on;
figure(3); plot(t,y(:,3),'linewidth',2); ylabel('Gamma'); grid on;
figure(4); plot(t,y(:,4)); ylabel('H'); grid on;
figure(5); stairs(t,y(:,5)); ylabel('Alpha'); grid on;
figure(6); stairs(t,y(:,6)); ylabel('Wx'); grid on;

% n = 1:10:300;  % discretization counter
% for tn = 1 : length(y)
%     if n(ismember(n,tn)) == tn
%        Ualpha(tn) = - sqrt(aa*ba)*y(tn,5)/ba;          
%        Uwx(tn) = - sqrt(aw*bw)*y(tn,6)/bw;  % 
%     else
%        Ualpha(tn) = Ualpha(tn-1);
%        Uwx(tn) =  Uwx(tn-1);    
%     end
% end
% 
% figure(7); stairs(t,Ualpha); ylabel('Alpha'); grid on;
% figure(8); stairs(t,Uwx); ylabel('Wx'); grid on;

% figure(9); % discret; alpha
%         plot(t(1:end),alpha1(1:end-1),'LineWidth',2);hold all
% %         stairs(1:16:length(alphaz),abs(alphaz(1:1:10)),'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'})
% %         yticks([0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.345]);
% %         ax.YLim = [0.00 0.35];
% %         ax.XLim = [0.00 150];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\alpha [Рад]') 
% 
% figure(10); % discret wx euler and z
%         plot(t(1:end),wx1(1:end-1),'LineWidth',2); hold all
% %         stairs(1:16:length(wxz),abs(wxz(1:1:10)),'LineWidth',2); 
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
% %         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
% %         ax.YLim = [0.00 0.14];
% %         ax.XLim = [0.00 150];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\omega [1/C]');  
% 
% figure(11); % Ualphaz
%         stairs(kk(1:end),Ualpha(1:end),'linewidth',2); hold all
%         grid on; box on; xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% % %         xticks([0 50 100 150 200 250 300]) 
% %         yticks([-0.04 -0.03 -0.02 -0.01 0])
% %         yticklabels({'-0.04' '-0.03' '-0.02' '-0.01' '0'})
% %         %         axis([0 300 -0.1 0])  
% %         ax.YLim = [-0.04 0];
% %         ax.XLim = [0 150];
% figure(12); % Uwxz
%         stairs(kk(1:end),Uwx(1:end),'linewidth',2); hold all
%         grid on; box on; xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% % %         xticks([0 25 50 75 100 125 150]) 
% %         yticks([-0.05 -0.04 -0.03 -0.02 -0.01 0])
% %         yticklabels({'-0.05' '-0.04' '-0.03' '-0.02' '-0.01' '0'})
% % %         axis([0 150 -2 0]) 
% %         ax.YLim = [-0.05 0];
% %         ax.XLim = [0 150];
%  

