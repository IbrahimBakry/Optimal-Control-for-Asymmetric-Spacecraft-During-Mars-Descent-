function dydt = eqn(t,y)
V = y(1); theta = y(2); gamma = y(3); h = y(4); alpha1 = y(5); wx1 = y(6); thetal = y(7);

mv = [174]; rv = [1.15]; Lv = [1.5]; Ixv= [90]; Izv= [80]; g0 = 3.9; rho = 0.019;  % Mars Atmospher
r = rv; S = pi*r^2; L = Lv; m = mv; Ix = Ixv; Iz = Izv; Ixd= Ix/Iz; Rmars = 3396000;

% Assumptions
 Cxv = 0.04; epsilon = 0.01; % epsilon 0.02 changes the shape of convergence
mzn0 = -0.01;  mz1 = -0.3;   my0f = 0.3;    mz0f = 0.3; 
 Cx1 = 0.9;      Cy1 = 0.9; dzdash = 0.05; dydash = 0.05; Cya = 0.05;
 
% Controller
aa = 1; ba = 0.05; aw = 1; bw = 0.05;  % quadratic condition parameters
  
% Using Euler           
     [rho,~] = marsatmoshper(h); % Mars Atmospher Density Model
           q = 0.5*(rho)*V^2;
           w = sqrt(-mzn0*q*S*L/Iz);
%            g = g0*(Rmars/Rmars+h)^2;
          g = g0;
      
      dVdt = -Cxv*q*S/m - g*sin(theta);
  dthetadt = (Cya*cos(gamma)*q*S - m*g*cos(theta)*(1-V^2/(h+Rmars))/(V*m));
  dgammadt = (wx1*cos(alpha1)+tan(theta)*(Cya*q*S*sin(gamma)/(m*V)));
      dhdt = V*sin(theta);     

     ma1 = (-my0f + Cx1*dzdash)*(w^2/mz1); %-(my0f/mz1)*w^2 + (Cx1/mz1)*w^2*dzdash;
     ma2 = (-mz0f - Cx1*dydash)*(w^2/mz1);
  madash = sqrt(ma1^2 + ma2^2)/w^2;
  
    mxa1 = -(Cy1/mz1)*w^2*dydash;
    mxa2 = -(Cy1/mz1)*w^2*dzdash;
 mxadash = sqrt(mxa1^2 + mxa2^2)/w^2;
 
  theta1 = asin(ma1/(madash*w^2));
  theta2 = asin(-mxa1/(mxadash*w^2));

   Ualpha = - sqrt(aa*ba)*alpha1/ba;          
   Uwx = - sqrt(aw*bw)*wx1/bw; 

dalpha1dt = epsilon*(0.5*madash*w*cos(thetal+theta1) + Ualpha) ;          
   dwx1dt = epsilon*(mxadash*w^2*sin(thetal+theta2)/Ixd + Uwx) ;          
dthetaldt = epsilon*(1-0.5*Ixd)*wx1-w;

dydt = [dVdt; dthetadt; dgammadt; dhdt; dalpha1dt; dwx1dt; dthetaldt];
end
