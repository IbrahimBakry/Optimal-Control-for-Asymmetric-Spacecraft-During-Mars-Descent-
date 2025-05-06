

mzn = - (C1*sin(alpha) + C2*sin(2*alpha));
mz1 = mzn;% when dx = 0;

p1 = 1.5; p2 = 1.5; p3 = 0.5; p4 = 0.5;
p3k = 1; p5k = 0.5; p4k = 1; p6k = 0.5;

myf01 = (-mz1/(1-Ixd))*(p3/p1)*Ixy; % Мои
mzf01 = (-mz1/(1-Ixd))*(p4/p2)*Ixz; % Мои

myf02 = (-mz1/(1-Ixd))*(p3k/p5k)*Ixy; % Куркина
mzf02 = (-mz1/(1-Ixd))*(p4k/p6k)*Ixz; % Куркина
