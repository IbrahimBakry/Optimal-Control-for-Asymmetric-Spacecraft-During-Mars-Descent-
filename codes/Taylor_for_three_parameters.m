% Taylor for Three parameters
syms x y z a1 a2 b1 b2 c1 c2
f = sqrt((a1*x+b1*y+c1*z)^2+(a2*x+b2*y+c2*z)^2);
T = taylor(f,[x,y,z],'ExpansionPoint',[1,1,1], ...
            'OrderMode','Absolute','Order',2);
