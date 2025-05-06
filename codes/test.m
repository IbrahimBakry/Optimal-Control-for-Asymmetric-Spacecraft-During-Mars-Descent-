
syms A B C aw aa bw ba fw

eqn = [2*A*fw - (1/bw)*A^2 - (1/(4*ba))*C^2 == -aw,...
       -(1/ba)*B^2 - (1/(4*bw))*C^2 == -aa,...
       -(1/bw)*A - (1/ba)*B == -fw];
   
   solve(eqn,[A B C])

