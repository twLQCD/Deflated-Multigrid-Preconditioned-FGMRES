function [A] = Diracop(N,theta2,mass)

N2 = N^2;

P1_plus = [1,1;1,1]/2; P1_minus = [1,-1;-1,1]/2;
P2_plus = [1,-1i;1i,1]/2; P2_minus = [1,1i;-1i,1]/2;

U1=cos(theta2)+sqrt(-1)*sin(theta2);

A = Dirac_W(mass,N,N2,U1,P1_plus,P2_plus,P1_minus,P2_minus);