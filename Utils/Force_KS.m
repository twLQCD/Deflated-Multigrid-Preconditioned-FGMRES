function [pf,pg]=Force_KS(N2,kp,km,theta1,eta,chi,U1);
sgn=(-1).^(0:N2-1); % For square lattice it is easy!
eta=reshape(eta,1,N2);
chi=reshape(chi,1,N2);
% U^+_shift eta & chi
eta_s1=       U1(:,1)'.*eta;
eta_s2=-sgn.*(U1(:,2)'.*eta);
chi_s1=       U1(:,1)'.*chi;
chi_s2= sgn.*(U1(:,2)'.*chi);
% Fermion force
p1=-conj(eta(kp(:,1))).*chi_s1+conj(chi(kp(:,1))).*eta_s1;
p2=-conj(eta(kp(:,2))).*chi_s2+conj(chi(kp(:,2))).*eta_s2;
pf(:,1)=2.*real(sqrt(-1)*p1)';
pf(:,2)=2.*real(sqrt(-1)*p2)';
% Gauge boson force
for mu=1:2;
  if mu==1,nu=2; else nu=1; end
  alpha1=theta1(kp(:,mu),nu)-theta1(kp(:,nu),mu)-theta1(:,nu);
  alpha2=theta1(km(:,nu),nu)-theta1(km(:,nu),mu)-theta1(kp(km(:,nu),mu),nu);
  pg(:,mu)=sin(theta1(:,mu)+alpha1)+sin(theta1(:,mu)+alpha2);
end
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
