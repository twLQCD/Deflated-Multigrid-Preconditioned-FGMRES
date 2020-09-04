function [pf,pg,eta_s1]=Force_W(N2,kp,km,theta1,eta,chi,U1,P1_plus,P2_plus,P1_minus,P2_minus);
eta=reshape(eta,2,N2);
chi=reshape(chi,2,N2);
% U^+_Pmu_plus_shift eta and U^+_Pmu_minus_shift chi
eta_s1=[U1(:,1),U1(:,1)]'.*(P1_plus*eta);
eta_s2=[U1(:,2),U1(:,2)]'.*(P2_plus*eta);
chi_s1=[U1(:,1),U1(:,1)]'.*(P1_minus*chi);
chi_s2=[U1(:,2),U1(:,2)]'.*(P2_minus*chi);

sU12 = size([U1(:,1),U1(:,1)]);
setas1 = size(eta_s1);
% Fermion force
for k=1:2;
  p1(k,:)=conj(eta(k,kp(:,1))).*chi_s1(k,:)+conj(chi(k,kp(:,1))).*eta_s1(k,:);
  p2(k,:)=conj(eta(k,kp(:,2))).*chi_s2(k,:)+conj(chi(k,kp(:,2))).*eta_s2(k,:);
end
pf(:,1)=2.*real(sqrt(-1)*p1(1,:)+p1(2,:)*sqrt(-1))';
pf(:,2)=2.*real(sqrt(-1)*p2(1,:)+p2(2,:)*sqrt(-1))';
% Gauge boson force
for mu=1:2;
  if mu==1,nu=2; else nu=1; end
  alpha1=theta1(kp(:,mu),nu)-theta1(kp(:,nu),mu)-theta1(:,nu);
  alpha2=theta1(km(:,nu),nu)-theta1(km(:,nu),mu)-theta1(kp(km(:,nu),mu),nu);
  pg(:,mu)=sin(theta1(:,mu)+alpha1)+sin(theta1(:,mu)+alpha2);
end
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
