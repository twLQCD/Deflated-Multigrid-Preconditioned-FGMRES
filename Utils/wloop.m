function wlp = wloop(N,N2,kp,km,theta1);
% Computation of Wilson loops
q=1:N2;
for mu=1:2;
  tleg1(q,mu)=theta1(q,mu);
  tleg2(q,mu)=theta1(kp(q,mu),mu);
  tleg3(q,mu)=theta1(kp(kp(q,mu),mu),mu);
  tleg4(q,mu)=theta1(kp(kp(kp(q,mu),mu),mu),mu);
  tleg5(q,mu)=theta1(kp(kp(kp(kp(q,mu),mu),mu),mu),mu);
  tleg6(q,mu)=theta1(kp(kp(kp(kp(kp(q,mu),mu),mu),mu),mu),mu);
  tleg7(q,mu)=theta1(kp(kp(kp(kp(kp(kp(q,mu),mu),mu),mu),mu),mu),mu);
  tleg8(q,mu)=theta1(kp(kp(kp(kp(kp(kp(kp(q,mu),mu),mu),mu),mu),mu),mu),mu);
  tleg9(q,mu)=theta1(kp(kp(kp(kp(kp(kp(kp(kp(q,mu),mu),mu),mu),mu),mu),mu),mu),mu);
  tleg10(q,mu)=theta1(kp(kp(kp(kp(kp(kp(kp(kp(kp(q,mu),mu),mu),mu),mu),mu),mu),mu),mu),mu);
  wleg1(q,mu)=tleg1(q,mu);
  wleg2(q,mu)=wleg1(q,mu)+tleg2(q,mu);
  wleg3(q,mu)=wleg2(q,mu)+tleg3(q,mu);
  wleg4(q,mu)=wleg3(q,mu)+tleg4(q,mu);
  wleg5(q,mu)=wleg4(q,mu)+tleg5(q,mu);
  wleg6(q,mu)=wleg5(q,mu)+tleg6(q,mu);
  wleg7(q,mu)=wleg6(q,mu)+tleg7(q,mu);
  wleg8(q,mu)=wleg7(q,mu)+tleg8(q,mu);
  wleg9(q,mu)=wleg8(q,mu)+tleg9(q,mu);
  wleg10(q,mu)=wleg9(q,mu)+tleg10(q,mu);
end
wtheta=wleg1(q,1)+wleg1(kp(q,1),2)-wleg1(q,2)-wleg1(kp(q,2),1);
wlp(1)=mean(cos(wtheta));
wtheta=wleg2(q,1)+wleg2(kp(q,1),2)-wleg2(q,2)-wleg2(kp(q,2),1);
wlp(2)=mean(cos(wtheta));
wtheta=wleg3(q,1)+wleg3(kp(q,1),2)-wleg3(q,2)-wleg3(kp(q,2),1);
wlp(3)=mean(cos(wtheta));
wtheta=wleg4(q,1)+wleg4(kp(q,1),2)-wleg4(q,2)-wleg4(kp(q,2),1);
wlp(4)=mean(cos(wtheta));
wtheta=wleg5(q,1)+wleg5(kp(q,1),2)-wleg5(q,2)-wleg5(kp(q,2),1);
wlp(5)=mean(cos(wtheta));
wtheta=wleg6(q,1)+wleg6(kp(q,1),2)-wleg6(q,2)-wleg6(kp(q,2),1);
wlp(6)=mean(cos(wtheta));
wtheta=wleg7(q,1)+wleg7(kp(q,1),2)-wleg7(q,2)-wleg7(kp(q,2),1);
wlp(7)=mean(cos(wtheta));
wtheta=wleg8(q,1)+wleg8(kp(q,1),2)-wleg8(q,2)-wleg8(kp(q,2),1);
wlp(8)=mean(cos(wtheta));
wtheta=wleg9(q,1)+wleg9(kp(q,1),2)-wleg9(q,2)-wleg9(kp(q,2),1);
wlp(9)=mean(cos(wtheta));
wtheta=wleg10(q,1)+wleg10(kp(q,1),2)-wleg10(q,2)-wleg10(kp(q,2),1);
wlp(10)=mean(cos(wtheta));
wlp=[wlp(1);wlp(2);wlp(3);wlp(4);wlp(5);wlp(6);wlp(7);wlp(8);wlp(9);wlp(10)];
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
