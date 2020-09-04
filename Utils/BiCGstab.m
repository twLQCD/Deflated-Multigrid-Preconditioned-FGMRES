function [z,rr,mvps]=BiCGstab(A,b,z0,tol,nmax)
n=1; b=b(:); z=z0(:); y0=rand(size(b)); y0=y0/norm(y0);
r=b-A*z; s=r; psi(1)=0; delta = y0'*r;
rnorm=norm(r); rho=rnorm; rr=[]; mvp = 0; mvps = [];
while ((rnorm>tol)&(n<=nmax));
  As_n=A*s; 
  mvp = mvp + 1;
  phi(n)=y0'*As_n/delta;
  omega=1/phi(n);
  w=r-omega*As_n;
  Aw=A*w;
  mvp = mvp + 1;
  mvps = [mvps,mvp];
  chi=(Aw'*w)/(Aw'*Aw);
  zelos(n)=1/chi;
  r=w-chi*Aw;
  z=z+omega*s+chi*w;
  dtmp=y0'*r;
  psi(n+1)=-omega*dtmp/(delta*chi);
  delta=dtmp;
  s=r-psi(n+1)*(s-chi*As_n);
  n=n+1;
  rnorm=norm(r);
  rr=[rr,rnorm];
end
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
