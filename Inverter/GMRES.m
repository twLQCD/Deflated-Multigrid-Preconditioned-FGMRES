function [z,rt,mvps] = GMRES(A,b,z0,tol,nmax);
i=1; b=b(:); z=z0(:);
r=b-A*z; rnorm=norm(r); rho=rnorm; rt=[];
V(:,i)=r/rho; mvps = 0;
while ((rnorm/rho>tol)&&(i<=nmax));
  V(:,i+1)=A*V(:,i);
  mvps = mvps + 1;
  H(1:i,i)=V(:,1:i)'*V(:,i+1);
  V(:,i+1)=V(:,i+1)-V(:,1:i)*H(1:i,i);
  H(i+1,i)=norm(V(:,i+1));
  V(:,i+1)=V(:,i+1)/H(i+1,i);
  
  c = (H(1:i+1,1:i)\[rho;zeros(i,1)]);
%   z=z0+V(:,1:i)*(H(1:i+1,1:i)\[rho;zeros(i,1)]);
z=z0+V(:,1:i)*c;
  r=b-A*z; rnorm=norm(r); rt=[rt;rnorm/rho];
  i=i+1;
end
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
