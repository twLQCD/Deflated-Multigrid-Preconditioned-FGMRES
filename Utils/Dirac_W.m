function A1=Dirac_W(mass,N,N2,U1,P1_plus,P2_plus,P1_minus,P2_minus);
% Constructs Wilson-Dirac operator             
v1=sparse(reshape(U1(:,1),N,N));
v2=sparse(reshape(U1(:,2),N,N));
clear U1
U=sparse(2*N2,2*N2); L=sparse(2*N2,2*N2);
p=[N,1:N-1];
for k=1:N; %loop over `time'
  e_k=sparse(N,1); e_k(k) = 1; I = sparse(diag(e_k)); T=sparse(I(:,p));
  % Upper triangular
  U=U+kron(diag(v1(k,:)),kron(T,P1_plus));
  U=U+kron(T,kron(diag(v2(:,k)),P2_plus));
  % Lower triangular
  L=L+kron(diag(v1(k,:)),kron(T,P1_minus));
  L=L+kron(T,kron(diag(v2(:,k)),P2_minus));
end
A1=(mass+2)*speye(2*N2)-(L'+U);
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
