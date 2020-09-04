function A1=Dirac_KS(mass,N,N2,U1);
% Constructs Kogut-Susskind Dirac operator
v1=reshape(U1(:,1),N,N);
v2=reshape(U1(:,2),N,N);
clear U1
U=sparse(N2,N2);
p=[N,1:N-1]; eta=sparse(diag((-1).^(0:N-1)));
for k=1:N; %loop over `time'
  e_k=sparse(N,1); e_k(k) = 1; I = sparse(diag(e_k)); T=sparse(I(:,p));
  % Upper triangular
  U=U+kron(diag(v1(k,:)),T);
  U=U+kron(T,diag(v2(:,k))*eta);
end
A1=mass*speye(N2)+(U-U');

% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
