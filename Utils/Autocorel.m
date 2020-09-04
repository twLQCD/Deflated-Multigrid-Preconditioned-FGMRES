function y=Autocorel(x,t);
% Computes autocorrelation function of time-shifted data
x=x(:);
N=max(size(x));
X=[];
for i=1:t;
  y=x(i:N-t+i);
  X=[X,y];
end
C=cov(X);
y=C(:,1)/C(1,1);
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
