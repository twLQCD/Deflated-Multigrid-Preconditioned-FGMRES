function err=Binning(x,t);
% Compute errors by blocking data
x=x(:);
N=max(size(x));
err=std(x)/sqrt(N);
for i=2:t;
  n=i*floor(N/i);
  y=mean(reshape(x(1:n),i,n/i));
  err=[err;std(y)/sqrt(n/i)];
end
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
