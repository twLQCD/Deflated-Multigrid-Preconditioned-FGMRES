function [a,b]=cdot5(u,v)
% gamma5 inner product u'*gamma5*v
% here for Schwinger model: gamma5->sigma3
u=u(:);
v=v(:);
N=max(size(u));
u=reshape(u,2,N/2)';
v=reshape(v,2,N/2)';
a1=u(:,1)'*v(:,1);
a2=u(:,2)'*v(:,2);
a=conj(a1-a2);
b=conj(a1+a2);
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
