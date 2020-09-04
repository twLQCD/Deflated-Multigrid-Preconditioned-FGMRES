function [z0,mgmvps,coarsemvps,totalintmvps,totalfinemvps] = MGprecV2_intsolve(A,Ahat,Ahh,P,Ph,v,z0,outer)

smoothtol = 1e-2;
nmax = 2;
intparams;
coarseparams;
totalfinemvps = 0;
totalintmvps = 0;

%pre-smooth the error
[z,~,smoothmvps] = GMRES(A,v,z0,smoothtol,nmax); 
% switched to minres
%[z,~,~,smoothmvps,~] = minres(A,v,smoothtol,nmax,[],[],z0);
res = v-A*z;
totalfinemvps = totalfinemvps + smoothmvps + 1;
% totalfinemvps = presmooth + 1;


%restrict residual to intermediate grid
rint = P'*res;


%smooth the error on intermediate grid
%comment out below for minres test
[err,~,presmooth] = GMRES(Ahat,rint,zeros(size(Ahat,1),1),smoothtol,nmax); 
[xs,mvptotal,~,~,~] = gmresrestarted(Ahat,rint,err,20,.2*norm(rint),10);
%compute residual
rsmooth = rint-Ahat*xs;

totalintmvps = totalintmvps + presmooth + mvptotal + 1;


%restrict to coarse grid
rc = Ph'*rsmooth;

%solve on coarse grid
if (outer == 1)
[xc,cflag,relres,coarsemvps,resvec] = pcg(Ahh'*Ahh,Ahh'*rc,coarsetol,coarsecyclim);
cflag;
% save(['residual_cg' num2str(outer) '.mat'],'resvec');
else 
  [xc,cflag,relres,coarsemvps,resvec] = pcg(Ahh'*Ahh,Ahh'*rc,matchdefl,coarsecyclim);
  cflag;
%   if ( outer <= 10)
%   save(['residual_cg' num2str(outer) '.mat'],'resvec');
%   end
end
coarsemvps = 2*coarsemvps;

%prolong to intermediate grid
errc = Ph*xc;

%add the error to the intermediate solve
xint = xs + errc;

%smooth on intermediate grid
% comment out below for minres test
[xsmooth,~,postsmooth] = GMRES(Ahat,rint,xint,smoothtol,nmax);
%[xsmooth,~,~,postsmooth,~] = minres(Ahat,rint,smoothtol,nmax,[],[],xint);

totalintmvps = totalintmvps + postsmooth;

%prolongate the error back to fine level
err = P*xsmooth;

%add the error to solution on fine grid
z0 = z + err;

%comment out below for minres test
[z,~,postsmooth] = GMRES(A,v,z0,smoothtol,nmax);
%[z,~,~,postsmooth,~] = minres(A,v,smoothtol,nmax,[],[],z0);
totalfinemvps = totalfinemvps + postsmooth;


%post smooth on fine grid, added 10/17/19
% res = v-A*z0;
% [z,~,~,postsmooth] = GMRES(A,v,z0,smoothtol,nmax);
 
 z0 = z;
% totalfinemvps = totalfinemvps + postsmooth + 1;
% z0 = z0 + z;
% end add


mgmvps = totalfinemvps + (size(Ahat,1)/size(A,1))*totalintmvps + (size(Ahh,1)/size(A,1))*coarsemvps;