function [z0,mgmvps,coarsemvps,totalintmvps,totalfinemvps,coarsedots] = MGprecDeflated_intsolve(A,Ahat,Ahh,P,Ph,v,z0,outer,rhsnum)


smoothtol = 1e-2;
nmax = 2;
intparams;
coarseparams;
totalfinemvps = 0;
totalintmvps = 0;
coarsedots = 0;

%pre-smooth the error
[z,~,smoothmvps] = GMRES(A,v,z0,smoothtol,nmax);
res = v-A*z;
totalfinemvps = totalfinemvps + smoothmvps + 1;
% totalfinemvps = presmooth + 1;


%restrict residual to intermediate grid
rint = P'*res;


%smooth the error on intermediate grid
[err,~,presmooth] = GMRES(Ahat,rint,zeros(size(Ahat,1),1),smoothtol,nmax);
%rps = rint-Ahat*err;
%partial solve
[xs,mvptotal,~,~,~] = gmresrestarted(Ahat,rint,err,20,.2*norm(rint),10);
%compute residual
rsmooth = rint-Ahat*xs;
totalintmvps = totalintmvps + presmooth + mvptotal + 1;


%restrict to coarse grid
rc = Ph'*rsmooth;

%solve on coarse grid
%commented out for Dr. Morgan data 10/26/19
if ((outer == 1) && (rhsnum == 1))
   [xc,~,~,hk,vk,rnout,mvps,dots] = gmresdrEIG(Ahh,rc,mdefl,kdefl,zeros(size(Ahh,1),1),rtoldefl,cyclimdeflfirst);
   %save(['residual_' num2str(outer) '.mat'],'rnout');
   coarsemvps = mvps(length(mvps))
   save('hk.mat','hk');
   save('vk.mat','vk');
   coarsedots = coarsedots + dots;
else
   load('hk.mat');
   load('vk.mat');
   [xc,gdr,projmvps,dots] = gmresproj(Ahh,rc,mdefl,kdefl,vk,hk,projtol,cyclimdefl);
%    if ( outer <= 10)
%    save(['residual_' num2str(outer) '.mat'],'gdr');
%    end
   coarsemvps = projmvps;
   coarsedots = coarsedots + dots;
end
%end comment out for Dr. Morgan data 10/26/19

% [xc,~,~,~,~,~,mvps] = gmresdrEIG(Ahh,rc,mdefl,kdefl,zeros(size(Ahh,1),1),rtoldefl,cyclimdefl); % put in for Dr. Morgan Data
% coarsemvps = mvps(length(mvps));

%prolong to intermediate grid
errc = Ph*xc;

%add the error to the intermediate solve
xint = xs + errc;

%smooth on intermediate grid
[xsmooth,~,postsmooth] = GMRES(Ahat,rint,xint,smoothtol,nmax);
totalintmvps = totalintmvps + postsmooth;

%prolongate the error back to fine level
err = P*xsmooth;

%add the error to solution on fine grid
z0 = z + err;

[z,~,postsmooth] = GMRES(A,v,z0,smoothtol,nmax);
totalfinemvps = totalfinemvps + postsmooth;


%post smooth on fine grid, added 10/17/19
% res = v-A*z0;
% [z,~,~,postsmooth] = GMRES(A,v,z0,smoothtol,nmax);
 
 z0 = z;
% totalfinemvps = totalfinemvps + postsmooth + 1;
% z0 = z0 + z;
% end add


mgmvps = totalfinemvps + (size(Ahat,1)/size(A,1))*totalintmvps + (size(Ahh,1)/size(A,1))*coarsemvps;