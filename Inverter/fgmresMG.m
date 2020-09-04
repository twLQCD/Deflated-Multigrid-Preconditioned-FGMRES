function [x,mvptotal,rnout,rnormvec,time,totalcoarsemvps,totalintmvps,finemvps,outeriter,cdots] = fgmresMG(A,Ah,Ahh,P,Ph,b,m,rtol,cyclim,deflag,rhsnum)

n = size(A,1);
cycle = 1; j = 1; outeriter = 1;
mvptotal = [];
mvp = 0;
finemvps = 0;
totalintmvps = 0;
totalcoarsemvps = 0;
cdots = 0;

rncycle = [];
rnout = [];
rnormvec = [];
rnp = 0;
vsrv = b;

Z = sparse(n,m+1);

x = zeros(size(A,1),1);
xiter = zeros(size(A,1),1);
r = b;
norminit = norm(r);
rn = norminit;
c(1,1) = norm(r);
c(2:m+1,1) = zeros(m,1);
rnormiter = norm(r);


tic();
while( (cycle <= cyclim) && ( rnormiter > rtol ))
    
    v(:,1) = r/rn;
        
    while ( (j <= m) && ( rnormiter > rtol))
        
        %multigrid prec
        % in A, Ah, Ahh, P, Ph, v(:,j)
        % out: z, mvps
        if (deflag == 0)
%           [z,mgmvps,coarsemvps,intmvps,smoothmvps] = MGprecV2(A,Ah,Ahh,P,Ph,v(:,j),x,outeriter);
            [z,mgmvps,coarsemvps,intmvps,smoothmvps] = MGprecV2_intsolve(A,Ah,Ahh,P,Ph,v(:,j),x,outeriter);
            cdots = [];
            dots = [];
        elseif (deflag == 1)
          %[z,mgmvps,coarsemvps,intmvps,smoothmvps] = MGprecDeflatedV2(A,Ah,Ahh,P,Ph,v(:,j),x,outeriter,rhsnum);
          %%commented out for stagg test -TW 11/29/19
          [z,mgmvps,coarsemvps,intmvps,smoothmvps,dots] = MGprecDeflated_intsolve(A,Ah,Ahh,P,Ph,v(:,j),x,outeriter,rhsnum);
        end
        Z(:,j) = z;
        
        %w = A*v(:,j);
        w = A*z;
        %mvp = mvp + 1;
        totalintmvps = totalintmvps + intmvps;
        totalcoarsemvps = totalcoarsemvps + coarsemvps;
        finemvps = finemvps + smoothmvps + 1;
        mvp = mvp + mgmvps + 1;
        mvptotal = [mvptotal,mvp];
        cdots = cdots + dots;

        
        for i = 1:j
            h(i,j) = v(:,i)'*w;
            w = w - h(i,j)*v(:,i);
        end
        
        h(j+1,j) = norm(w);
        v(:,j+1) = w/h(j+1,j);

        %for residual at every iteration
        d(1:j,1) = (h(1:j+1,1:j)) \ c(1:j+1) ;
        srv(1:j+1,1) = c(1:j+1)-(h(1:j+1,1:j))*d(1:j,1);
        rnormiter = norm(srv(1:j+1,1));

        rnormvec = [rnormvec,rnormiter];
        
        if (rnormiter < rtol)
            x = x + Z(:,1:j)*d(1:j);
            fprintf("\n Linear equations have converged in fgmres \n")
        end
        
        outeriter = outeriter + 1;
        
        j = j + 1;
        

    end
    
  if ( rnormiter > rtol)
  d(1:m,1) = h(1:m+1,1:m) \ c(1:m+1,1) ;
  srv(1:m+1,1) = c(1:m+1)-(h(1:m+1,1:m)*d(1:m,1));


  %x = x + v(:,1:m)*d(1:m);
  x = x + Z(:,1:m)*d(1:m);
  end
  r = b-A*x; 

  rnouter = norm(r)
  rn = norm(r);
  
  
  c(1,1) = norm(r);
  c(2:m+1,1) = zeros(m,1);
  rnout = [rnout,rn];
    
    
    j = 1;
    cycle = cycle + 1;
end
time = toc();