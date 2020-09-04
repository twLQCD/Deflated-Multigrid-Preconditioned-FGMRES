function [x,mvptotal,mvps,rnout,rnormvec,dots] = gmresrestarted(A,b,x,m,rtol,cyclim)

n = size(A,1);
cycle = 1; j = 1;
m1 = m;
mvptotal = 0;
mvps = [];
dots = 0;
rncycle = [];
rnout = [];
rnp = 0;
vsrv = b;
rnormvec = [];
evalnorm = [];
rnormiter = norm(b);

%x = zeros(size(A,1),1);
r = b;
norminit = norm(r);
rn = norminit;
rndiffvec = zeros(2,1);
c(1,1) = norm(r);
c(2:m+1,1) = zeros(m,1);


while( (cycle <= cyclim) && ( rnormiter > rtol ))
    
    v(:,1) = r/rn;
  
    
    mvp = 0;
    
    while ( j <= m)
        
        w = A*v(:,j);
        
        mvp = mvp + 1;
        mvptotal = mvptotal + 1;
        mvps = [mvps,mvptotal];

        
        for i = 1:j
            h(i,j) = v(:,i)'*w;
            dots = dots + 1;
            w = w - h(i,j)*v(:,i);
        end
        
        h(j+1,j) = norm(w);
        v(:,j+1) = w/h(j+1,j);

        %for residual at every iteration
        d(1:j,1) = (h(1:j+1,1:j)) \ c(1:j+1) ;
        srv(1:j+1,1) = c(1:j+1)-(h(1:j+1,1:j))*d(1:j,1);
        rnormiter = norm(srv(1:j+1,1));
        rnormvec = [rnormvec,rnormiter];
        
        j = j + 1;
        
    end
    

  d(1:m,1) = h(1:m+1,1:m) \ c(1:m+1,1) ;
  srv(1:m+1,1) = c(1:m+1)-(h(1:m+1,1:m)*d(1:m,1));

  x = x + v(:,1:m)*d(1:m);
  r = b-A*x; 

  rn = norm(r);
 
%   if ( cycle > 1)
%   if ( (evalnorm(cycle)/evalnorm(cycle-1)) >= 0.90)
%       A = Af - rho*speye(n,n) + b*b';
%       x = zeros(n,1);
%       r = b/norm(b);
%   end
%   end
  
  vsrv = v(:,1:m+1)*srv(1:m+1,1);
  nvsrv = norm(vsrv);
  
  c(1,1) = norm(r);
  c(2:m+1,1) = zeros(m,1);
  rnout = [rnout,rn];
    
    
    j = 1;
    cycle = cycle + 1;
end


    
