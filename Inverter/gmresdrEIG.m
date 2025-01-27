function [x,th,evector,hk,vk,rnout,mvps,dots] = gmresdrEIG(A,b,m,k,x,rtol,cyclim)
% gmresdr   solves systems of linear equations using a deflated restarted 
%           method. 
%
% Synopsis: 
%     [x,vk,hk,gdr] = gmresdr(A,b)
%     [x,vk,hk,gdr] = gmresdr(A,b,m,k)
%     [x,vk,hk,gdr] = gmresdr(A,b,m,k,rtol,cyclim)
%
% Input:
%   A      = an nxn non-hermitian matrix
%   b      = single right hand side
%   m      = size of Krylov Subspace (number of cycles before restarting)
%            Default = 10
%   k      = number of deflated eigenvalues/eigenvectors
%            Default = 5
%   rtol   = relative residual nececcary
%            Default = 1e-6
%   cyclim = maximum number of cycles to perform
%            Default = 100
%
% Output:
%   x   = solution vector to Ax=b
%   vk  = k approximate eigenvectors
%   hk  = k approcimate eigenvalues (located on diagonal)
%   gdr = residual at each iteration
%   
% Notes: In order to use gmresdr with gmresproj the relative residual
% should be run excessivly high to gain more accurate eigenvalue/vector
% information. This will lead to much faster convergence of gmresproj.

% if (nargin < 7)
%   fid=1;
%   if (nargin < 5)
%     rtol = 1e-6;
%     cyclim = 10;
%     if (nargin < 3)
%       m = 10;
%       k = 5;
%     end
%   end
% end

n = size(A,1);

%x = zeros(n,1);
cycle = 1;
j = 1;
rnout = [];
norm1 = [];
mvp = 0;
mvps = [];
m0 = m;
normtmp = 1;
dots = 0;

rninit = norm(b(:,1));
rn = rninit;
linnorm = rn;
r = b;
vn = norm(r);
v = r/vn;
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);


while ( (rn > rtol) && (cycle <= cyclim) )
%while ( (cycle <= cyclim) )
  while ( j <= m ) 
    wv = A*v(:,j);
   %wv = A*gamma5(v(:,j),n,1);
   %wv = A*epsmult(v(:,j),n,16);
    mvp = mvp + 1;
    mvps = [mvps,mvp];
    
    %FOR NORMAL EQNS
%      wv = A'*wv;
%      mvps = mvps + 1;
    %END NORMAL EQNS
    
    vnf = norm(wv);

    for i = 1:j
      h(i,j) = v(:,i)'*wv;
       wv = wv - h(i,j) * v(:,i);
       dots = dots + 1;
    end
    vn = norm(wv);

    %--------reorthogonalization section-------%
    if( vn < 1.1*vnf )
      for i = 1:j
        dot = v(:,i)'*wv;
        dots = dots + 1;
        wv = wv - dot * v(:,i);
        h(i,j) = h(i,j) + dot;
      end
      vn = norm(wv);
    end
    %--------------------------------------------------%

    h(j+1,j) = vn;
    v(:,j+1) = wv/h(j+1,j);

  % output data per cycle
  % disp(['iteration: ',num2str(j)]);

    j = j + 1;
  end
  
  
  j = 1;
  d(1:m,1) = h(1:m+1,1:m) \ c(1:m+1) ;
  srv(1:m+1,1) = c(1:m+1)-(h(1:m+1,1:m)*d(1:m,1));

  %----Set up and solve linear equations.-----%
  x(:,1) = x(:,1) + v(:,1:m)*d(1:m);
  r = v(:,1:m+1)*srv; %commented out for testing TW 3/23/19
  rn = norm(r);
  if (rn < rtol)
      fprintf("\n Linear equations have converged in gmres-dr \n")
  end
  linnorm = norm(A*x-b);
  %linnorm = norm(A*gamma5(x,n,1)-b)
  
  %for Normal eqns
  %linnorm = norm(A'*A*x-b)
  
  rnout = [rnout,rn];
  
  j = j+1;

  hh = h(1:m,1:m);
  em = zeros(m,1);
  em(m) = 1;
  ff = hh' \ em;

  hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
  hh=full(hh); [g,dd] = eig(hh,'nobalance'); 
% no balance needed!!!-WW; 03/25/2014!

  dd = diag(dd);
  dabs = abs(dd);
  [dabs,ind] = sort(dabs);
  gout = g(:,ind);
  ddout = dd(ind);
  th = dd(ind);
  gg = g(:,ind(1:k));

  for i=1:k
    rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i); 
    tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
    tvn = norm(tv);
    rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
    tha(cycle,i) = th(i);
    rhoa(cycle,i) = rho(i);
  end

  [rnasort,ind] = sort(rna(cycle,:));
  rna(cycle,1:k)=rna(cycle,ind(1:k));
  gg = gg(:,ind(1:k));
  th = th(ind);
  

  greal = gg;

  greal(m+1,1:k) = zeros(1,k);

% Chris 

  beta = h(m+1,m);
  punty = zeros(m,1);
  punty = -beta*ff;
  greal(1:m,k+1) = punty(1:m,1);
  greal(m+1,k+1) = 1.0;

%  greal(:,k+1) = srv;

% end Chris


  [gon,rr] = qr(greal(:,1:k+1),0);
  hnew = gon'*h*gon(1:m,1:k);
  h(k+1,:) = zeros(1,m);

  j = 1;
  rtolev = 1e-11;
  while ( j<= k && rna(cycle,j) <= rtolev)
    hnew(j+1:k+1,j) = zeros(k-j+1,1);
    j = j + 1;
  end



  evector = v(:,1:m) * greal(1:m,1:k);
  normtmp = norm(A*evector(:,1)-th(1)*evector(:,1));
  normtmp = norm(A*evector(:,k)-th(k)*evector(:,k));
%   if ( normtmp <= 1e-3)
%       cycle = cyclim;
%   end


  h(1:k+1,1:k) = hnew;


  c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
  c(k+2:m+1,1) = zeros(m-k,1);


  work = v*gon;
  v(:,1:k+1) = work;

  %section for just reorthog. one vector, v_{k+1}
  for i = 1:k
    dot = v(:,i)'*v(:,k+1) ;
    dots = dots + 1;
    v(:,k+1) = v(:,k+1) - dot * v(:,i);
  end

  v(:,k+1) = v(:,k+1)/norm(v(:,k+1));

  % output data per cycle
  %fprintf(fid,'Cycle: %d  Rel Res Norm: %12.8g\n',cycle,rn/rninit);
  vout = v(:,1:m);
  j = k + 1;
  cycle = cycle + 1;
  
end

hk = h(1:k+1,1:k);
vk=v(:,1:k+1);

%semilogy(rna)                                                               
% 
% if (rn/rninit < rtol) && (cycle-1 <= cyclim)
%   fprintf(fid,'gmresdrEIG(%d,%d) converged in %d cycles with a relative residual of %12.8g\n', ...
%           m,k,cycle-1,rn/rninit);
% else
%   fprintf(fid,'gmresdrEIG(%d,%d) stoped after %d cycles without converging to the desired tolerence of %12.8g\n', ...
%           m,k,cycle-1,rn/rninit);
%   fprintf(fid,'a relative residual of %12.8g has been reached\n',rn/rninit);
% end

return
