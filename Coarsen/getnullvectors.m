function [nullvec,totalmvps] = getnullvectors(A,num,tol)

% Input:

% A          : the l-level Wilson-Dirac operator
% num        : the number of near null vectors desired
% tol        : the tolerance for the near null vector
%              calculation


% Output:

% nullvec    : the array containing num near null vectors
% totalmvps  : the total number of matrix vector products for
%              the generation of near null vectors 

n = size(A,1);
nmax = 250;
z = sparse(n,1);
nullvec = zeros(n,num);
totalmvps = 0;

% generate near null vectors
for i = 1:num
    z0 = sqrt(0.5)*(randn(n,1)+1i*randn(n,1));
    b = -A*A'*z0;
    [z,~,~,iter] = pcg(A*A',b,tol,nmax,[],[]);
    nullvec(:,i) = z0 + z;
    totalmvps = 2*iter + totalmvps;
end

% othonormalize with modified gram schmidt
for i = 1:num
    for j = 1:i-1
    dot = nullvec(:,j)'*nullvec(:,i);
    nullvec(:,i) = nullvec(:,i) - dot * nullvec(:,j);
    end
  nullvec(:,i) = nullvec(:,i)/norm(nullvec(:,i));
end



