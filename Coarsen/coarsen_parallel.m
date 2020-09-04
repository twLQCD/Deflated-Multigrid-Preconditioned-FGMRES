function [Ahat,P,finemvps] = coarsen_parallel(A,N,l,nsub,...
                             gridsz,numi,num,tol,numprocs)

% Inputs:

% A           : l-level Wilson-Dirac operator
% N           : lattice dim
% l           : the level of coarsening
% nsub        : number of subgrids in one dim. Since lattice
%               is isotropic, nsub = N/gridsz
% gridsz      : the size of the subgrids in one dim
% num         : the number of desired near nullvectors for
%               this coarsening
% numi        : the number of near null vectors for the
%               previous coarsening if l > 1
% tol         : the tolerance for the null vectors
% numprocs    : the number of processors desired for
%               assembling the prolongator, limited 
%               to number of cores per node
           
% Outputs:

% Ahat       : l+1-level coarse operator
% P          : prolongator matrix from l-1 to l grid
% finemvps   : matrix vector product counter

% below retains appropriate degrees of freedom for
% l-level coarsening
if ( l == 1 )
    nv = 1;
elseif ( l > 1)
    nv = numi;
end

% preallocations
ch = 2; %chirality factor, corresponds to spin d.o.f
cdparray = zeros(ch*nv*N^2,num);
cdmarray = zeros(ch*nv*N^2,num);
numgrids = N^2/gridsz^2;
rowrest = num*numgrids*ch;
R = sparse(size(A,1),rowrest);

[nullvec,finemvps] = getnullvectors(A,num,tol);

[cdparray,cdmarray] = chiral_double(nullvec,N,num,nv,ch);

[grids] = getsubgrids(N,nsub,gridsz,numgrids);

parpool('local',numprocs)

spmd(numprocs)
    
    ij = (numgrids*num*ch)/numprocs;
    jk = numgrids/numprocs;

    k = (labindex-1)*ij + 1;

tic()
for j = (labindex-1)*jk+1:(labindex*jk)
    
    for i = 1:num
    
    nullp = reshape(cdparray(:,i),ch*nv*N^2,1);
    nullp = reshape(nullp,ch,nv,N^2);
    nullm = reshape(cdmarray(:,i),ch*nv*N^2,1);
    nullm = reshape(nullm,ch,nv,N^2);
    
        
    rp = zeros(ch,nv,N^2);
    rm = zeros(ch,nv,N^2);
    rp(:,:,grids(:,j)) = nullp(:,:,grids(:,j));
    rm(:,:,grids(:,j)) = nullm(:,:,grids(:,j));
    rp = reshape(rp,ch*nv*N^2,1);
    rm = reshape(rm,ch*nv*N^2,1);
                
    R(:,k) = rp;
    R(:,k+1) = rm;
        
    if ( mod(k+1,2*num) == 0 )
      for ii = ((k+1)-2*num)+1:(k+1)
         for jj = ((k+1)-2*num)+1:(ii-1)
             dot = R(:,jj)'*R(:,ii);
             R(:,ii) = R(:,ii) - dot*R(:,jj);
         end
             R(:,ii) = R(:,ii)/norm(R(:,ii));
      end
    end    

       k = k + 2;
    end %i
end %j

end %spmd
Rout = R(:);

P = sparse(size(A,1),rowrest);
for i = 1:numprocs
    P = P + cell2mat(Rout(i));
end

delete(gcp('nocreate'))
Ahat = P'*A*P;

clear null cdp cdm nullvec rp rm nullp nullm grids Rout R
