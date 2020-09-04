function [grids] = getsubgrids(N,nsub,gridsz,numgrids)
% Input:

% N          : dimension of lattice
% nsub       : number of grids in one dim
% gridsz     : is the size of grid
% numgrids   : the total number of grids

% Output:

% grids      : array of vectors containing
%              partitioned lattice points
                                                  
grids = sparse(gridsz^2,numgrids);
vec = [1:N^2];
vec = reshape(vec, N, N);
forcell = gridsz*ones(1,nsub);

j = N;
for i = 1:N
    newvec(i,:) = vec(:,j);
    j = j - 1;
end

B = mat2cell(newvec, forcell, forcell);

k = 1;
for i = nsub:-1:1
    for j = 1:nsub
        sgrid = cell2mat(B(i,j));
        sgrid = sgrid';
        ix = [gridsz:-1:1];
        sgrid = sgrid(:,ix);
        sgrid = reshape(sgrid,gridsz^2,1);
        grids(:,k) = sgrid;
        k = k + 1;
    end
end


        