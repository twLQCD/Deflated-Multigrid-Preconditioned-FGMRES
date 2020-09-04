function [cdparray,cdmarray] = chiral_double(nullvec,N,...
                                             num,nv,ch)

% Inputs:

% nullvec   : the array containing "num" near null vectors
% N         : the dimension of the lattice 
% num       : the number of near null vectors
% nv        : factor to retain appropriate degrees of freedom
% ch        : the chiral doubling factor. 

% Outputs: 

% cdparray  : array containing "num" positively chirally
%             doubled near null vectors
% cdmarray  : array containing "num" negatively chirally
%             doubled near null vectors


% form chiral projectors
sig = [1 0; 0 -1];
pp = (1/2)*(eye(ch,ch) + sig);
pm = (1/2)*(eye(ch,ch) - sig);


k = 1;
for i = 1:num
    cdp = zeros(ch,nv,(N^2));
    cdm = zeros(ch,nv,(N^2));
    null = reshape(nullvec(:,i),ch,nv,(N^2));
    for j = 1:((N^2))
        cdp(:,:,j) = pp*null(:,:,j);
        cdm(:,:,j) = pm*null(:,:,j);
    end
    cdp = reshape(cdp,ch*nv*N^2,1);
    cdm = reshape(cdm,ch*nv*N^2,1);
    cdparray(:,k) = cdp;
    cdmarray(:,k) = cdm;
    k = k + 1;
end