%outer solver parameters
finetol = 1e-8;
mMG = 8;
finecyclim = 10;

%coarsening parameters
N = 16; %size of lattice in one dim
l = 1;  % level of coarsening
gridsz = 4; %size of subgrids within the lattice
nsub = N/gridsz; %number of grids in one dimension
num = 8; %number of setup vectors (near null or near evec)
typevec = 1; %1 gives near null vecs, 2 gives near evecs
nntol = 1e-4; %tolerance to compute the near null vecs 

% critical mass for 128 latt = -0.0676
% critical mass for 256 latt = -0.0676
% critical mass for 64 latt  = -0.0770 this is not accurate
% use avg for critical mass: m = -0.0771
%mass = -0.0706; %mass for Dirac op
%mass = -0.0666;

deflag = 1; % 0 for no deflation on coarse level
            % 1 for deflating coarse level
rhsnum = 1;