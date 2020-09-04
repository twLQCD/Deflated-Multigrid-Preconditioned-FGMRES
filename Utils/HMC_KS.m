function [A2,Plaq, Q_top, Wloop, theta2, stat] = HMC_KS(theta1,iconf)
% HMC for 4 flavours of KS fermions
beta = 6; %Coupling constant
mass = 0.35; %bare mass
warmuplim = 0;
warmup = 1;
N = 384; N2 = N^2;
% Make nearest neighbours list
for j2=0:N-1
for j1=0:N-1
  k       = 1 +            j1 +            j2*N;
  kp(k,1) = 1 + mod(j1+1  ,N) +            j2*N;
  kp(k,2) = 1 +            j1 + mod(j2+1  ,N)*N;
  km(k,1) = 1 + mod(j1-1+N,N) +            j2*N;
  km(k,2) = 1 +            j1 + mod(j2-1+N,N)*N;
end
end
% Starting configuration theta1
if iconf==0
  theta1=2*pi*rand(N2,2);
end
% Start Monte Carlo
ntest = 0; Plaq = []; Wloop = []; Q_top = []; stat=[];
NMC = 100; NMD=10; deltat = 0.05;
for mc = 1:NMC
  fprintf("\n in Monte Carlo")
  tic()
  U1=cos(theta1)+sqrt(-1)*sin(theta1); % Current configuration U1
  % get fermion matrix
  A1=Dirac_KS(mass,N,N2,U1);
  eta0=rand(N2,1)+sqrt(-1)*rand(N2,1); % Gaussian noise
  phi=A1*eta0; % Update pseudofermion field
  P=randn(N2,2); % Refresh momenta
% Compute H1
  theta=theta1(:,1)+theta1(kp(:,1),2)-theta1(:,2)-theta1(kp(:,2),1);
  H1 = norm(P)^2/2 - beta*sum(cos(theta)) + norm(eta0)^2;
% Propose theta2 using MD evolution
  %chi=A1'\eta0;
  [chi,~,~,~,~] = gmresrestarted(A1',eta0,100,1e-8,250);
  theta2=theta1;
  % Advance momenta half step
  [pf,pg]=Force_KS(N2,kp,km,theta1,eta0,chi,U1);
  pdot=-beta*pg+pf;
  P=P+pdot*deltat/2;
  % MD loop
  for md=1:NMD-1
    theta2=theta2+P*deltat; % Advance angles full step
    U2=cos(theta2)+sqrt(-1)*sin(theta2); % Update U2
    % get fermion matrix
    A2=Dirac_KS(mass,N,N2,U2);
    %eta=A2\phi; chi=A2'\eta; % Compute eta, chi
    [eta,~,~,~,~] = gmresrestarted(A2,phi,100,1e-8,250);
    [chi,~,~,~,~] = gmresrestarted(A2',eta,100,1e-8,250);
    % Advance momenta full step
    [pf,pg]=Force_KS(N2,kp,km,theta2,eta,chi,U2);
    pdot=-beta*pg+pf;
    P=P+pdot*deltat;
    md
  end
  theta2=theta2+P*deltat; % Advance angles full step
  U2=cos(theta2)+sqrt(-1)*sin(theta2); % New configuration U2
  % get fermion matrix
  A2=Dirac_KS(mass,N,N2,U2);
  %eta=A2\phi; chi=A2'\eta; % Compute eta, chi
  [eta,~,~,~,~] = gmresrestarted(A2,phi,100,1e-8,250);
  [chi,~,~,~,~] = gmresrestarted(A2',eta,100,1e-8,250);
  % Advance momenta half step
  [pf,pg]=Force_KS(N2,kp,km,theta2,eta,chi,U2);
  pdot=-beta*pg+pf;
  P=P+pdot*deltat/2;
% Compute H2
  theta=theta2(:,1)+theta2(kp(:,1),2)-theta2(:,2)-theta2(kp(:,2),1);
  [x,~,~,~,~] = gmresrestarted(A2,eta0,100,1e-8,250);
  %H2 = norm(P)^2/2 - beta*sum(cos(theta)) + norm(A2\eta0)^2;
  H2 = norm(P)^2/2 - beta*sum(cos(theta)) + norm(x)^2;
% Metropolis test
  R = min([1,exp(-(H2-H1))]);
  random=rand;
  istat = [mc, random, R, H2-H1];
  stat = [stat; istat];
  test=0;
  if (warmup <= warmuplim)
      theta1 = theta2;
      warmup = warmup + 1;
  else
  if random<R, test=1; end
  if test==1
    theta1=theta2;
    ntest = ntest + 1;
    % Compute Plaquette & Topological charge
    theta=theta1(:,1)+theta1(kp(:,1),2)-theta1(:,2)-theta1(kp(:,2),1);
    plaq=mean(cos(theta)); Plaq = [Plaq;plaq];
    q_top=sum(sin(theta))/(2*pi); Q_top = [Q_top;q_top];
    [mc, ntest/mc, plaq, q_top] % Display some results
    % Compute Wilson loops
    wlp=wloop(N,N2,kp,km,theta1);
    Wloop = [Wloop;wlp'];
  end
  end
  mc
  toc()
end
% Copyright (C) 2006 Artan Borici.
% This program is a free software licensed under the terms of the GNU General Public License
