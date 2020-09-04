fineparams;
fine = zeros(2,1);
avgcoarse = zeros(2,1);
avgint = zeros(2,1);
measuret = zeros(2,1);
iters = zeros(2,1);

n = size(Af,1);

%initializations
b = ((sqrt(0.5)*(randn(n,1) + 1i*randn(n,1)))); %rand gauss dist rhs

tic();
[x,mvptotal,rnout,rnormvec,time,coarsemvps,intmvps,finemvps,iter,cdots] = fgmresMG(Af,Ahat,Ahh,P,Ph,b,mMG,finetol,finecyclim,deflag,rhsnum); %changed to Af'Af for testing 10/17/19
time2 = toc();
iters(1) = iter;
measuret(1) = time2;
fine(1) = finemvps;
avgcoarse(1) = coarsemvps/iter;
avgint(1) = intmvps/iter;

subplot(1,2,1)
hold on
set(gca,'yscale','log')
plot(mvptotal,rnormvec,'--b*')
xlabel('Fine Equivalent Mvps')
ylabel('Residual Norm')

deflag = 0;
tic();
[xd,mvptotal,rnout,rnormvec,time,coarsemvps,intmvps,finemvps,iter,cdots] = fgmresMG(Af,Ahat,Ahh,P,Ph,b,mMG,finetol,finecyclim,deflag,rhsnum); %changed to Af'Af for testing 10/17/19
time2 = toc();
iters(2) = iter;
measuret(2) = time2;
plot(mvptotal,rnormvec,'--r*')
ylim([10^-8 10^2])
legend('Deflated Coarse Solve','No Deflation')

fine(2) = finemvps;
avgcoarse(2) = coarsemvps/iter;
avgint(2) = intmvps/iter;

subplot(1,2,2)
hold on
X = categorical({'Fine Mvps','Avg Int Mvps','Avg Coarse Mvps','Time','Outer Iterations'});
Y = [fine,avgint,avgcoarse,measuret,iters];
bar(X,Y');
legend('Deflated','No Deflation')
legend('location','best')