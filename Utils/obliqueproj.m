function [obproj] = obliqueproj(P,A,Ah,evec)

obproj = zeros(size(evec,2),1);

for i = 1:size(evec,2)
    prod = P'*A*evec(:,i);
    prod = Ah \ prod;
    prod = P*prod;
    prod = evec(:,i)- prod;
    obproj(i) = norm(prod);
end

scatter([1:size(evec,2)],obproj,'filled','m')
    