function [col] = localcolinearity(R,evec)

col = zeros(size(evec,2),1);
for i = 1:size(evec,2)
    prod = R'*evec(:,i);
    prod = R*prod;
    col(i) = norm(evec(:,i)-prod);
end

scatter([1:size(evec,2)],col,'filled','r')