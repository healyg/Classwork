function [coeff, score, eigs, ncoeff, nscore, nev] = pcnew(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My own personal PCA function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = size(x);
dev = zeros(s(1),s(2));

for i = 1:s(2)
    t = x(:,i);
    m = mean(t);
    for j = 1:s(1)
        dev(j,i) = x(j,i) - m;
    end
end
S = cov(dev); 

[coeff, lambdas] = eig(S);

eigs = zeros(length(lambdas),1);

for i = 1:length(lambdas)
    for j = 1:length(lambdas)
        if i == j
            eigs(i) = lambdas(i,j);
        end
    end
end

score = coeff*dev';
score = score';
ss = size(score);
sp = size(coeff);
ncoeff = zeros(sp(1),sp(2));
nev = sort(eigs,1,'descend');
nscore = zeros(ss(1),ss(2));

for i = 1:length(eigs) 
    ind = nev == eigs(i);
    ncoeff(:,ind) = coeff(:,i);
    nscore(:,ind) = score(:,i);
end
end
