function [Wout,ind] = train_cl(W,inp,lr)
%competitive learning rule

nfeatures = size(W,2);
ncat = size(W,1);
    
out=W*inp;
[mx,ind]=max(out);
W(ind,:)=W(ind,:)+lr*(inp'-W(ind,:));
Wout=W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));

while all(ind == ind(1))
    W=rand(ncat,nfeatures);
    Wout=W./repmat(sqrt(sum(W.^2,2)),1,nfeatures);
    out=Wout*inp;
    [mx,ind]=max(out);
end

end