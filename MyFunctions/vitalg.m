function [mlhiddenstates,prob] = vitalg(states,obvstates,trans_mat,emprobs,initprob,truobvs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the Viterbi algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = length(truobvs);

mlhiddenstates = blanks(l);

nprob = initprob(1);
cprob = initprob(2);
for j = 1:length(obvstates)
    if isequal(truobvs(1),obvstates(j))
        s = j;
    end
end
[prob,curstate] = max([(log(nprob) + log(emprobs(1,s))),(log(cprob) + log(emprobs(2,s)))]);

mlhiddenstates(1)=states(curstate);

for i = 2:l
    if curstate == 1
        op = 2;
    else 
        op = 1;
    end 
    
    s = 0;
    for j = 1:length(obvstates)
        if isequal(truobvs(i),obvstates(j))
            s = j;
        end
    end
    [prob,new] = max([(prob + log(emprobs(curstate,s))+log(trans_mat(curstate,curstate))), ...
        (prob + log(emprobs(op,s))+log(trans_mat(op,curstate)))]);
    if new == 2
        curstate = op;
    end
    mlhiddenstates(i) = states(curstate);
end



end