function [hidstates, obvs] = hide(states,obvstates,trans_mat,emprobs,initprob,l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function....
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initstate = randsample(states, 1, true, initprob);

if isequal(initstate, states(1))
    curstate = 1;
else
    curstate = 2; 
end

hidstates = blanks(l);
obvs = blanks(l);

for i = 1:l
    % find the observed states
    
    omatrix = emprobs(curstate,:);
    r = rand; 
    if r<omatrix(1)
        obvs(i) = obvstates(1);
    elseif r<omatrix(2)+omatrix(1)
        obvs(i) = obvstates(2);
    elseif r <omatrix(3)+omatrix(2)+omatrix(1) 
        obvs(i) = obvstates(3);
    elseif r <omatrix(4)+omatrix(3)+omatrix(2)+omatrix(1) 
        obvs(i) = obvstates(4);
    end
    
    % find the hidden states
    
    trans = trans_mat(:,curstate);
    r = rand;
    if r<trans(curstate)
        hidstates(i) = states(curstate);
    else
        hidstates(i) = states(curstate);
        if curstate == 1
            curstate = 2;
        else 
            curstate = 1;
        end 
    end 
end 




end