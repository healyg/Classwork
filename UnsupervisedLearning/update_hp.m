function [state]=update_hp(inp,W,state,prob,thr)

% updates the Hopfield function
s = size(state);
for i = 1:s(1)
    for j = 1:s(2)
        r = rand; 
        if r<=prob
            if inp(i,j) >=thr
                state(i,j) = 1;
            else 
                state(i,j) = -1;
            end
        else
        end
    end
end

end





