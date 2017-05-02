function [smat,tmat,alignval,alignstr] = scormat(seq1, seq2, gappen) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function to find the scoring matrix between two sequences, in
% anticipation of running a Needleman-Wunsch Algorithm.
% 
% INPUTS
% 
% seq1, seq2 - the two sequences, of lengths m and n, that are to be compared
% 
% gappen - the penatly for choosing a gap over a diagonal move
% 
% OUTPUTS 
% 
% smat - an mxn scoring matrix
% 
% tmat - an mxn traceback matrix, where each cell indicates where the value
% came from (U = up, L = left, D = diagonal)
%
% alignval - the sequence alignment value 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


l1 = length(seq1)+1;
l2 = length(seq2)+1;
smat = zeros(l1,l2);
tmat = zeros(l1,l2);
tmatvals = [45 0 90];

%set the outer bounderies according to the gap penalty 

for i = 2:l1
    smat(i,1) = smat(i-1,1) - gappen;
end

for i = 2:l2
    smat(1,i) = smat(1,i-1) - gappen;
end

% create the scoring matrix, assuming seq is an array of chars 
if ischar(seq1) && ischar(seq2)
    for i = 2:l1
        for j = 2:l2
            if seq1(i-1) == seq2(j-1)
                smat(i,j) = smat(i-1,j-1) + 1;
                tmat(i,j) = 45;
            else
                diag = smat(i-1,j-1) - 1;
                up = smat(i-1,j) - 1;
                left = smat(i,j-1) - 1;
                [score,ind] = max([diag up left]);
                smat(i,j) = score;
                tmat(i,j) = tmatvals(ind);
            end
        end
    end
end

alignval = smat(i,j);
alignstr = repmat(char(0),2,l1+l2);
c=1;

while (i~=1 || j~=1) && c<(l1+l2+1) 
    if tmat(i,j)==45
        i = i-1;
        j = j-1;
        alignstr(1,c)=seq1(i);
        alignstr(2,c)=seq2(j);
    elseif tmat(i,j)==90
        j = j-1;
        alignstr(1,c)='_';
        alignstr(2,c)=seq2(j);
    else 
        i = i-1;
        alignstr(1,c)=seq1(i);
        alignstr(2,c)='_';
    end
    alignval = alignval + smat(i,j);
    c=c+1;
end
alignstr = fliplr(alignstr);
end