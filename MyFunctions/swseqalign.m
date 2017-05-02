function [smat,tmat,alignval,alignstr] = swseqalign(seq1, seq2, gappen, mmpen) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function to find the scoring matrix between two sequences and
% their subsequent alignment using the Smith-Waterman Method
% 
% INPUTS
% 
% seq1, seq2 - the two sequences, of lengths m and n, that are to be compared
% 
% gappen - the penatly for choosing a gap over a diagonal move
%
% mmpen - the penalty for a mismatch
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
% alignstr - the actual best alignment 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l1 = length(seq1)+1;
l2 = length(seq2)+1;
smat = zeros(l1,l2);
tmat = zeros(l1,l2);
tmatvals = [45 0 90];

for i = 2:l2
    tmat(1,i) = 90;
end

% create the scoring matrix
for i = 2:l1
    for j = 2:l2
        if seq1(i-1) == seq2(j-1)
            smat(i,j) = smat(i-1,j-1) + 1;
            tmat(i,j) = 45;
        else
            diag = smat(i-1,j-1) - mmpen;
            up = smat(i-1,j) - gappen;
            left = smat(i,j-1) - gappen;
            [score,ind] = max([diag up left]);
            if score>0
                smat(i,j) = score;
            else 
                smat(i,j) = 0;
            end
            tmat(i,j) = tmatvals(ind);
        end
    end
end

%calculate the alignment score and the alignment string
[m,i] = max(smat);
[m2,j] = max(m);
i=i(j);
alignval = m2;
alignstr = repmat(char(0),3,l1+l2);
c=1;

while (i~=1 || j~=1) && c<(l1+l2+1) && smat(i,j)~=0
    if tmat(i,j)==45
        i = i-1;
        j = j-1;
        alignstr(1,c)=seq1(i);
        alignstr(3,c)=seq2(j);
        if seq1(i)==seq2(j)
            alignstr(2,c)='|';
        end
    elseif tmat(i,j)==90
        j = j-1;
        alignstr(1,c)='-';
        alignstr(3,c)=seq2(j);
    else 
        i = i-1;
        alignstr(3,c)='-';
        alignstr(1,c)=seq1(i);
    end
    c=c+1;
end
alignstr = fliplr(alignstr);
end