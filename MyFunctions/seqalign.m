function [smat,tmat,alignval,alignstr] = seqalign(seq1, seq2, gappen, mmpen,algorithm) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function to find the scoring matrix between two sequences and
% their subsequent alignment
% 
% INPUTS
% 
% seq1, seq2 - the two sequences, of lengths m and n, that are to be compared
% 
% gappen - the penatly for choosing a gap over a diagonal move
%
% mmpen - the penalty for a mismatch
% 
% algorithm - choose the algorith preferred for sequence alignment. NW =
% Needleman-Wunsch, SW = Smith-Waterman.
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

if nargin < 1, error('Sequence alignment requires at least two input arguments'); end
if nargin < 2, error('Sequence alignment requires at least two input arguments'); end
if nargin < 3, gappen = 1; end
if nargin < 4, mmpen = 1; end
    
l1 = length(seq1)+1;
l2 = length(seq2)+1;
smat = zeros(l1,l2);
tmat = zeros(l1,l2);
tmatvals = [45 0 90];

%Needleman-Wunsch
%set the outer bounderies according to the gap penalty 

for i = 2:l1
    smat(i,1) = smat(i-1,1) - gappen;
end

for i = 2:l2
    smat(1,i) = smat(1,i-1) - gappen;
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
            smat(i,j) = score;
            tmat(i,j) = tmatvals(ind);
        end
    end
end


%calculate the alignment score and the alignment string
alignval = smat(i,j);
alignstr = repmat(char(0),3,l1+l2);
c=1;

while (i~=1 || j~=1) && c<(l1+l2+1) 
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