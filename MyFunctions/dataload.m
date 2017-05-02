function [it] = dataload(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = 'cell_11l10.mat';
s2 = 'cell_11l18.mat';
s3 = 'cell_12r08.mat';
s4 = 'cell_14l07.mat';
s5 = 'cell_14l13.mat';
s6 = 'cell_11l04.mat';
s7 = 'cell_11l07.mat';
s8 = 'cell_11l14.mat';
s9 = 'cell_11l15.mat';
s10 = 'cell_11l16.mat';
s11 = 'cell_11l19.mat';
s12 = 'cell_12r03.mat';
s13 = 'cell_12r05.mat';
s14 = 'cell_12r06.mat';
s15 = 'cell_12r11.mat';
s16 = 'cell_12r12.mat';
s17 = 'cell_12r17.mat';
s18 = 'cell_12r18.mat';
s19 = 'cell_12r20.mat';
s20 = 'cell_12r21.mat';
s21 = 'cell_12r22.mat';
s22 = 'cell_14l01.mat';
s23 = 'cell_14l02.mat';
s24 = 'cell_14l03.mat';
s25 = 'cell_14l04.mat';
s26 = 'cell_14l05.mat';
s27 = 'cell_14l06.mat';
s28 = 'cell_14l08.mat';
s29 = 'cell_14l09.mat';
s30 = 'cell_14l10.mat';
s31 = 'cell_14l11.mat';
s32 = 'cell_14l12.mat';
s33 = 'cell_14l15.mat';
s34 = 'cell_14l16.mat';
s35 = 'cell_14l19.mat';
s36 = 'cell_12r15.mat';

it = struct;

if n >= 1 
    cell1 = load('cell_11l10.mat');
    s = size(cell1.data);
    if s(2)>14 
        cell1.data(:,15:s(2),:) = [];
    end
    it.first = cell1.data;
end
if n >= 2
    cell2 = load('cell_11l18.mat');
    it.second = cell2.data;
end
if n >= 3
    cell3 = load('cell_12r08.mat');
    if s(2)>14 
        cell1.data(:,15:s(2),:) = [];
    end
    it.third = cell3.data;
end 
if n >= 4
    cell4 = load('cell_14l07.mat');
    it.fourth = cell4.data;
end 
if n >= 5
    cell5 = load('cell_14l13.mat');
    it.fifth = cell5.data;
end
if n>=6
    cell6 = load('cell_11l04.mat');
    it.sixth = cell6.data;
end
if n>=7
    cell7 = load('cell_11l07.mat');
    it.seventh = cell7.data;
end
if n>=8
    cell8 = load('cell_11l14.mat');
    it.eighth = cell8.data;
end
if n>=9
    cell9 = load('cell_11l15.mat');
    it.ninth = cell9.data;
end
if n>=10
    cell1 = load('cell_11l16.mat');
    it.tenth = cell1.data;
end
if n>=11
    cell1 = load('cell_11l19.mat');
    it.eleventh = cell1.data;
end
if n>=12
    cell1 = load('cell_12r03.mat');
    it.twelth = cell1.data;
end
if n >=13
    cell1 = load('cell_12r05.mat');
    it.thirteenth = cell1.data;
end 
if n>=14
    cell1 = load('cell_12r06.mat');
    it.fourteenth = cell1.data;
end 
if n>=15
    cell1 = load('cell_12r11.mat');
    it.fifteenth = cell1.data;
end
if n>=16
    cell1 = load('cell_12r12.mat');
    it.sixteenth = cell1.data;
end
if n>=17
    cell1 = load('cell_12r17.mat');
    it.seventeenth = cell1.data;
end
if n>=18
    cell1 = load('cell_12r18.mat');
    it.eighteenth = cell1.data;
end
if n>=19
    cell1 = load('cell_12r20.mat');
    it.nineteenth = cell1.data;
end
if n>=20
    cell1 = load('cell_12r21.mat');
    it.twentieth = cell1.data;
end
if n>=21
    cell1 = load('cell_12r22.mat');
    it.twentyfirst = cell1.data;
end
if n>=22
    cell1 = load('cell_14l01.mat');
    it.twentysecond = cell1.data;
end
if n>=23
    cell1 = load('cell_14l02.mat');
    it.twentythird = cell1.data;
end
if n>=24 
    cell1 = load('cell_14l03.mat');
    it.twentyfourth = cell1.data;
end
if n>=25
    cell1 = load('cell_14l04.mat');
    it.twentyfifth = cell1.data;
end
if n>=26
    cell1 = load('cell_14l05.mat');
    it.twentysixth = cell1.data;
end
if n>=27
    cell1 = load('cell_14l06.mat');
    it.twentyseventh = cell1.data;
end
if n>=28
    cell1 = load('cell_14l08.mat');
    it.twentyeigth = cell1.data;
end
if n>=29
    cell1 = load('cell_14l09.mat');
    it.twentyninth = cell1.data;
end
if n>=30
    cell1 = load('cell_14l10.mat');
    it.thirtieth = cell1.data;
end 
if n>=31
    cell1 = load('cell_14l11.mat');
    it.thirtyfirst = cell1.data;
end
if n>=32
    cell1 = load('cell_14l12.mat');
    it.thirtysecond = cell1.data;
end
if n>=33
    cell33 = load('cell_14l15.mat');
    it.thirtythird = cell33.data;
end
if n>=34
    cell33 = load('cell_14l16.mat');
    it.thirtyfourth = cell33.data;
end
if n>=35
    cell33 = load('cell_14l19.mat');
    it.thirtyfifth = cell33.data;
end
if n>=36
    cell33 = load('cell_12r15.mat');
    it.thirtysixth = cell33.data;
end
if n>=37
    error('Only 33 possible cells... for now')
end
