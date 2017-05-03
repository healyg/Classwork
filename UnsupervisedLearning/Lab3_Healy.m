%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Script for Lab 3 (Neural Nets) for Computational Approaches to
% Cognitive Neuroscience, Spring Quarter 2017
% 
% 19APR2017
% 
% Garrett Healy 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

inp = [0.1 0.8 0.1 0.9 0.2 0.7;0.7 0.9 0.8 0.8 0.75 0.9];

ncat = 2;
nfeatures = 2;

w=rand(ncat,nfeatures);
w=w./repmat(sqrt(sum(w.^2,2)),1,nfeatures);
lr=0.05;

figure
axis([0 1 0 1]);
hold on
plot(inp(1,:),inp(2,:),'b+');

wout = w;

for i = 1:1000
    plot([0 wout(1,1)], [0 wout(1,2)],'color','r');
    plot([0 wout(2,1)], [0 wout(2,2)],'color','g'); 
    [wout,ind] = train_cl(wout,inp,lr);
end
    
plot([0 wout(1,1)], [0 wout(1,2)],'bo');
plot([0 wout(2,1)], [0 wout(2,2)],'bo'); 

disp(ind);

% dead unit solution
% I added a while loop to deal with the dead unit issue, however one value
% is far enough away that it eludes this solution and attaches to that one
% value while remaining uncompetitive. 

%% Hopfield

clear all

inp = [1 -1;1 -1;-1 1;-1 1];
newinp = inp(1,:);
W=0;
for j=1:size(inp,2)
    W=W+inp(:,j)*inp(:,j)';
end
W=W-diag(diag(W));

p = 0.7;
thr = 0;

newstate = [0 0 0 0]';
for i = 1:1000
    newinp = inp(:,2)+W*newstate;
    newstate = update_hp(newinp,W,newstate,p,thr);
end
%% test
Test=[1 1 0 0;0 0 -1 -1]';
t = Test;tn=t;
t2 = Test(:,2);t2n=t2;
W=0;
for j=1:size(inp,2)
    W=W+inp(:,j)*inp(:,j)';
end
W=W-diag(diag(W));
p = 0.7;

newstate2 = [0 0 0 0;0 0 0 0]';
for i = 1:1000
    tn = t+W*newstate2;
    newstate2 = update_hp(tn,W,newstate2,p,thr);
end
%%
newstate3 = [0 0 0 0]';
for i = 1:1000
    t2n = t2+W*newstate3;
    newstate3 = update_hp(t2n,W,newstate3,p,thr);
end

%% Project (Greebles) Competitive 

clear all
close all 

bad = xlsread('BadGreeblesTraining.xls');
good = xlsread('GoodGreeblesTraining.xls');
test = xlsread('GreeblesTest.xls');

bad = bad';
good = good';
test = test';

inp = [good bad];

figure 
hold on; grid on;
plot3(bad(1,:),bad(2,:),bad(3,:),'+');
plot3(good(1,:),good(2,:),good(3,:),'+');
view(3)
hold off

ncat = 2;
nfeatures = 3;

w=rand(ncat,nfeatures);
w=w./repmat(sqrt(sum(w.^2,2)),1,nfeatures);
lr=2;

wout = w;

hold on 
for i = 1:1000
    plot3([0 wout(1,1)], [0 wout(1,2)], [0 wout(1,3)],'r');
    plot3([0 wout(2,2)], [0 wout(2,2)], [0 wout(2,3)],'g'); 
    view(3);
    [wout,ind] = train_cl(wout,inp,lr);
end
hold off

one = length(find(ind==1));
two = length(find(ind==2));
disp(one/two);

% test

inpt = test;
out = wout*inpt;
[mx,indt]=max(out);

onet = length(find(indt==1));
twot = length(find(indt==2));
disp(onet/twot);

%% test classification
inpt = test;

figure 
hold on; grid on;
plot3(test(1,:),test(2,:),test(3,:));
view(3)
hold off

ncat = 2;
nfeatures = 3;

wt=rand(ncat,nfeatures);
wt=wt./repmat(sqrt(sum(wt.^2,2)),1,nfeatures);
lr=0.05;

wout = wt;

hold on 
for i = 1:100
    plot3([0 wout(1,1)], [0 wout(1,2)], [0 wout(1,3)],'r');
    plot3([0 wout(2,2)], [0 wout(2,2)], [0 wout(2,3)],'g'); 
    view(3);
    [wout,indt] = train_cl(wout,inpt,lr);
end
hold off

onet = length(find(indt==1));
twot = length(find(indt==2));
disp(onet/twot);

 
%% Greebles Hopfield

m = max(good,[],2);
good = good./m;
m = max(bad,[],2);
bad = bad./m;
m = max(test,[],2);
test = test./m;testn=test;

avgood = mean(good,2);
avbad = mean(bad,2);
avfeats = [avgood avbad];
avfeats = [1 -1;-1 1;1 -1];

W=0;
for j=1:size(avfeats,2)
    W=W+avfeats(:,j)*avfeats(:,j)';
end
W=W-diag(diag(W));

p = 0.5;
thr = 0;

greeblestate = zeros(3,length(test));
for i = 1:5000
    testn = test+W*greeblestate;
    greeblestate = update_hp(testn,W,greeblestate,p,thr);
end

gs = zeros(1,length(test));
for i = 1:length(greeblestate)
    if isequal(greeblestate(:,i),avfeats(:,1))
        gs(i) = 1;
    elseif isequal(greeblestate(:,i),avfeats(:,2))
        gs(i) = 2;
    else 
        gs(i)=0;
    end
end

one = length(find(gs==1));
two = length(find(gs==2));
disp(one/two);


% It appears, through both the competitive learning and the Hopfield
% modules, that Greebles are more likely to be good than bad. The models
% differ in their estimation of this, though. The competitive model claims
% the ratio is 4:1 good to bad, while the Hopfield model oscillates around
% this number, but does not predict it exactly as it never fully reaches
% equilibrium. 
% I do not have very much confidence in this classification as a means of
% determining whether someone is good or bad - however I am confident it
% can point in the right direction to a distinction, while maybe not going all the way. 















