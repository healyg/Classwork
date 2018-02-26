% pr29_2

clear;
close all;

C=1e-10;
R=2e8;
RL=2e7;
L=2e6;
gam=(1/R)+(RL*C/L);
eps=(1/L)*(1+(RL/R));
B=[1 RL/L];
A=[C gam eps];
w=0:1000;
freqs(B,A,w)
