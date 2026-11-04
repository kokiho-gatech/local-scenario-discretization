%Systems Partitioning
%Author: Ho



clear all
% rng(0) 
% Note that rng(0) does not guarantee reproducibility because
% timeit function runs the function for an unidentified number of times.
N=3;


% cases: 0: no analysis (evaluation only); 1: deterministic + fully flexible + proposed method; 2: proposed method only
% sceanrios: number of scenarios per subsystem (vector if multiple cases tested)
% time: 0: timeit not used; 1: timeit used
% eval: 0: evaluation not needed; 1: evaluation needed
% x1,x,costref: x1 for evaluation (if evaluation only); x for initial guess of evaluation, costref for reference 

cases=1;
scenario=8;
time=0;
eval=1;
x1=[];
x=[];
costref=[];

parameter.C1=[4;1;2];
parameter.C2=[5;2;3];
parameter.D1=[1;1;1];
parameter.D2min=parameter.D1;
parameter.D2max=[2,4,2];
parameter.Dji=[0,0.1,0.03;0.3,0,0.01;0.03,0.01,0]; %Dji(i,j) is i->j (opposite of the notation)
parameter.alpha=0.9;

ThreeSubsystemSub(parameter,1);
