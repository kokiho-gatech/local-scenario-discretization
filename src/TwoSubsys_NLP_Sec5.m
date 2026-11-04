%Two Subsystem NLP case
%Author: Ho
  
clear all
% rng(0)
% Note that rng(0) does not guarantee reproducibility because
% timeit function runs the function for an unidentified number of times.
N=2;
parameter.C1=[4;1];
parameter.C2=[5;2];
parameter.Dji=[0,0.1;0.3,0]; %Dji(i,j) is i->j (opposite of the notation)
parameter.D1=[1;1];
parameter.D2min=parameter.D1;
parameter.D2max=[2,4];
parameter.alpha=0.9;


% cases: 0: no analysis (evaluation only); 1: deterministic + fully flexible + proposed method; 2: proposed method only
% sceanrios: number of scenarios per subsystem (vector if multiple cases tested)
% time: 0: timeit not used; 1: timeit used
% eval: 0: evaluation not needed; 1: evaluation needed
% x1,x,costref: x1 for evaluation (if evaluation only); x for initial guess of evaluation, costref for reference 
cases=1;
scenario=[2,4,8,16];
time=1;
eval=1;
x1=[];
x=[];
costref=[];

[costdet,x1det,xdet,nvardet,tdet,...
    costconvk,x1convk,xconvk,nvarconvk,tconvk,...
    costnewk,x1newk,nvarnewk,tnewk,...
    costnewevalk,gap]=solve(N,parameter,cases,scenario,time,eval,x1,x,costref);
fprintf('Deterministic: cost: %.4f\n',costdet)
fprintf('Deterministic: x1(A): %.4f\n',x1det{1})
fprintf('Deterministic: x1(B): %.4f\n',x1det{2})
fprintf('Deterministic: nvar: %d\n',nvardet)
fprintf('Deterministic: time: %.4f\n',tdet)
fprintf('Fully Flexible: cost: %.4f & %.4f & %.4f & %.4f\n',costconvk(1),costconvk(2),costconvk(3),costconvk(4))
fprintf('Fully Flexible: x1(A): %.4f & %.4f & %.4f & %.4f\n',x1convk{1}(1),x1convk{1}(2),x1convk{1}(3),x1convk{1}(4))
fprintf('Fully Flexible: x1(B): %.4f & %.4f & %.4f & %.4f\n',x1convk{2}(1),x1convk{2}(2),x1convk{2}(3),x1convk{2}(4))
fprintf('Fully Flexible: nvar: %d & %d & %d & %d\n',nvarconvk(1),nvarconvk(2),nvarconvk(3),nvarconvk(4))
fprintf('Fully Flexible: time: %.4f & %.4f & %.4f & %.4f\n',tconvk(1),tconvk(2),tconvk(3),tconvk(4))
fprintf('Proposed Method: cost: %.4f & %.4f & %.4f & %.4f\n',costnewk(1),costnewk(2),costnewk(3),costnewk(4))
fprintf('Proposed Method: x1(A): %.4f & %.4f & %.4f & %.4f\n',x1newk{1}(1),x1newk{1}(2),x1newk{1}(3),x1newk{1}(4))
fprintf('Proposed Method: x1(B): %.4f & %.4f & %.4f & %.4f\n',x1newk{2}(1),x1newk{2}(2),x1newk{2}(3),x1newk{2}(4))
fprintf('Proposed Method: nvar: %d & %d & %d & %d\n',nvarnewk(1),nvarnewk(2),nvarnewk(3),nvarnewk(4))
fprintf('Proposed Method: time: %.4f & %.4f & %.4f & %.4f\n',tnewk(1),tnewk(2),tnewk(3),tnewk(4))
fprintf('Proposed Method: cost evaluated with fully flexible scenarios: %.4f & %.4f & %.4f & %.4f\n',costnewevalk(1),costnewevalk(2),costnewevalk(3),costnewevalk(4))
fprintf('Proposed Method: gap: %.2f \\%%& %.2f \\%%& %.2f \\%%& %.2f \\%%\n',100*gap(1),100*gap(2),100*gap(3),100*gap(4))
save('outputSec5.mat')