
% Three-subsystem nexus subfunction
% Inputs:
% parameter: parameters
% out: 1: save; 0: not save
% Outputs:
% gappercentlocal, gappercentpartition: optimality gap in percent for local scenario discretization and the best partitioning cases, respectively
% count: count of successful partitioning cases
function [gappercentlocal,gappercentpartition,count]=ThreeSubsystemSub(parameter,out)

N=3;
cases=1;
eval=1;
x1=[];
x=[];
costref=[];
scenario=8;
time=0;

parameter.C1orig=parameter.C1;
parameter.C2orig=parameter.C2;
parameter.D1orig=parameter.D1;
parameter.D2minorig=parameter.D2min;
parameter.D2maxorig=parameter.D2max;
parameter.Djiorig=parameter.Dji;


[costdet,x1det,xdet,nvardet,tdet,...
    costconvk,x1convk,xconvk,nvarconvk,tconvk,...
    costnewk,x1newk,nvarnewk,tnewk,...
    costnewevalk,gap,feasible,deltaji]=solve(N,parameter,cases,scenario,time,eval,x1,x,costref);

x0=xconvk{1};
x1conv=zeros(N,1);
x1new=zeros(N,1);
for i=1:N
    x1conv(i)=x1convk{i};
    x1new(i)=x1newk{i};
end


CSAB_C=deltaji(3,1)+deltaji(3,2)+deltaji(2,3)+deltaji(1,3);
CSBC_A=deltaji(1,2)+deltaji(1,3)+deltaji(2,1)+deltaji(3,1);
CSAC_B=deltaji(2,1)+deltaji(2,3)+deltaji(1,2)+deltaji(3,2);
CSA_B_C=deltaji(1,2)+deltaji(1,3)+deltaji(2,1)+deltaji(2,3)+deltaji(3,1)+deltaji(3,2);


%% Disconnecting weak coupling

parameter.C1orig=parameter.C1;
parameter.C2orig=parameter.C2;
parameter.D1orig=parameter.D1;
parameter.D2minorig=parameter.D2min;
parameter.D2maxorig=parameter.D2max;
parameter.Djiorig=parameter.Dji;

%AB-C:
x1AB_C=zeros(3,1);
nvarnewkAB_C=zeros(2,1);
sub=1:2;
[x1AB_C(1:2),nvarnewkAB_C(1)]=solvesub(sub,parameter,scenario,xdet);

sub=3;
[x1AB_C(3),nvarnewkAB_C(2)]=solvesub(sub,parameter,scenario,xdet);

%evaluate
cases=0;
eval=1;
[~,~,~,~,~,...
    ~,~,~,~,~,...
    ~,~,~,~,...
    costnewevalkAB_C,gapAB_C,feasibleAB_C]=solve(N,parameter,cases,scenario,time,eval,x1AB_C,x0,costconvk);


%BC-A:
x1BC_A=zeros(3,1);
nvarnewkBC_A=zeros(2,1);
sub=2:3;
[x1BC_A(2:3),nvarnewkBC_A(1)]=solvesub(sub,parameter,scenario,xdet);

sub=1;
[x1BC_A(1),nvarnewkBC_A(2)]=solvesub(sub,parameter,scenario,xdet);

%evaluate
cases=0;
eval=1;
[~,~,~,~,~,...
    ~,~,~,~,~,...
    ~,~,~,~,...
    costnewevalkBC_A,gapBC_A,feasibleBC_A]=solve(N,parameter,cases,scenario,time,eval,x1BC_A,x0,costconvk);

%AC-B:
x1AC_B=zeros(3,1);
nvarnewkAC_B=zeros(2,1);
sub=[1,3];
[x1AC_B([1,3]),nvarnewkAC_B(1)]=solvesub(sub,parameter,scenario,xdet);

sub=2;
[x1AC_B(2),nvarnewkAC_B(2)]=solvesub(sub,parameter,scenario,xdet);

%evaluate
cases=0;
eval=1;
[~,~,~,~,~,...
    ~,~,~,~,~,...
    ~,~,~,~,...
    costnewevalkAC_B,gapAC_B,feasibleAC_B]=solve(N,parameter,cases,scenario,time,eval,x1AC_B,x0,costconvk);

%A-B-C
x1A_B_C=zeros(3,1);
nvarnewkA_B_C=zeros(3,1);
sub=1;
[x1A_B_C(1),nvarnewkA_B_C(1)]=solvesub(sub,parameter,scenario,xdet);

sub=2;
[x1A_B_C(2),nvarnewkA_B_C(2)]=solvesub(sub,parameter,scenario,xdet);

sub=3;
[x1A_B_C(3),nvarnewkA_B_C(3)]=solvesub(sub,parameter,scenario,xdet);


%evaluate
cases=0;
eval=1;
[~,~,~,~,~,...
    ~,~,~,~,~,...
    ~,~,~,~,...
    costnewevalkA_B_C,gapA_B_C,feasibleA_B_C]=solve(N,parameter,cases,scenario,time,eval,x1A_B_C,x0,costconvk);

fprintf('CS: -- & -- & %.4f & %.4f & %.4f & %.4f\n',CSAB_C,CSBC_A,CSAC_B,CSA_B_C);
fprintf('SSmax: -- & -- & %d & %d & %d & %d\n',max(nvarnewkAB_C),max(nvarnewkBC_A),max(nvarnewkAC_B),max(nvarnewkA_B_C));

for i=1:3
    fprintf('x1{%d}: %.4f & %.4f & %.4f & %.4f & %.4f & %.4f\n',i,x1conv(i),x1new(i),x1AB_C(i),x1BC_A(i),x1AC_B(i),x1A_B_C(i));
end
fprintf('Gap: -- &  %.4f\\%% &  %.2f\\%% &  %.2f\\%% &  %.2f\\%% &  %.2f\\%% \n\n',100*gap*feasible,100*gapAB_C*feasibleAB_C,100*gapBC_A*feasibleBC_A,100*gapAC_B*feasibleAC_B,100*gapA_B_C*feasibleA_B_C)

gappercentlocal=100*gap*feasible;

CSvec=[CSAB_C,CSBC_A,CSAC_B,CSA_B_C];
[min1, idx1] = min(CSvec);
gapvec=[100*gapAB_C*feasibleAB_C,100*gapBC_A*feasibleBC_A,100*gapAC_B*feasibleAC_B,100*gapA_B_C*feasibleA_B_C];
[min1, idx2] = min(gapvec);

gappercentpartition=gapvec(idx1);

if abs(gapvec(idx1)-min(gapvec))<0.01 %chosen partitioning option returns best gap within tolerance
    count=1;
else
    count=0;
end
if out==1
    save('output.mat');
end
end


% repeat solve for each corresponding subproblems
% Input:
% sub: subproblems
% parameter: parameters
% time: 0: timeit not used; 1: timeit used
% Output:
% x1: 1st stage decision for the proposed method
% nvarnewk: numebr of variables
function [x1,nvarnewk]=solvesub(sub,parameter,scenario,xdet)
parameter.C1=parameter.C1orig(sub);
parameter.C2=parameter.C2orig(sub);
parameter.D1=parameter.D1orig(sub);
parameter.D2min=parameter.D2minorig(sub);
parameter.D2max=parameter.D2maxorig(sub);
parameter.Dji=parameter.Djiorig(sub,sub);
k=0;
N=length(parameter.C1orig);
for i=1:N
    if ismember(i,sub)==1
        k=k+1;
        for j=1:N
            if ismember(j,sub)==0
                %offset demand to account for coupling j(outside of sub)->i(in sub) using xdet
                parameter.D1(k)=parameter.D1(k)+parameter.Djiorig(j,sub(k))*xdet(3*(j-1)+1);
                parameter.D2min(k)=parameter.D2min(k)+parameter.Djiorig(j,sub(k))*xdet(3*(j-1)+3);
                parameter.D2max(k)=parameter.D2max(k)+parameter.Djiorig(j,sub(k))*xdet(3*(j-1)+3);
            end
        end
    end
end
cases=2;
time=0;
eval=0;
x=[];
x1=[];
costref=[];
[~,~,~,~,~,...
    ~,~,~,~,~,...
    ~,x1newk,nvarnewk,~,...
    ~,~]=solve(length(sub),parameter,cases,scenario,time,eval,x,x1,costref);
x1=zeros(length(sub),1);
for i=1:length(sub)
    x1(i)=x1newk{i};
end
end
