%Computational experiment with random cases

%Author: Ho

%Note: the current solve.m uses linprog for LP, but linprog can fail
%occasionally. Gurobi is recommended for robust results.
%We confirm for loop=1:300, the results betwen linprog and Gurobi match,
%but this is not guaranteed.

clear all

N=3;

k=0;
countite=zeros(N,2);
gaplocalite=zeros(N,1);

%Case 1: all connected 
%Case 2: C decoupled
for decouple=[1,0.1]
    rng(0)
    nloop=1000;
    gappercentlocal=zeros(nloop,1);
    gappercentpartition=zeros(nloop,1);
    countk=zeros(nloop,1);
    gaplocal=0;
    k=k+1;
    for loop=1:nloop
        N=3;
        parameter.C1=rand(N,1);
        parameter.C2=parameter.C1.*(1+rand(N,1));
        parameter.D1=rand(N,1);
        parameter.D2min=parameter.D1;
        parameter.D2max=parameter.D1.*(1+rand(N,1));
        parameter.Dji=0.5*rand(N,N); %Dji(i,j) is i->j (opposite of the notation)
        for i=1:N
            parameter.Dji(i,i)=0;
            for j=1:3
                if i==3||j==3
                    parameter.Dji(i,j)=parameter.Dji(i,j)*decouple; %decouple subsystem 3
                end
            end
        end
        parameter.alpha=1;
        [gappercentlocal(loop),gappercentpartition(loop),countk(loop)]=ThreeSubsystemSub(parameter,0);
        loop
    end
    successrate=mean(countk(:))*100;
    gappercentaverage=mean(gappercentlocal);
    gappercentstd=std(gappercentlocal);
    gappercentpartitionaverage=mean(gappercentpartition);
    gappercentpartitionstd=std(gappercentpartition);
    fprintf('Success Rate for all:  %.4f\\%% \n',successrate)
    fprintf('Full Gap Mean:  %.4f\\%% \n',gappercentaverage)
    fprintf('Full Gap Std:  %.4f\\%% \n',gappercentstd)
    fprintf('Partitioned Gap Mean: %.4f\\%% \n',gappercentpartitionaverage)
    fprintf('Partitioned Gap Std:  %.4f\\%% \n',gappercentpartitionstd)
    filename = sprintf('computationalResults_%.1f.mat',decouple);
    save(filename);
end




