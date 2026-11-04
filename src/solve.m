%Main solve function
%Problem: defined in the NSF CAREER proposal by Ho with nonliner cost
%Local scenario discretization using local sensitivity analysis
%generalized to multiple subsystems
%Author: Ho
%Crate Date: 4/30/2025
%Modify Date: 9/7/2025


% solve function
% Input:
% N: number of subsystems
% parameter: parameters
% cases: 0: no analysis (evaluation only); 1: deterministic + fully flexible + proposed method; 2: proposed method only
% sceanrios: number of scenarios per subsystem (vector if multiple cases tested)
% time: 0: timeit not used; 1: timeit used
% eval: 0: evaluation not needed; 1: evaluation needed
% x1,x,costref: x1 for evaluation (if evaluation only); x for initial guess of evaluation, costref for reference
% Output:
% costdet,x1det,xdet,nvardet,tdet: Output about deterministic case
% costconvk,x1convk,xconvk,nvarconvk,tconvk: Output for conventional fully flexible case
% costnewk,x1newk,nvarnewk,tnewk: Output about proposed method case
% costnewevalk,gap,feasible: evaluation of new method, optimality gap, and if the solution was feasible or not
% deltaji: relative normalized sensitivity (for partitioning purpose)
function [costdet,x1det,xdet,nvardet,tdet,...
    costconvk,x1convk,xconvk,nvarconvk,tconvk,...
    costnewk,x1newk,nvarnewk,tnewk,...
    costnewevalk,gap,feasible,deltaji]=solve(N,parameter,cases,scenario,time,eval,x1,x,costref)

x1det=cell(N,1);
if cases==1
    % deterministic approach
    [xdet,fval,nvar]=deterministic(parameter,N);
    costdet=fval;
    x1det{1}=xdet(1);
    if N>=2
        x1det{2}=xdet(4);
    end
    if N>=3
        x1det{3}=xdet(7);
    end
    nvardet=nvar;
else
    xdet=[];
    costdet=[];
    x1det=[];
    nvardet=[];
end

% Computational time
if time==1
    f=@()deterministic(parameter,N);
    tdet=timeit(f);
else
    tdet=[];
end

k=0;
xconvk=cell(length(scenario),1);
costconvk=zeros(length(scenario),1);
x1convk=cell(length(scenario),1);
nvarconvk=zeros(length(scenario),1);
tconvk=zeros(length(scenario),1);
costnewk=zeros(length(scenario),1);
x1newk=cell(length(scenario),1);
nvarnewk=zeros(length(scenario),1);
tnewk=zeros(length(scenario),1);
costnewevalk=zeros(length(scenario),1);
gap=zeros(length(scenario),1);

if cases==1
    costref=zeros(length(scenario),1);
end

for S=scenario
    for n=1:N
        parameter.values{n}=linspace(parameter.D2min(n),parameter.D2max(n),S);
%         parameter.values{n}=linspace(parameter.D1(n),parameter.D2max(n),S);
%         parameter.values{n}=linspace(0,parameter.D2max(n),S);
%         mu=(parameter.D2max(n)-1)/2+1;
%         sigma=(parameter.D2max(n)-1)/6;
%         parameter.values{n}=sort(mu + sigma * randn(S,1))';
    end
    k=k+1;
    if cases==1
        % Conventional approach (fully flexible)
        % all variables consider all combinations of these scenarios
        %         tic
        [x1,x,fval,nvar]=fullyflexible(S,parameter,N);
        xconvk{k}=x;
        %         toc
        costconvk(k)=fval;
        costref(k)=fval;
        x1convk{1}(k)=x1(1);
        if N>=2
            x1convk{2}(k)=x1(2);
        end
        if N>=3
            x1convk{3}(k)=x1(3);
        end
        nvarconvk(k)=nvar;
    end
    % Computational time
    if time==1
        f=@()fullyflexible(S,parameter,N);
        tconvk(k)=timeit(f);
    end
    if cases==1||cases==2
        method=1; % proposed method (=1 always)
        %         tic
        [x1,fval,nvar,deltaji]=proposedapproach(S,parameter,N,method);
        %         toc
        costnewk(k)=fval;
        x1newk{1}(k)=x1(1);
        if N>=2
            x1newk{2}(k)=x1(2);
        end
        if N>=3
            x1newk{3}(k)=x1(3);
        end
        nvarnewk(k)=nvar;
    end
    % Computational time
    if time==1
        f=@()proposedapproach(S,parameter,N,method);
        tnewk(k)=timeit(f);
    end

    % evaluate x1 from the propseod method using fully flexible formulation
    if eval==1
        [fval,~,feasible]=evaluate(S,parameter,N,x1,x);
        if feasible==1
            costnewevalk(k)=fval;
            gap(k)=abs(costnewevalk(k)-costref(k))/costref(k);
        else
            costnewevalk(k)=NaN;
            gap(k)=NaN;
        end
    end
end
end


% solving the problem using the deterministic approach
% Input:
% parameter: parameter
% N: number of subsystems
% Output:
% xdet: variables
% fval: cost
% nvar: Number of variables
function [xdet,fval,nvar]=deterministic(parameter,N)
for n=1:N
    parameter.values{n}=parameter.D2max(n);
end
parameter.allcomb=combvec(parameter.values{:});

[fLP,A,b,lb,nvar,ncon]=definefAblb(ones(N,1),parameter,N,ones(N,N));

if parameter.alpha~=1
    f=@(x)costfunction(x,ones(N,1),parameter,N);
    options=optimoptions('fmincon','Display','none','MaxFunctionEvaluations',200*nvar);
    [xdet,fval,exitflag]=fmincon(f,ones(N*3,1),A,b,[],[],lb,[],[],options);
    if exitflag<0
        disp('infeasibility in deterministic case');
    end
else
%     options=optimoptions('linprog','Display','none');
%     [xdet,fval,exitflag]=linprog(fLP,A,b,[],[],lb,[],options);
    [xdet,fval,exitflag]=gurobi_linprog(fLP,A,b,[],[],lb,inf*ones(size(lb,1),1));

    if exitflag<0
        disp('infeasibility in deterministic case');
    end
end
end


% solving the problem using the conventional approach (fully flexible)
% Input:
% S: # of scenarios in A or B demand
% parameter: parameter
% N: number of subsystems
% Output:
% x1,x: variables
% fval: cost
% nvar: Number of variables
function [x1,x,fval,nvar]=fullyflexible(S,parameter,N)

% find initial guess by solving the deterministic optimization
[xdet,fval,nvar]=deterministic(parameter,N);

% solve the fully flexible problem using xdet as initial guess
[x1,x,fval,nvar]=problem(S,parameter,N,xdet);
end


% solving the problem using the proposed approach
% Input:
% S: # of scenarios in A or B demand
% parameter: parameter
% N: number of subsystems
% method: scenario discretization method (1: sensitivity; 2: SAA)
% Output:
% x1: variables
% fval: cost
% nvar: number of variables
% deltaji: relative normalized sensitivity
function [x1,fval,nvar,deltaji]=proposedapproach(S,parameter,N,method)
Sji=eye(N);
deltaji=eye(N);

% find initial guess by solving the deterministic optimization
[xdet,fval,nvar]=deterministic(parameter,N);

for n=1:N
    for m=1:N
        if m~=n
            %normalize with deterministic demand and index_x2sn
            deltaji(m,n)=parameter.Dji(m,n)/parameter.D2max(n)*xdet(index_x2sn(m,ones(N,1)));
            %deltaji is the relative normalized sensitivity of incoming coupling vs. the subsytem's own demand
            Sji(m,n)=min(ceil(deltaji(m,n)*S),S);
        end
    end
end

%optimize with SA scenarios for DS1 and SB scenarios for DS2, but
%SAB scenarios for couplingAB and SBA scnearios for couplingBA
[x1,fval,nvar]=problem_new(S,Sji,parameter,N,xdet);
end

% solving the problem using scenario discretization
% Input:
% S: # of scenarios in each subsystem
% parameter: parameter
% N: number of subsytems
% xdet: deterministic solution (for initial guess)
% Output:
% x1,x: variables
% fval: cost
% nvar: number of variables
function [x1,x,fval,nvar]=problem(S,parameter,N,xdet)

parameter.allcomb=combvec(parameter.values{:});
S_vec=ones(N,1)*S^N;%conventional approach: all variables have S^N scenarios

[fLP,A,b,lb,nvar,ncon]=definefAblb(S_vec,parameter,N,S^N*ones(N,N));

if parameter.alpha~=1
    f=@(x)costfunction(x,S_vec,parameter,N);
    x0=ones(nvar,1);
    for n=1:N
        x0([index_x1n(n,S_vec),index_dx2sn(n,S_vec),index_x2sn(n,S_vec)])=[xdet(3*(n-1)+1);ones(S_vec(n),1)*xdet(3*(n-1)+2);ones(S_vec(n),1)*xdet(3*(n-1)+3)];
    end
    options=optimoptions('fmincon','Display','none','MaxFunctionEvaluations',200*nvar);
    [x,fval,exitflag]=fmincon(f,x0,A,b,[],[],lb,[],[],options);
    if exitflag<0
        disp('infeasibility in fully flexible case');
    end
else
%     options=optimoptions('linprog','Display','none');
%     [x,fval,exitflag]=linprog(fLP,A,b,[],[],lb,[],options);
    [x,fval,exitflag]=gurobi_linprog(fLP,A,b,[],[],lb,inf*ones(size(lb,1),1));
    if exitflag<0
        disp('infeasibility in fully flexible case');
    end
end

x1=zeros(nvar,1);
for n=1:N
    x1(n)=x(index_x1n(n,S_vec));
end
end


% solving the problem using scenario discretization and reduced number of
% coupling variables
% Input:
% S: # of scenarios in each subsystem
% Sji: contribution to # of scenarios by couling from i to j
% parameter: parameter
% N: number of subsystems
% xdet: deterministic solution for initial guess
% Output:
% x1: variables
% fval: cost
% nvar: number of variables
function [x1,fval,nvar]=problem_new(S,Sji,parameter,N,xdet)
parameter.allcomb=zeros(N,S^N);
S_vec=zeros(N,1);
for n=1:N
    S_vec(n)=S*prod(Sji(:,n));
    cycles=prod(Sji(n+1:N,n));
    repeats=prod(Sji(1:n-1,n));
    parameter.allcomb(n,1:S_vec(n))=repmat(repelem(parameter.values{n},repeats),1,cycles);
end

[fLP,A,b,lb,nvar,ncon]=definefAblb(S_vec,parameter,N,Sji);

if parameter.alpha~=1
    f=@(x)costfunction(x,S_vec,parameter,N);
    x0=ones(nvar,1);
    for n=1:N
        x0([index_x1n(n,S_vec),index_dx2sn(n,S_vec),index_x2sn(n,S_vec)])=[xdet(3*(n-1)+1);ones(S_vec(n),1)*xdet(3*(n-1)+2);ones(S_vec(n),1)*xdet(3*(n-1)+3)];
    end
    options=optimoptions('fmincon','Display','none','MaxFunctionEvaluations',200*nvar);
    [x,fval,exitflag]=fmincon(f,x0,A,b,[],[],lb,[],[],options);
    if exitflag<0
        disp('infeasibility in proposed method');
    end
else
%     options=optimoptions('linprog','Display','none');
%     [x,fval,exitflag]=linprog(fLP,A,b,[],[],lb,[],options);
    [x,fval,exitflag] =gurobi_linprog(fLP,A,b,[],[],lb,inf*ones(size(lb,1),1));
    if exitflag<0
        disp('infeasibility in proposed method');
    end
end
x1=zeros(N,1);
for n=1:N
    x1(n)=x(index_x1n(n,S_vec));
end
end


% Matrix to convert an n-by-1 vector into an m-by-1 vector by percentiles
% Input
% m: final vector size
% n: initial vector size
% Output
% matrix: the matrix to be multiplied
function matrix=resizematrix(m,n)
if m==n
    matrix=eye(m);
else
    matrix=zeros(m,n);
    for i=1:m
        idx=ceil(i/m*n);
        matrix(i,idx)=1;
    end
end
end


% define constraints and LP costs
% Input:
% S_vec: scenarios for each subsystem
% parameter: parameter
% N: number of subsystems
% Sji: contribution to # of scenarios by couling from i to j
% Output
% fLP: cost coefficient (if linear)
% A,b: constraints
% lb: lower bound
% nvar: number of variables
% ncon: number of constraints
function [fLP,A,b,lb,nvar,ncon]=definefAblb(S_vec,parameter,N,Sji)

C1=parameter.C1;
C2=parameter.C2;
Dji=parameter.Dji;  %Dji(i,j) is i->j
D1=parameter.D1;
D2s=parameter.allcomb';

% variables: x1A,dx2sA,x2sA,x1B,dx2sB,x2sB
% variable index
% x1n: (n-1)+sum(2*S_vec(1:n-1))+1
% dx2sn: n+sum(2*S_vec(1:n-1))+1:n+sum(2*S_vec(1:n-1))+S_vec(n)
% x2sn: n+sum(2*S_vec(1:n-1))+S_vec(n)+1:n+sum(2*S_vec(1:n))
nvar=N+2*sum(S_vec(:));
ncon=N+2*sum(S_vec(:));

if parameter.alpha~=1
    fLP=0;
else
    %linear cost function
    %C=C1A*x1A+Exp[C2A*dx2sA]+C1B*x1B+Exp[C2B*dx2sB]+...
    fLP=zeros(nvar,1);
    for n=1:N
        fLP(index_x1n(n,S_vec))=C1(n);
        fLP(index_dx2sn(n,S_vec))=C2(n)*ones(S_vec(n),1)/S_vec(n);
        fLP(index_x2sn(n,S_vec))=zeros(S_vec(n),1);
    end
end

%constraints
A=zeros(ncon,nvar);
b=zeros(ncon,1);
lb=-inf*ones(nvar,1);

ncon_temp=0; % number of constraints so far
for n=1:N
    A0=zeros(S_vec(n),nvar);
    %dx2sA=max(x2sA-x1A,0) --converted--> dx2sA>=x2sA-x1A, dx2sA>=0
    %dx2sB=max(x2sB-x1B,0) --converted--> dx2sB>=x2sB-x1B, dx2sB>=0

    A0(1:S_vec(n),index_x1n(n,S_vec))=-ones(S_vec(n),1);
    A0(1:S_vec(n),index_dx2sn(n,S_vec))=-eye(S_vec(n));
    A0(1:S_vec(n),index_x2sn(n,S_vec))=eye(S_vec(n));

    A(ncon_temp+1:ncon_temp+S_vec(n),:)=A0;
    b(ncon_temp+1:ncon_temp+S_vec(n))=zeros(S_vec(n),1);
    lb(index_dx2sn(n,S_vec))=zeros(S_vec(n),1);
    ncon_temp=ncon_temp+S_vec(n); % number of constraints so far
end

for n=1:N
    %x1A>=D1A+DAB*coupling1AB (B->A), where coupling1AB=x1B
    %x1B>=D1B+DBA*coupling1BA (A->B), where coupling1BA=x1A
    A1=zeros(1,nvar);
    for m=1:N
        if n==m
            A1(index_x1n(m,S_vec))=-1;
        else
            A1(index_x1n(m,S_vec))=Dji(m,n);
        end
    end
    A(ncon_temp+1,:)=A1;
    b(ncon_temp+1)=-D1(n);
    ncon_temp=ncon_temp+1; % number of constraints so far
end

for n=1:N
    %x2sA>=D2sA+DAB*coupling2AB (B->A), where coupling2AB=resize(x2sB)
    %x2sB>=D2sB+DBA*coupling2BA (A->B), where coupling2BA=resize(x2sA)
    A2=zeros(S_vec(n),nvar);
    for m=1:N
        if n==m
            A2(:,index_x2sn(m,S_vec))=-eye(S_vec(n));
        else
            %The following option can assumes Si*Sji(m,n) variables are
            %conveyed from m to n (better performance and good for integrated
            %optimization)
            A2(:,index_x2sn(m,S_vec))=Dji(m,n)*resizematrix(S_vec(n),S_vec(m));
            %The following option can assumes only Sji(m,n) variables are
            %conveyed from m to n (worse performance but good for MDO)
%             A2(:,index_x2sn(m,S_vec))=Dji(m,n)*repeat_matrix(S_vec(n)/Sji(m,n),Sji(m,n))*resizematrix(Sji(m,n),S_vec(m));
        end
    end
    A(ncon_temp+1:ncon_temp+S_vec(n),:)=A2;
    b(ncon_temp+1:ncon_temp+S_vec(n))=-D2s(1:S_vec(n),n);
    ncon_temp=ncon_temp+S_vec(n); % number of constraints so far
end
end


% cost function
% Input:
% S_vec: scenarios for each subsystem
% parameter: parameter
% N: number of subsystems
% Output
% f: objective function
function f=costfunction(x,S_vec,parameter,N)

C1=parameter.C1;
C2=parameter.C2;

%cost function
%C=C1A*x1A+Exp[C2A*dx2sA]+C1B*x1B+Exp[C2B*dx2sB]+...
f=0;
for n=1:N
    f=f+costsub(C1(n),x(index_x1n(n,S_vec)),parameter);
    f=f+(ones(S_vec(n),1)/S_vec(n))'*costsub(C2(n),x(index_dx2sn(n,S_vec)),parameter);
end
end


% cost subfunction
% Input:
% C: cost parameter
% dxs: deployment or additional deployment at a stage
% parameter: parameter
% Output
% out: stage cost
function out=costsub(C,dxs,parameter)
out=C*dxs.^parameter.alpha;
end


% variable index subfunctions ([x1A,dx2sA,x2sA,x1B,dx2sB,x2sB. ...])
% Input:
% n: subsystem index
% S_vec: number of scenarios for each subsystem
% Output:
% for x1n: (n-1)+sum(2*S_vec(1:n-1))+1
% for dx2sn: n+sum(2*S_vec(1:n-1))+1:n+sum(2*S_vec(1:n-1))+S_vec(n)
% for x2sn: n+sum(2*S_vec(1:n-1))+S_vec(n)+1:n+sum(2*S_vec(1:n))
function out=index_x1n(n,S_vec)
out=(n-1)+sum(2*S_vec(1:n-1))+1;
end
function out=index_dx2sn(n,S_vec)
out=n+sum(2*S_vec(1:n-1))+1:n+sum(2*S_vec(1:n-1))+S_vec(n);
end
function out=index_x2sn(n,S_vec)
out=n+sum(2*S_vec(1:n-1))+S_vec(n)+1:n+sum(2*S_vec(1:n));
end


% Evaluate x1 from the propseod method using fully flexible formulation
% Input:
% S: # of scenarios in each subsystem
% parameter: parameter
% N: number of subsytems
% x1: variables from proposed method
% xconvk: fully flexible case for initial guess
% Output:
% fval: cost
% nvar: number of variables
% feasible: 1: feasible did not relax; NaN: did relax (for LP only)
function [fval,nvar,feasible]=evaluate(S,parameter,N,x1,xconvk)

parameter.allcomb=combvec(parameter.values{:});
S_vec=ones(N,1)*S^N;%conventional approach: all variables have S^N scenarios

[fLP,A,b,lb,nvar,ncon]=definefAblb(S_vec,parameter,N,S^N*ones(N,N));

%fix x1 values to the solution one wants to evaluate
ub=Inf*ones(length(lb),1);
for n=1:N
    lb(index_x1n(n,S_vec))=x1(n);
    ub(index_x1n(n,S_vec))=x1(n);
end

feasible=1;

if parameter.alpha~=1
    f=@(x)costfunction(x,S_vec,parameter,N);
    %use the conventional fully flexible solution as initial guess
    x0=xconvk;
    options=optimoptions('fmincon','Display','none','MaxFunctionEvaluations',200*nvar);
    [x,fval,exitflag]=fmincon(f,x0,A,b,[],[],lb,ub,[],options);
    %if infeasible,
    if exitflag<0
        disp('infeasibility in evaluation (expected)');
        feasible=NaN;
    end
else
%     options=optimoptions('linprog','Display','none');
%     [x,fval,exitflag]=linprog(fLP,A,b,[],[],lb,ub,options);
    [x,fval,exitflag] =gurobi_linprog(fLP,A,b,[],[],lb,ub);
    %if infeasible,
    if exitflag<0
        disp('infeasibility in evaluation (expected)');
        feasible=NaN;
    end
end
end


% repeat_matrix to make vector A repetaed n times 
% Inputs:
% n: number of repetitions per element
% m: number of elements in vector A
% Outputs:
% M: matrix
function M=repeat_matrix(n, m)
    M=zeros(n * m, m);
    for i=1:m
        M((i-1)*n+(1:n),i)=1;
    end
end
