% Subfunction to solve LP problem using Gurobi
% Inputs:
% f: objective
% A,b,Aeq,beq: inequality constraints, equality constraints
% lb,ub:lower bound and upper bound
% Outputs:
% x: optimal solution
% fval: optimal value
% exitflag: 1 if sucess; -1 if failure

function [x,fval,exitflag]=gurobi_linprog(f,A,b,Aeq,beq,lb,ub)
    model.obj=f;
    model.A=sparse([A;Aeq]); 
    model.rhs=[b;beq];
    model.sense=[repmat('<',size(A,1),1);repmat('=',size(Aeq,1),1)];
    model.lb=lb;
    model.ub=ub;
    model.modelsense='min';
    params.outputflag=0; 

    result=gurobi(model,params);

    x=[];
    fval=[];

    if isfield(result,'x')
        x=result.x;
        fval=result.objval;
        exitflag=1; %success
    else
        exitflag=-1; %failure
    end
end
