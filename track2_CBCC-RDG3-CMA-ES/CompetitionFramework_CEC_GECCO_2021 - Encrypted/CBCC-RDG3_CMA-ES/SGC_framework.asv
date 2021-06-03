function [fit_and_p,FVr_bestmemit, fitMaxVector] ...
    = SGC_framework(caseStudyData, otherParameters,low_habitat_limit,up_habitat_limit)
% set random seed
rand('state', sum(100*clock)); 
randn('state', sum(100*clock));

% number of fitness evaluations
Max_FEs = 100000; % should be 100000

% set size dim
lb_eq_ub = find(low_habitat_limit == up_habitat_limit);
D = 940 - length(lb_eq_ub);
lb_Neq_ub = find(low_habitat_limit ~= up_habitat_limit);
% change lb and ub
low_habitat_limit_new = low_habitat_limit(lb_Neq_ub);
up_habitat_limit_new  = up_habitat_limit(lb_Neq_ub);
% There start RDG3 algorithm 
opts.lbound  = low_habitat_limit_new;
opts.ubound  = up_habitat_limit_new;
opts.dim     = D;

function y = fnc(x)
    x_new = up_habitat_limit;
    x_new(lb_Neq_ub) = x;
    [y, ~] = feval(otherParameters.fnc, x_new,caseStudyData,otherParameters);
end


fprintf('Start RDG3\n');
    [seps, nonseps, FEs] = RDG3(@fnc, caseStudyData, otherParameters, opts);
fprintf(1, 'End RDG3 \n Grupping FEs = %d', FEs);


FEs = Max_FEs - FEs;

% trace the fitness 
[solution]  = cbcc(@fnc, D, low_habitat_limit_new, up_habitat_limit_new, FEs, caseStudyData, otherParameters, seps, nonseps);
fit_and_p = solution.fit_and_p;
Fvr = up_habitat_limit;
Fvr(lb_Neq_ub) = solution.sol;
FVr_bestmemit = Fvr;
fitMaxVector = solution.fitVector;

end