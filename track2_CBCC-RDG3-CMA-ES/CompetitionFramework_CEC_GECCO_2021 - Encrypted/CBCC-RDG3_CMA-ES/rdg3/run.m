% Author: Yuan SUN
% Email:  yuan.sun@rmit.edu.au 
%         suiyuanpku@gmail.com
%
% ------------
% Description:
% ------------
% This file is the entry point for running the recursive differential 
% gropuing 3 algorithm for decomposing CEC'2013 benchmark problems.

clear;
func_num = 1;
% set the dimensionality and upper and lower bounds of the search space
D = 940;
DEparameters %Function defined by the participant
No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
DB=2;
[caseStudyData, ~]=callDatabase(DB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set other parameters
otherParameters =setOtherParameters(caseStudyData,No_solutions,DB);
otherParameters.one=setOtherParameters(caseStudyData,1,DB);
[lb,ub] = setVariablesBounds(caseStudyData,otherParameters, DB);
    opts.lbound  = lb;
    opts.ubound  = ub;
    opts.dim     = D;

    global initial_flag;
    initial_flag = 0;
    fprintf('Start' );
    [seps, nonseps, FEs] = RDG3('fitFunc', func_num, opts);
    fprintf('End');
    filename = sprintf('/MATLAB Drive/Solution_Compilation/rdg3/results/fitFunc');
    save (filename, 'seps', 'nonseps', 'FEs','-v7');

