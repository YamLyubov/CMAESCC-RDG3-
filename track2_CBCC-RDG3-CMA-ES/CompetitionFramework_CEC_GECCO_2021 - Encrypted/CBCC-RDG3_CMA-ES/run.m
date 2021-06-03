clear;
% set random seed
rand('state', sum(100*clock)); 
randn('state', sum(100*clock));
%warning('off' ,'Octave:divide-by-zero');

% number of independent runs
runs = 1;

% number of fitness evaluations
Max_FEs = 10000; % 100000

% for the benchmark functions initialization
global initial_flag;

% load the FEs used by RDG3 in the decomposition process
decResults = sprintf('/MATLAB Drive/Solution_Compilation/rdg3/results/fitFunc');
load (decResults);
FEs = Max_FEs - FEs;

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
    lowerB = lb;
    upperB = ub;

    VTRs = [];
    bestval = [];
    for runindex = 1:runs
    % trace the fitness
    fprintf('Run %02d\n', runindex);
    
    initial_flag = 0;
   

    fnc =  otherParameters.fnc;  
    [val]  = cbcc('fitFunc', 940, lowerB, upperB, FEs);
    bestval = [bestval; val];
    end
    
    
    % Save results
    bestval = sort(bestval);
    % the best results of each independent run
    filename = sprintf('/MATLAB Drive/Solution_Compilation/final_results/best_fitFunc');
    [fid, message] = fopen(filename, 'w');
    fprintf(fid, '%e\n', bestval);
    fclose(fid);

    % mean
    filename = sprintf('/MATLAB Drive/Solution_Compilation/final_results/mean_fitFunc.txt');
    [fid, message] = fopen(filename, 'w');
    fprintf(fid, '%e\n', mean(bestval)); 
    fclose(fid);

    %median
    filename = sprintf('/MATLAB Drive/Solution_Compilation/final_results/median_fitFunc.txt');
    [fid, message] = fopen(filename, 'w');
    fprintf(fid, '%e\n', median(bestval));
    fclose(fid);

    % std
    filename = sprintf('/MATLAB Drive/Solution_Compilation/final_results/std_fitFunc.txt');
    [fid, message] = fopen(filename, 'w');
    fprintf(fid, '%e\n', std(bestval));
    fclose(fid);
    

