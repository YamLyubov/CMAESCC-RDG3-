function [solution] = cbcc(func,dim,Lbound,Ubound,FEMax, caseStudyData, otherParameters, seps, nonseps)

    % trace for the fitness value
    groupingFE = 100000 - FEMax; % FEs used in decomposition
    
    countFE = 0;
    maxFECycle = 50;
  
    % grouping cell
    [allGroups] = grouping(seps, nonseps);
    numGroups = size(allGroups, 2);
        
    % Initialization for CMA-ES
    xmean = Lbound+(Ubound-Lbound)/2;    
    sigma = 0.3*(Ubound(1,:)-Lbound(1,:)); % coordinate wise standard deviation (step size)
    
    xbest=xmean;
    bestval= func(xmean);
    
    solution.fit_and_p = [bestval 0];
    solution.sol = xmean;
    solution.fitVector = [bestval];
    
    
    options=cell(numGroups,1);
    for i = 1:numGroups                                      
        % Initialize dynamic (internal) strategy parameters and constants
        options{i}.subpc = zeros(length(allGroups{i}),1);
        options{i}.subps = zeros(length(allGroups{i}),1);
        options{i}.subB = eye(length(allGroups{i}),length(allGroups{i}));
        options{i}.subD = ones(length(allGroups{i}),1);
        options{i}.subC = options{i}.subB * diag(options{i}.subD.^2) * options{i}.subB';
        options{i}.subinvsqrtC = options{i}.subB * diag(options{i}.subD.^-1) * options{i}.subB';
        options{i}.subeigeneval = 0;
        options{i}.subcounteval = 0;
    end
    
    
    % Initialization of fitness improvement array
    fitImpAccum = zeros(1,1*numGroups);
    
    % optimization using CBCC  
    cycle = 0;
    display = 1; % 1 display results; 0 not display results 
    
    while (countFE < FEMax)
        reminder = 100000 - (groupingFE + countFE);
        cycle = cycle + 1;   
        maxVal = max(fitImpAccum);
        idx    = find(fitImpAccum == maxVal);  
        for i = 1:length(idx)    
            groupIdx = idx(i);     
            dimIdx = allGroups{groupIdx};
             FEused = 0;
             [xbestnew,bestvalnew,xmean,sigma,FEused,options{groupIdx}] = ...
             cmaes(func,dim,dimIdx,xbest,xmean,sigma,Lbound(1,:),... 
             Ubound(1,:),maxFECycle,options{groupIdx}, caseStudyData, otherParameters,reminder);       
            if bestvalnew < bestval
                fitimp = (bestval-bestvalnew)/bestval;
                fitImpAccum(idx(i)) = (fitImpAccum(idx(i))+fitimp)/2;   
                xbest = xbestnew;
                bestval = bestvalnew;
            else
                fitImpAccum(idx(i)) = fitImpAccum(idx(i))/2;
            end
            
            countFE = countFE + FEused;                       
        end
        solution.sol = xbest;
        solution.fit_and_p = [bestval 0];
        solution.fitVector = [solution.fitVector bestval];
        if(display == 1 && mod(cycle,1)==0)
           fprintf(1, 'Cycle = %d, bestval = %e, fitImp = %e, component = %d, \n', cycle, bestval, fitImpAccum(idx(i)), groupIdx);
           fprintf('FEs = %d', countFE + groupingFE);
        end
        
        
    end
    fprintf(1, 'Cycle = %d, bestval = %e, component = %d, \n', cycle, bestval, groupIdx);
end