function [bestmem,bestval,xmean,sigma,FECount,options]=cmaes(func, dim, dim_index,bestmem,xmean,sigma,lbound, ubound,oneitermax,options, caseStudyData, otherParameters, reminder)
    
    subdims = length(dim_index);
    subxmean = xmean(dim_index)';
    subsigma = sigma(dim_index)';
  
    % Strategy parameter setting: Selection
    sublambda   = 4+floor(3*log(subdims));  % population size, offspring number
    if sublambda*oneitermax > reminder
        oneitermax = floor(oneitermax/2);
        if sublambda*oneitermax > reminder
            FECount = reminder;
            bestmem = bestmem;
            bestval = func(bestmem);
            xmean = xmean;
            sigma = sigma;
            oprionst = options;
            return
        end
    end
    
   
    submu       = sublambda/2;                          % number of parents/points for recombination
    subweights  = log(submu+1/2)-log(1:submu)';         % muXone array for weighted recombination
    submu       = floor(submu);
    subweights  = subweights/sum(subweights);           % normalize recombination weights array
    submueff    = sum(subweights)^2/sum(subweights.^2); % variance-effectiveness of sum w_i x_i
      
    % Strategy parameter setting: Adaptation
    subcc    = (4 + submueff/subdims) / (subdims + 4 + 2*submueff/subdims); % time constant for cumulation for C
    subcs    = (submueff+2) / (subdims+submueff+5);  % t-const for cumulation for sigma control
    subc1    = 2 / ((subdims+1.3)^2+submueff);    % learning rate for rank-one update of C
%   subcmu   = 2 * (submueff-2+1/submueff) / ((subdims+2)^2+submueff);
    subcmu   = min(1-subc1, 2 * (submueff-2+1/submueff) / ((subdims+2)^2+submueff));  % and for rank-mu update
    subdamps = 1 + 2*max(0, sqrt((submueff-1)/(subdims+1))-1) + subcs; % damping for sigma 
                                                          % usually close to 1                                                    
    subchiN = subdims^0.5*(1-1/(4*subdims)+1/(21*subdims^2));
    
    
    % Initialize dynamic (internal) strategy parameters and constants
    subpc = options.subpc;
    subps = options.subps;
    subB = options.subB;
    subD = options.subD;
    subC = options.subC;
    subinvsqrtC = options.subinvsqrtC;
    subeigeneval = options.subeigeneval;
    subcounteval=options.subcounteval;
    
    
  % ---------------------------------------------------------------------

  % -------------------- Generation Loop --------------------------------
  % Call the purecmaes loop in turn for each subcomponent in each
  % generation

  FECount = 0;  % the next 40 lines contain the 20 lines of interesting code 
  GenCount = 0;
  
  while GenCount < oneitermax
      % Start a new generation
      GenCount = GenCount + 1;
      % Repeat for each subcomponent
      subx = zeros(subdims, sublambda);
      x = zeros(dim, sublambda);
      val = zeros(sublambda, 1);
      % Generate and evaluate lambda offspring
      for k=1:sublambda,
          x(:,k)=bestmem;
          subx(:,k) = subxmean + subsigma .* (subB * (subD.* randn(subdims,1))); % m + sig * Normal(0,C) 
          
          % Checking lb ub
          ubound_dim = ubound(dim_index);
          indexup=find(subx(:,k)>ubound_dim');
          subx(indexup,k)=ubound_dim(indexup)';
          lbound_dim = lbound(dim_index);
          indexlow=find(subx(:,k)<lbound_dim');
          subx(indexlow,k)=lbound_dim(indexlow)';
          
          x(dim_index,k) = subx(:,k); % this is only for grouping{i}
          val(k) = func(x(:,k)');
          FECount = FECount+1;
          subcounteval=subcounteval+1;
      end

      % Sort by fitness and compute weighted mean into xmean
      [val, xindex] = sort(val);  % minimization
      
      subxold = subxmean;
      subxmean = x(dim_index,xindex(1:submu)) * subweights;  % recombination, new mean value

      % Cumulation: Update evolution paths
      subps = (1-subcs) * subps + sqrt(subcs*(2-subcs)*submueff) * subinvsqrtC * ((subxmean-subxold) ./ subsigma); 
      subhsig = sum(subps.^2)/(1-(1-subcs)^(2*subcounteval/sublambda))/subdims < 2 + 4/(subdims+1);
      subpc = (1-subcc) * subpc + subhsig * sqrt(subcc.*(2-subcc)*submueff) * ((subxmean-subxold) ./ subsigma); 

      % Adapt covariance matrix C
      subartmp = repmat(1./subsigma, 1, submu) .* ((x(dim_index,xindex(1:submu)) - repmat(subxold,1,submu)));  % mu difference vectors
      subC = (1-subc1-subcmu) * subC ...                   % regard old matrix  
             + subc1 * (subpc * subpc' ...                % plus rank one update
             + (1-subhsig) * subcc*(2-subcc) * subC) ... % minor correction if hsig==0
             + subcmu * subartmp * diag(subweights) * subartmp'; % plus rank mu update 

      % Adapt step size sigma
      subsigma = subsigma * exp((subcs/subdamps)*(norm(subps)/subchiN - 1)); 

       % Update B and D from C
       if subcounteval - subeigeneval > sublambda/(subc1+subcmu)/subdims/10  % to achieve O(N^2)
          subeigeneval = subcounteval;
          subC = triu(subC) + triu(subC,1)'; % enforce symmetry
          
          subC(subC==Inf)=0;
          subC(isnan(subC))=0;
          
          [subB,subD] = eig(subC);           % eigen decomposition, B==normalized eigenvectors
          subD = sqrt(diag(abs(subD)));        % D contains standard deviations now
          subinvsqrtC = subB * diag(subD.^-1) * subB';
        end
    
  end
  
  bestmem(dim_index)=subx(:,xindex(1));
  bestval=val(1);
  xmean(dim_index)=subxmean';
  sigma(dim_index)=subsigma';
  
  options.subpc = subpc;
  options.subps = subps;
  options.subB  = subB;
  options.subD  = subD;
  options.subC  = subC;
  options.subinvsqrtC  = subinvsqrtC;
  options.subeigeneval = subeigeneval;
  options.subcounteval = subcounteval;  
 end
