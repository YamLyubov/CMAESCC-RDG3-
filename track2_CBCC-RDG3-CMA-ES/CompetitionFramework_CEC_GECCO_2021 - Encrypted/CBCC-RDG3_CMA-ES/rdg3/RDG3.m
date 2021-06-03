% RDG3 - This function runs the recursive differential grouping 3
%        procedure to decompose a problem
%
% -------
% Inputs:
% -------
%    fun        : the function suite for which the interaction structure 
%                 is going to be identified in this case benchmark problem
%
%    fun_number : the function number.
%
%    options    : this variable contains the options such as problem
%                 dimensionality, upper and lower bounds.
%
% --------
% Outputs:
% --------
%    seps       : a vector of all separable variables.
%    nongroups  : a cell array containing all non-separable groups.
%    FEs        : the total number of fitness evaluations used.


function [seps, nongroups, FEs] = RDG3(fun, caseStudyData, otherParameters, options)
    lb        = options.lbound;
    dim       = options.dim;
    seps      = [];
    nongroups = {};
    FEs       = 0; 
   
    p1  = lb; 
    
    y1  = fun(p1);
    FEs = FEs+1;
    
    xremain = [1:1:dim];
    sub1 = xremain(1);
    sub2 = xremain(2:end);
    
    tn = 50;
      
    while size(xremain,2) > 0 
        xremain = [];
        [sub1_a,xremain,FEs] = INTERACT(fun,sub1,sub2,p1,y1,FEs,xremain,options, caseStudyData, otherParameters);
        if length(sub1_a) ~= length(sub1) && length(sub1_a) < tn             
            sub1 = sub1_a;
            sub2 = xremain;
            if size(xremain,2) == 0
                nongroups = {nongroups{1:end}, sub1};
                break;
            end                  
        else 
            if length(sub1_a) == 1
               seps = [seps;sub1_a];
            else
               nongroups = {nongroups{1:end},sub1_a};
            end 
            
            if length(xremain) > 1
                sub1 = xremain(1);
                xremain(1) = [];
                sub2 = xremain(1:end);
            else if length(xremain) == 1
                seps = [seps;xremain(1)];
                break;
                end                       
            end   
        end
    end 