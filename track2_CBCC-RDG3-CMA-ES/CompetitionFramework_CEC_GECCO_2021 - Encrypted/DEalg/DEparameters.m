% Author:           Rainer Storn, Ken Price, Arnold Neumaier, Jim Van Zandt
% Modified by FLC \GECAD 04/winter/2017

deParameters.I_NP= 5; % population in DE
deParameters.F_weight= 0.8; %Mutation factor
deParameters.F_CR= 0.9; %Recombination constant
deParameters.I_itermax= 1000; % number of max iterations/gen
deParameters.I_strategy   = 2; %DE strategy

deParameters.I_bnd_constr = 1; %Using bound constraints 
% 1 repair to the lower or upper violated bound 
% 2 rand value in the allowed range
% 3 bounce back

%
