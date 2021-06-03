deParameters.I_bnd_constr =2; %Using bound constraints /is possible to change direct in DE
% 1 repair to the lower or upper violated bound
% 2 rand value in the allowed range
% 3 bounce back

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings
deParameters.adaptActivated=0; %Change this to one for adaptive DE
deParameters.I_strategy=3; 
% 1 FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
% 2 FM_bm=repmat(FVr_bestmemit,I_NP,1);
% 3: Activates three strategies
deParameters.I_strategyVersion=2;  
    % I_strategyVersion==1; %Emulates Vortex algorithm
    % I_strategyVersion==2; %HyDE-DF
    %  I_strategyVersion==3; %HyDE

deParameters.I_itermax= 999; %For testbed1
deParameters.I_itermax= 9900; %For testbed2

deParameters.I_NP=10;
deParameters.F_weight=0.5; %0.5
deParameters.F_CR=0.5; %0.9