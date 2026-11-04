%(1)Deterministic
%(2)Fully Flexible
%(3)Disconnecting weak coupling (using method 1)
%With parameter values from water-electricity-rice nexus (to be refined)
%using LP 

%Author: Ho, Patel



clear all
%rng(0)

%Water Electricity Rice 
%Note: units normalized so the demand stay in similar
%orders of magnitude
%Note: linear costs

parameter.C1 = [2.16; 0.163; 2]; % USD per m^3 , kWh, kg
parameter.C2 = parameter.C1 * 1.06^15;
parameter.D1 = [1385; 2940; 80]; % m^3, kWh, kg
parameter.D2min = parameter.D1;
parameter.D2max = parameter.D1 .* [1.2; 2; 1.4];
parameter.alpha = 1;

parameter.Dji = [0, 0.147, 1e-10;    % dAA dBA dCA
                 0.057, 0, 1e-10;  % dAB dBB dCB
                 2.5, 0.05, 0];     % dAC dBC dCC

ThreeSubsystemSub(parameter,1);



% parameter.Dji = [0, 0.34, 1e-10;    % dAA dBA dCA
%                  0.057, 0, 1e-10;  % dAB dBB dCB
%                  0.25, 0.05, 0];     % dAC dBC dCC 
% 
% 
% 
% parameter.Dji = [0, 0.34, 1e-10;    % dAA dBA dCA
%                  0.027, 0, 1e-10;  % dAB dBB dCB
%                  2.5, 0.05, 0];     % dAC dBC dCC 
% 


