function [s, t_sim, bg_sim, bi_sim]= sim_clin_data_for_mh(k,tout,mud,sigd)
%Additions for data generation by R Allen and T Rieger
% This function is for simulation of the model without generating simulated data with comparison to previously generated
% data distributions. 
%
% Inputs:
%	k - parameter vector for simulation (5x1)
%	tout - vector of simulation times requested back, max(tout) is final time
%	mud - Vector of E[X] for the multivariate data distribution we are scoring against (2x1)
% 	sigd - Matrix of Cov[X] for the multivariate data distribution we are scoring against (2x2)
%
% Outputs:
%	s - score based on position of simulation within the multivariate normal PDF parameterized by mud and sigd
%	t_sim - vector of output times, this vector should be the same as tout
%	bg_sim - vector of A1 (or blood glucose) vs. time, same length as t_sim
%	bi_sim - vector of B1 (or blood insulin) vs. time, same length as t_sim
%

%%Modified from Duffull et. al:
    %Damping example with Data
    % Using a simple model that resembles an insulin-glucose - 3-state model
    % This produces Figure 1 (panel with acceptable)
    %
    % S Duffull 4th May, 


%% glucose parameters/variables
GlucoseDose = 10;
ka          = k(1); % absorption rate constant for glucose
k_gluc      = log(2)/k(2); % \hour
BLG         = k(3); % baseline glucose
RinGluc     = BLG*k_gluc; %Production rate for glucose

%% insulin parameters
%t_half_ins  = k(4);  %hours
k_ins       = log(2)/k(4); % /hour 
BLI         = k(5); % normalised baseline value
RinIns      = BLI * k_ins; % Production rate for insulin

%% initial conditions
y0=[GlucoseDose BLG BLI];

%% ode solution - this is just for comparison
options = odeset('RelTol',1e-6);
[t_sim,y] = ode15s(@drugpk, tout, y0, options, ka, k_gluc, RinGluc, BLG, k_ins, RinIns, BLI);

bg_sim = y(:,2);
bi_sim = y(:,3);

x = [bg_sim;bi_sim]'; % Just for packaging for mvnpdf

s = mvnpdf(x,mud,sigd);


%% dydt function
function dydt = drugpk(~, y,ka, k_gluc, RinGluc, BLG, k_ins, RinIns, BLI)

GI = 1+(y(2)-BLG)^2/((y(2)-BLG)^2+1);
IG = 1-(y(3)-BLI)/((y(3)-BLI)+1);

dydt = [-ka*y(1)
    RinGluc*IG+ka*y(1)-k_gluc*y(2)
    RinIns*GI-k_ins*y(3)];
end

end
