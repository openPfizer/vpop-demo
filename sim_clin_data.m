function [T,BG_data,BI_data]= sim_clin_data(N_of_subjects,CV,rsquared,auto_decay)
% Generation of simulated clinical data from "acceptable" profile.
% Additions for data generation by R Allen and T Rieger
% This function generates the simulated data that we will attempt to match.
% We generated simulated profiles by starting from an "acceptable" profile and adding
% a user-specified amount of noise.
%
% Inputs:
%	N_of_subjects - number of individual trajectories ,
%	CV coefficient of variation of marginal distributions, to define sigma
% 	rsquared  - correlation between the two variables
% 	auto-decay - constant for decay rate of auto-correlation of time-series. 
%
%	Outputs:
%		T - vector of times that BG (A1) and BI (B1) are collected at
%		BG_data - Matrix of all generated subjects vs. time for variable A1
%		BI_data - Matrix of all generated subjects vs. time for variable B1
%

% Modified from Duffull et. al:
    %Damping example with Data
    % Using a simple model that resembles an insulin-glucose - 3-state model
    % This produces Figure 1 (panel with acceptable)
    %
    % S Duffull 4th May, 

%% glucose parameters/variables
GlucoseDose = 10;
ka = 1; % absorption rate constant for glucose
t_half_gluc = 1; % hours
k_gluc = log(2)/t_half_gluc; % \hour
BLG = 5; % baseline glucose [mM]
RinGluc = BLG*k_gluc; %Production rate for glucose

%% insulin parameters
t_half_ins = 0.1;  %hours
k_ins = log(2)/t_half_ins; % /hour 
BLI = 1; % normalised baseline value
RinIns = BLI * k_ins; % Production rate for insulin

%% initial conditions
y0 = [GlucoseDose BLG BLI];

%% ode solution - this is just for comparison
options = odeset('RelTol',1e-6);

[t,y]=ode15s(@drugpk,[0 48], [y0], options, ka,k_gluc, RinGluc, BLG, k_ins, RinIns, BLI);

BG = y(:,2);
BI = y(:,3);

%% Created simulated profiles:

sigma_BG = BG(1)*CV;
sigma_BI = BI(1)*CV;
sigma_BIBG = rsquared^0.5*(sigma_BI*sigma_BG);

cov =[sigma_BI^2 sigma_BIBG; sigma_BIBG sigma_BG^2];
T =[0 1  2 4 8 12];

sample = gen_clin_data(t,T,BI,BG,N_of_subjects, cov, auto_decay);

BG_data = squeeze(sample(:,:,2));
BI_data = squeeze(sample(:,:,1));

%% dydt function
function dydt = drugpk(t, y, ka,k_gluc, RinGluc, BLG, k_ins, RinIns, BLI)

GI=1+(y(2)-BLG)^2/((y(2)-BLG)^2+1);
IG=1-(y(3)-BLI)/((y(3)-BLI)+1);

dydt=[-ka*y(1)
    RinGluc*IG+ka*y(1)-k_gluc*y(2)
    RinIns*GI-k_ins*y(3)];
end


%% generate data function
function [sample] = gen_clin_data(t,T,BI,BG,N_of_subjects,cov,auto_decay)

BI_base = interp1(t,BI,T);
BG_base = interp1(t,BG,T);

sigma = cov;
mu = [BI_base; BG_base]';
R = chol(sigma);

sample  = zeros(numel(T), N_of_subjects,2);

for j = 1:N_of_subjects
    for i = 1:numel(T)
        if i ==1
            sample(i,j,:) = mu(i,:) + randn(1,2)*R;
        else
            prev_sample = squeeze(sample(i-1,j,:));
            prev_mu     = squeeze(mu(i-1,:));
            delta = prev_sample'-prev_mu;
            auto_corr = exp(-auto_decay*(T(i)-T(i-1)));
            mu_new   = mu(i,:) + delta*auto_corr; % add the offet from the previous time point with decay
            
            tmp = [-1 -1];
            
            while any(tmp<=0)    
                tmp = mu_new + randn(1,2)*R;
            end
            
            if any(tmp<=0)
                disp(tmp);
            end

            sample(i,j,:) = tmp;
            
        end
    end
end

end

end




