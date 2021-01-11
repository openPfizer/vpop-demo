function [p_pp,pp_yield] = mh_generate_pps(num_pps,mu,sigma,cv_proposal,tsim,p_bnds)
% function mh_generate_pps
%   Metropolis-Hastings search for plausible patients.
%
% Edited from version published in Rieger et al.
% This function executes the main M-H loop, generating plausible patients
% in accordance with the distribution defined as an input (assumed multivariate normal).
%
%	Inputs:
%		num_pps - number of PPs desired
%		mu - vector of E[X] for the data distribution
%		sigma - covariance matrix for the data distribution
%		cv_proposal - coefficient of variation of the proposal distribution 
%		tsim - vector of times to simulate, must correspond to times captured in mu vector
%		p_bnds - Matrix of parameter bounds, (nparam x 2), col 1 - lb, col 2 - ub, must be positive
%
%	Outputs:
%		p_pp - matrix of found plausible patient parameters (nparam x num_pps)
%		pp_yield - metric of the efficiency of generation
%

% Function to call for simulation of each M-H iteration:
f = @(x)sim_clin_data_for_mh(x,tsim,mu,sigma);

% Initialize the looping variables:
k = 1;
iter = 1;
max_iter = num_pps*1000;

% Initial guess, geomean of the bounds:
p = (p_bnds(:,1).*p_bnds(:,2)).^0.5; % Geomean of bounds
num_p = numel(p);

% Set how wide for the next parameter proposal, too wide = random guessing, 
% too narrow = slow progress
sigma_proposal = p*cv_proposal; % <--- TUNABLE by the input value

%% Main MH loop:

s1 = f(p);
while k <= num_pps && iter <= max_iter
    % Propose new parameter set:
    
    for ii = 1:num_p
        q(ii) = icdf(truncate(makedist('Normal','mu',p(ii),'sigma', ...
            sigma_proposal(ii)) ...
            ,p_bnds(ii,1), p_bnds(ii,2)), rand());
    end
    
    %logq = log10(q);
    s2  = f(q);

    r = rand();
    if r < s2/s1 && s2 > eps
        p_pp(:,k) = q(:);
        p = q;
        %logp = log10(p);
        s1 = s2;
        k = k + 1;
    end
    iter = iter + 1;
end

% Check if we exceeded the preset maximum number of iterations:
if iter > max_iter
    error('Exceeded maximum iterations');
end

pp_yield = (k-1)/(iter-1); % calculate a yield metric

end % function mh_generate_pps
