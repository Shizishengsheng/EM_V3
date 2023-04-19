% This function is used to generate a random sample of Y
% define the dimension, # of valu
function [Y, tau, nu, mu, Psi] = GenData(p, n, missing_rate)
    mu = rand(p,1);
    Psi = rand(p);
    Psi = Psi * (Psi)';
    %nu = fix(0.9*n + 0.05*n*rand);
    nu = 5;
    tau = gamrnd(nu/2, 2/nu,[1,n]); 
    Y = zeros(p, n);
    for i = 1:n
        Y(:,i) = mvnrnd(mu,Psi/tau(i));
    end
    % adding the missing data
    missing_mask = rand(p, n) <= missing_rate;
    Y(missing_mask) = NaN;
    % adding the missing data
%     for i = 1:p
%         for j = 1:n
%             if rand() <= missing_rate
%                 Y(i,j) = NaN;
%             end
%         end
%     end
end