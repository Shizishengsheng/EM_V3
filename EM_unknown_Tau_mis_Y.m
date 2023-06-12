function [mu_hat, Psi_hat] = EM_unknown_Tau_mis_Y(Y, nu, iteration_time)
    %initialize
    [p, n] = size(Y);
    Y_ob = Y(:, all(~isnan(Y)));
    %calculate the initial mu and Psi with fully observed value.
%     if isempty(Y_ob) || length(Y_ob) == length(Y) || size(Y_ob,2) <=1
%     if length(Y_ob) == length(Y) || size(Y_ob,2) <=1
%         mu_k = rand(p,1);
%         Psi_k = rand(p);
%         Psi_k = Psi_k* Psi_k.'; 
%     else
%         mu_k = mean(Y_ob, 2);
%         Psi_k = cov(Y_ob');
%     end
    mu_k = rand(p,1);
    Psi_k = rand(p);
    Psi_k = Psi_k* Psi_k.'; 
    for iter = 1: iteration_time
        %E -step
        [S_tau,S_tau_Y,S_tau_Y_Y] = calculateStatistics(Y, mu_k, Psi_k, nu);
        %M-step
        mu_hat = S_tau_Y/S_tau;
        Psi_hat = (S_tau_Y_Y - ((S_tau_Y)*(S_tau_Y)')/S_tau)/n;

        if(iter == iteration_time)
            fprintf('not Converged! regenerating start point.\n')
            mu_k = rand(p,1);
            Psi_k = rand(p);
            Psi_k = Psi_k* Psi_k.'; 
            iter = 1;
        end
         % check the convegency
        if norm((mu_hat - mu_k),'fro') <= 0.001 && norm((Psi_hat - Psi_k),'fro') ...
                <= 0.005 
            fprintf('iteration ends in the %d-th round.\n',iter)
            break
        end

        % update the value that will be used in next iteration
        mu_k = mu_hat; 
        Psi_k = Psi_hat; 
    end
end
function [S_tau,S_tau_Y,S_tau_Y_Y] = calculateStatistics(Y, mu_k, Psi_k, nu)
    [p, n] = size(Y);
    S_tau_Y = zeros(p,1);
    S_tau_Y_Y = zeros(p,p);
    Psi_i_cnt = zeros(p,p);
    % firstly, get the expectation value of tau
    omega = zeros(1,n);
    for i = 1:n
        Yi = Y(:,i);
        % get the observed value
        mask_ob = ~isnan(Yi);
        mask_mis = ~mask_ob;
        % calculate delta_i(need check!!!)
        %此处加个判断，如果接近奇异值，则重新随机mu_k
        delta_i = (Yi(mask_ob) - mu_k(mask_ob))' * inv(Psi_k(mask_ob,mask_ob)) ...
            *(Yi(mask_ob) - mu_k(mask_ob));
        omega(i) = (nu + length(Yi(mask_ob))) / (nu + delta_i);

        % get the conditional mean to represents Yi_hat
        nan_indices = find(mask_mis);
        if isempty(nan_indices)
            Yi_hat = Yi;
            Psi_i = zeros(p);
        else
            % Fill the known value
            % Initialize Yi_hat
            Yi_hat = zeros(size(Yi)); 
            Yi_hat(mask_ob) = Yi(mask_ob);
%             mu_mis = mu_k(mask_mis) - Psi_k(mask_mis,mask_ob) * inv(Psi_k(mask_ob,mask_ob)) ...
%                 *(Yi(mask_ob)-mu_k(mask_ob));
            % correct: there should be plus rather than minus
            mu_mis = mu_k(mask_mis) + Psi_k(mask_mis,mask_ob) * inv(Psi_k(mask_ob,mask_ob)) ...
                *(Yi(mask_ob)-mu_k(mask_ob));
            % Fill the Nan value
            Yi_hat(nan_indices) = mu_mis; 
            %calculate Psi_i
            Psi_mis = (Psi_k(mask_mis,mask_mis) - Psi_k(mask_mis,mask_ob) ...
                * inv(Psi_k(mask_ob,mask_ob)) *Psi_k(mask_ob,mask_mis));
            %get the ronud result
            Psi_i = zeros(p);
            Psi_i(mask_mis,mask_mis) = Psi_mis;
        end
        S_tau_Y = S_tau_Y + omega(i) * Yi_hat;
        S_tau_Y_Y = S_tau_Y_Y + omega(i) * Yi_hat * Yi_hat';
        Psi_i_cnt = Psi_i_cnt + Psi_i;
    end

    %secondly, calculate the Statistics
    S_tau = sum(omega);
    S_tau_Y_Y = S_tau_Y_Y + Psi_i_cnt;
end






