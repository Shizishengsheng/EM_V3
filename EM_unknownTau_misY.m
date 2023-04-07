%% MLE \mu \Psi with missing \tau and missing Y
function [mu_hat, Psi_hat] = EM_unknownTau_misY(Y, nu, iteration_time)
    %initialize
    [p, n] = size(Y);
    Y_ob = Y(:, all(~isnan(Y)));
    %calculate the initial mu and Psi with fully observed value.
    if isempty(Y_ob)
        mu_k = rand(p,1);
        Psi_k = rand(p);
        Psi_k = Psi_k* Psi_k.'; 
    else
        mu_k = nanmean(Y_ob, 2);
        Psi_k = nancov(Y_ob');
%         mu_k = sum(Y_ob')'/size(Y_ob,2);
%         Psi_k = cov(Y_ob');
    end
    
    for iter = 1: iteration_time
        S_tau_Y = zeros(p,1);
        S_tau_Y_Y = zeros(p,p);
        omega = zeros(1,n);
        %step1: calculate the omega
        for i = 1: n
            mask_ob = ~isnan(Y(:,i));
            delta_i = (Y(mask_ob,i)- mu_k(mask_ob))' * inv(Psi_k(mask_ob,mask_ob)) * (Y(mask_ob,i)- mu_k(mask_ob));
            omega(i) = (nu+size(Y(mask_ob,i),i)) / (nu + delta_i);
        end

        %step2:calculate sufficeient staticstics
        S_tau = sum(omega);
        for i = 1:n
            Yi = Y(:,i);
            if(any(isnan(Yi)))
%                 [Yi_hat,Psi_i] = linearRegreesionY(mu_k,Psi_k,Yi,omega(i));
                [Yi_hat,Psi_i] = conditional_mean(mu_k,Psi_k,Yi);
                S_tau_Y = S_tau_Y + omega(i)*Yi_hat;
                S_tau_Y_Y = S_tau_Y_Y + omega(i) * (Yi_hat) * (Yi_hat)' + Psi_i;
            else
                S_tau_Y = S_tau_Y + omega(i) * Yi;
                S_tau_Y_Y = S_tau_Y_Y + omega(i) * (Yi) * (Yi)'; 
            end
        end
        %step 3: calculate the mu_hat, Psi_hat
        mu_hat = S_tau_Y/S_tau;
        Psi_hat = (S_tau_Y_Y - ((S_tau_Y)*(S_tau_Y)')/S_tau)/n;
%         check the rightness of calculation
%         if any(mu_hat > mu_k) || any(any(Psi_hat > Psi_k))
%             fprintf('calculation ERROR in %d-th round!!\n',iter)
%             break
%         end

        % check the convegency
        if norm((mu_hat - mu_k),'fro') <= 0.0001 && norm((Psi_hat - Psi_k),'fro') <= 0.0005
            fprintf('iteration ends in the %d-th round.\n',iter)
            break
        end

        % update the value that will be used in next iteration
        mu_k = mu_hat; 
        Psi_k = Psi_hat; 
    end
    
end

function [y_predict,Psi_i] = conditional_mean(mu_k,Psi_k,Yi)
    mask_ob = ~isnan(Yi);
    mask_mis = isnan(Yi);
    % Find the indice of NaN value
    nan_indices = find(mask_mis);
    if isempty(nan_indices)
        y_predict = Yi;
        return;
    end
    % Fill the known value
    y_predict = zeros(size(Yi)); % Initialize y_predict
    y_predict(mask_ob) = Yi(mask_ob);
    % Calculate the conditional mean and covariance
    mu_mis = mu_k(mask_mis) + Psi_k(mask_mis,mask_ob) * inv(Psi_k(mask_ob,mask_ob))...
        *(Yi(mask_ob)-mu_k(mask_ob));
    Psi_mis = (Psi_k(mask_mis,mask_mis) - Psi_k(mask_mis,mask_ob) * inv(Psi_k(mask_ob,mask_ob))...
        *Psi_k(mask_ob,mask_mis));
    Psi_mis = (Psi_mis + Psi_mis') / 2; 
    %used to calculation S_tau_Y_Y
    Psi_i = zeros(size(Yi,1));
    Psi_i(mask_mis,mask_mis) = Psi_mis;
    % Fill the Nan value
    y_predict(nan_indices) = mu_mis; 
end

function [y_predict,Psi_i] = linearRegreesionY(mu_k,Psi_k,Yi,omega_i)
    mask_ob = ~isnan(Yi);
    mask_mis = isnan(Yi);
    % Find the indice of NaN value
    nan_indices = find(mask_mis);
    if isempty(nan_indices)
        y_predict = Yi;
        return;
    end
    % Fill the known value
    y_predict = zeros(size(Yi)); % Initialize y_predict
    y_predict(mask_ob) = Yi(mask_ob);
    % Calculate the conditional mean and covariance
    mu_mis = mu_k(mask_mis) - Psi_k(mask_mis,mask_ob) * inv(Psi_k(mask_ob,mask_ob))...
        *(Yi(mask_ob)-mu_k(mask_ob));
    Psi_mis = (Psi_k(mask_mis,mask_mis) - Psi_k(mask_mis,mask_ob) * inv(Psi_k(mask_ob,mask_ob))...
        *Psi_k(mask_ob,mask_mis))/omega_i;
    Psi_mis = (Psi_mis + Psi_mis') / 2; 
    %used to calculation S_tau_Y_Y
    Psi_i = zeros(size(Yi,1));
    Psi_i(mask_mis,mask_mis) = Psi_mis;
    % Fill the Nan value
    y_predict(nan_indices) = mvnrnd(mu_mis, Psi_mis); 
end

% function f = calculateF(t,dir,mu,Psi,S_tau,S_tau_Y,S_tau_Y_Y,nu,Y)
%     f = zeros(1,size(t,2));
%     [p, n] = size(Y);
%     for i = 1: size(t,2)
%         mu_delta = mu + t(i).*(dir.*mu);
%         Psi_delta = Psi + t(i).*(dir.*Psi);
%         tmp = -n/2*log(det(Psi_delta));
%         sum_y = 0;
%         for j = 1:n
%             log_i = log(nu+(Y(:,j)-mu_delta)' * inv(Psi_delta) * (Y(:,j)-mu_delta));
%             sum_y  = sum_y + -(nu+p)/2 .*log_i;
%         end
%         f(i) = tmp+sum_y;
%     end
% end
% 
% function q = calculateQ(t,dir,mu,Psi,S_tau,S_tau_Y,S_tau_Y_Y,n)
%     q = zeros(1,size(t,2));
%     for i = 1: size(t,2)
%         mu_delta = mu + t(i).*(dir.*mu);
%         Psi_delta = Psi + t(i).*(dir.*Psi);
%         const_i = -n/2*log(det(Psi_delta));
%         tr = trace(-1/2*inv(Psi_delta)*S_tau_Y_Y) + (mu_delta')*inv(Psi_delta)*S_tau_Y...
%             -1/2*(mu_delta')*inv(Psi_delta)*mu_delta*S_tau;
%         q(i) = const_i+tr;
%     end
% end
