function [mu_hat, Psi_hat] = validationOfAlgorithm(Y, nu, iteration_time)
    %initialize
    [p, n] = size(Y);
    Y_ob = Y(:, all(~isnan(Y)));
    %calculate the initial mu and Psi with fully observed value.
    if isempty(Y_ob) || length(Y_ob) == length(Y)
        mu_k = rand(p,1);
        Psi_k = rand(p);
        Psi_k = Psi_k* Psi_k.'; 
    else
        mu_k = mean(Y_ob, 2);
        Psi_k = cov(Y_ob');
    end

    t= -2:0.1:2;
    f_store = [];
    q_store = [];

    for iter = 1: iteration_time
        %E -step
        [S_tau,S_tau_Y,S_tau_Y_Y,~] = calculateStatistics(Y, mu_k, Psi_k, nu);
        %draw Q and F
        if iter <=6
            fprintf("start calculate value, iter = %d\n",iter)
            [fValue,qValue] = calculationFAndQ(Y,Psi_k,mu_k,nu,t);
            f_store = [f_store fValue'];
            q_store = [q_store qValue'];
        end

%         figure;
%         if iter>=2 && iter <=6
%             subplot(5 ,1,iter-1);
%             plot(t,fValue);
%             hold on;
%             plot(t,qValue);
%             hold on;
%         end
        %M-step
        mu_hat = S_tau_Y/S_tau;
        Psi_hat = (S_tau_Y_Y - ((S_tau_Y)*(S_tau_Y)')/S_tau)/n;

         % check the convegency
        if norm((mu_hat - mu_k),'fro') <= 0.001 && norm((Psi_hat - Psi_k),'fro') ...
                <= 0.005 || iter == iteration_time
            fprintf('iteration ends in the %d-th round.\n',iter)
            break
        end

        % update the value that will be used in next iteration
        mu_k = mu_hat; 
        Psi_k = Psi_hat; 
    end

    figure;
    % There are k iterations
%     for k = 1:size(f_store,2)
    for k = 2:5
%         subplot(size(f_store,2),1,k);
        subplot(4,1,k-1);
        plot(t, f_store(:,k)', 'r', 'LineWidth', 2, 'DisplayName', 'f');
        hold on;
        plot(t, q_store(:,k)', 'b', 'LineWidth', 2, 'DisplayName', 'q');
        hold on;
        legend('f','q');
    end
end
function [S_tau,S_tau_Y,S_tau_Y_Y,Y_k] = calculateStatistics(Y, mu_k, Psi_k, nu)
    [p, n] = size(Y);
    Y_k = zeros(p,n);
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
        delta_i = (Yi(mask_ob) - mu_k(mask_ob))' *inv(Psi_k(mask_ob,mask_ob)) ...
            *(Yi(mask_ob) - mu_k(mask_ob));
        omega(i) = (nu + length(Yi(mask_ob))) / (nu + delta_i);
        omega(i) = omega(i);

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
            Y_k(:,i) = Yi_hat;
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
    % S_tau_Y_Y = S_tau_Y_Y;
end

function [fValue,qValue] = calculationFAndQ(Y,Psi,mu,nu,t)
    [p,n] = size(Y);
    fValue = zeros(1,size(t,2));
    qValue = zeros(1,size(t,2));
    for i = 1: size(t,2)
        mu_delta = mu + t(i)*randn(p,1);

%         fprintf("start generating\n")
%         while true
%             % 生成随机矩阵A
%             frac = 0.01;
%             Delta = frac*randn(p, p); % 从正态分布生成随机数填充A的元素
%             % 检查Ψ + t*A的正定性
%             if all(eig(Psi + t(i) * Delta) > 0)
%                 break; % 如果满足正定性条件，退出循环
%             end
%         end
        while true
            L = chol(Psi, 'lower'); % 对Psi进行Cholesky分解
            Delta = L' * randn(p); % 生成随机矩阵
            %%这里应该放到循环外
            Psi_delta= Psi + t(i) * Delta;
            if (all(eig(Psi_delta)) > 0)
                break;
            end
        end
%         fprintf("PD delta generated\n")
        Psi_delta = Psi + t(i) * Delta;
        [S_tau,S_tau_Y,S_tau_Y_Y,Y_k] = calculateStatistics(Y, mu_delta, Psi_delta, nu);
        const = - n / 2 * log(det(Psi_delta));
        tmp = 0;
        omega = zeros(1,n);
        %calculate f
        for j = 1:n
            delta = (Y_k(:,j) - mu_delta)' * inv(Psi_delta) * (Y_k(:,j) - mu_delta);
            tmp = tmp+ log(1 + delta/ nu);
            omega(j) = (nu + p) / (nu + delta);
        end
        fValue(i) = const -(nu + p)/2 * tmp;
        
        qValue(i) = const -1/2*trace(inv(Psi_delta) * S_tau_Y_Y) + trace(mu_delta' * inv(Psi_delta)...
            *S_tau_Y) -1/2*trace(mu_delta'*inv(Psi_delta)*mu_delta*S_tau);
    end
    idx = find(t == 0);
    qValue = qValue + (fValue(idx) - qValue(idx));
    
end






