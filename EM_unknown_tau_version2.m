%% part2 MLE of (mu, Psi) with unknown tau__version2
clear;clc;close all;
%code start
p=10;
n=1000;
missing_rate = 0;
max_iteration_times = 100;
[Y, tau, nu, mu, Psi] = GenData(p, n, missing_rate);

[mu_hat,Psi_hat,q_store,f_store,mu_store,Psi_store] = EM_unknown_tau_v2(Y, nu, max_iteration_times);
error_mu = norm((mu-mu_hat),'fro')
error_Psi = norm((Psi-Psi_hat),'fro')
disp(mu_store)
disp(Psi_store)


% [Y2] = MVTrand(n,mu,Psi,nu);
% [mu_hat2,Psi_hat2,q_store,f_store] = EM_unknown_tau_v2(Y2', nu, max_iteration_times);
% error_mu2 = norm((mu-mu_hat2),'fro')
% error_Psi2 = norm((Psi-Psi_hat2),'fro')

function [mu_hat,Psi_hat,q_store,f_store,mu_store,Psi_store] = EM_unknown_tau_v2(Y, nu, max_iteration_times)
    %init
    [p,n] = size(Y);
    mu_k = rand(p,1);
    Psi_k = rand(p);
    Psi_k = Psi_k * Psi_k';
    %draw the f function and q function figure
    t = -2:0.1:2;
    dir = rand;
    f_store = [];
    q_store = [];
    mu_store = [];
    Psi_store = [];
    for iter = 1:max_iteration_times
        %applyingEM to do estimation for paremeters
        [S_tau,S_tau_Y,S_tau_Y_Y] = calculateStatistics(Y, mu_k, Psi_k, nu);
        % mu_hat can be consider as mu_{k+1}
        % same to Psi_hat
        mu_hat = S_tau_Y / S_tau;
        Psi_hat = (S_tau_Y_Y - S_tau_Y * S_tau_Y' / S_tau) / n;
        
        %store the q, f value for figure drawing
%         f = calculateF(t,dir,mu_hat,Psi_hat,nu,Y);
        f = calculateF_V2(t,dir,mu_hat,Psi_hat,nu,Y);
        f = log10(f);
        f_store = [f_store f'];

%         f_2 = calculateF_V3(t,dir,mu_hat,Psi_hat,nu,Y);
%         f_2 = log10(f_2);
%         f_store_2 = [f_store_2 f_2'];

        q = calculateQ_V2(t,dir,mu_k,Psi_k,nu,Y);
        q = log10(q);
        idx = find(t == 0);
        q = q + (f(idx) - q(idx));
        q_store = [q_store q'];
        % check the convegency
        if norm((mu_hat - mu_k),'fro') <= 0.001 && norm((Psi_hat - Psi_k),'fro') <= 0.001
            fprintf('iteration ends in the %d-th round.\n',iter)
            break
        end
    
        % update the value that will be used in next iteration
        mu_k = mu_hat; 
        mu_store = [mu_store norm(mu_k)];
        Psi_k = Psi_hat; 
        Psi_store = [Psi_store norm(Psi_k,"fro")];
    end

    % when the estimation ends, start plotting with the store vectors
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
%         plot(t, f_store_2(:,k)', 'g', 'LineWidth', 2, 'DisplayName', 'q');
%         hold on;
%         legend('f','q','f_2');
        %rightness check
%         if any(f_store(:,k)' > q_store(:,k)')
%             for i = 1:size(f_store(:,k)',2)
%                 if f_store(i,k) > q_store(i,k)
%                     fprintf('data error for the %d iteration, with i = %d! \n',k,i);
%                     fprintf('f_store(i,k) = %d \n',f_store(i,k));
%                     fprintf('q_store(i,k) = %d \n',q_store(i,k));
%                 end
%             end
%         end
        %fprintf('end of plotting figure:%d! \n',k);
    end
    figure;
    subplot(2,1,1);
    plot(2:1:size(f_store,2),mu_store,'LineWidth', 2, 'DisplayName', '\mu value');
    legend('\mu value');
    title('\mu changes with the number of iterations');
    subplot(2,1,2);
    plot(2:1:size(f_store,2),Psi_store,'LineWidth', 2, 'DisplayName', '\Psi value');
    legend('\Psi value');
    title('\PsiE changes with the number of iterations');
end

function [S_tau,S_tau_Y,S_tau_Y_Y] = calculateStatistics(Y, mu, Psi, nu)
    [p,n] = size(Y);
    S_tau_Y = zeros(p,1);
    S_tau_Y_Y = zeros(p,p);
    % calculate the omega first
    omega = zeros(1,n);
    for i = 1:n
        delta_i = (Y(:,i) - mu)' * inv(Psi) * ((Y(:,i) - mu));
        omega(i) = (nu + p)/ (nu + delta_i);
    end
    % then calculate the statistics
    S_tau = sum(omega);
    for i = 1:n
        S_tau_Y = S_tau_Y + omega(i) * Y(:,i);
        S_tau_Y_Y = S_tau_Y_Y + omega(i) * Y(:,i) * Y(:,i)';
    end
end

% function q = calculateQ(t,dir,mu_k,Psi_k,nu,Y)
%     q = zeros(1,size(t,2));
%     [~,n] = size(Y);
%     for i = 1: size(t,2)
%         % transform the variables from theta to t
%         mu_delta = mu_k + t(i).*(dir.*mu_k);
%         Psi_delta = Psi_k + t(i).*(dir.*Psi_k);
% 
%         % calculation the (theta_k + t * delta) part
%         [S_tau,S_tau_Y,S_tau_Y_Y] = calculateStatistics(Y, mu_delta, Psi_delta, nu);
%         
%         % calculation Q value with new statistics and theta_k
%         const_i = -n/2*log(det(Psi_k));
%         tr = trace(-1/2*inv(Psi_k)*S_tau_Y_Y) + (mu_k')*inv(Psi_k)*S_tau_Y...
%             -1/2*(mu_k')*inv(Psi_k)*mu_k*S_tau;
%         q(i) = const_i+tr;
%     end
% end

% function f = calculateF(t,dir,mu,Psi,nu,Y)
%     f = zeros(1,size(t,2));
%     [p, n] = size(Y);
%     for i = 1: size(t,2)
%         % transform the variables from theta to t
%         mu_delta = mu + t(i).*(dir.*mu);
%         Psi_delta = Psi + t(i).*(dir.*Psi);
%         tmp = -n/2*log(det(Psi_delta));
%         sum_y = 0;
%         for j = 1:n
%             log_i = log(1+(1/nu)*(Y(:,j)-mu_delta)' * inv(Psi_delta) * (Y(:,j)-mu_delta));
%             sum_y  = sum_y -((nu+p)/2) *log_i;
%         end
%         f(i) = tmp+sum_y;
%     end
% end 

% calculate the majorization function
% step1. calculate the \theta + t * delta
% step2. calculate the tau within current k-th round (calculation_omega with \theta + t * delta)
% step3. calculate the pdf of each round
% step4. then current q value is sum(log(pdf))
function q = calculateQ_V2(t,dir,mu_k,Psi_k,nu,Y)
    q = zeros(1,size(t,2));
    for i = 1: size(t,2)
        % transform the variables from theta to t
        mu_delta = mu_k + t(i).*(dir.*mu_k);
        tmp = t(i).*(dir.* Psi_k);
        Psi_delta = Psi_k + tmp' * tmp;
%         Psi_delta = Psi_k + t(i).*(dir.*Psi_k);
%         Psi_delta = Psi_delta' * Psi_delta;
        omega = calculation_omega(Y,mu_delta,Psi_delta,nu);
        pdf = mvn_pdf(Y, mu_k, Psi_k, omega);
        q(i) = sum(log(pdf));
    end
end
function omega = calculation_omega(Y,mu,Psi,nu)
    [p,n] = size(Y);
    omega = zeros(1,n);
    for i = 1:n
        delta_i = (Y(:,i) - mu)' * inv(Psi) * ((Y(:,i) - mu));
        omega(i) = (nu + p)/ (nu + delta_i);
    end
end
function pdf = mvn_pdf(Y, mu, Psi, tau)
    [p,n] = size(Y);
    pdf = zeros(1,n);
    for i =1:n
        const = (2*pi)^(-p/2);
        pdf(i) =  const * det(Psi/tau(i))^(-1/2) * exp(-0.5*(Y(:,i)-mu)'*inv(Psi/tau(i))*(Y(:,i)-mu));
    end
end

% calculate the majorization function
% step1. calculate the \theta + t * delta
% step2. calculate the pdf of each round
% step3. then current f value is sum(log(pdf));
function f = calculateF_V2(t,dir,mu,Psi,nu,Y)
    f = zeros(1,size(t,2));
    for i = 1: size(t,2)
        % transform the variables from theta to t
%         mu_delta = rand(p,1);
        mu_delta = mu + t(i).*(dir.*mu);
        tmp = t(i).*(dir.*Psi);
        Psi_delta = Psi + tmp' * tmp;
%         Psi_delta = Psi + t(i).*(dir.*Psi);
%         Psi_delta = Psi_delta' * Psi_delta;

%         if ~(norm((Psi_delta - Psi_delta'),'fro') <= 0.00001)
%             disp(Psi_delta)
%             fprintf('Psi is not symmetic when t = %d\n',t(i))
%             fprintf('current dir = %d\n',dir)
%         end
%         if ~all(eig(Psi_delta) > 0)
%             disp('eig of Psi is')
%             disp(eig(Psi))
%             disp('error of Psi_delta!')
%             disp(eig(Psi_delta))
%         end
        f(i) = sum(log(mvt_pdf(Y,mu_delta,Psi_delta,nu)));
    end
end 
function pdf = mvt_pdf(Y,mu,Psi,nu)
    [p, n] = size(Y);
    pdf = zeros(1, n);
    for i = 1:n
        const_i = gamma((nu+p)/2) / (gamma(nu/2) * (nu*pi)^(p/2) * sqrt(det(Psi)));
        delta = (Y(:,i) - mu)' * inv(Psi) * (Y(:,i) - mu);
        pdf(i) = const_i * (1 + delta/nu)^(-(nu+p)/2);
    end
end

% function f = calculateF_V3(t,dir,mu,Psi,nu,Y)
%     f = zeros(1,size(t,2));
%     for i = 1: size(t,2)
%         % transform the variables from theta to t
%         mu_delta = mu + t(i).*(dir.*mu);
%         Psi_delta = Psi + t(i).*(dir.*Psi);
%         f(i) = sum(log(mvt_pdf_2(Y,mu_delta,Psi_delta,nu)));
%     end
% end 
% function pdf = mvt_pdf_2(Y,mu,Psi,nu)
%     [p,n] = size(Y);
%     pdf = zeros(1,n);
%     for i = 1:n
%         const_i = gamma((nu+p)/2)/((pi*nu)^(p/2)*gamma(nu/2));
%         delta = (Y(:,i) - mu)' * inv(Psi) * (Y(:,i) - mu);
%         pdf(i) = const_i * det(Psi)^(-0.5) * 1/(1+delta/nu)^((nu+p)/2);
%     end
% end

