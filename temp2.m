% clear;
% p = 10;
% n = 10000;
% mu = rand(p,1);
% Psi = rand(p);
% Psi = Psi* Psi'; 
% % nu = round((n- 50) * rand()+50);
% nu = fix(0.9*n + 0.05*n*rand);
% tau = gamrnd(nu/2, 2/nu,[1,n]); %?
% Y = zeros(p, n);
% for i = 1:n
%     Y(:,i) = mvnrnd(mu,Psi/tau(i));
% end
% % adding the missing data
% for i = 1:p
%     for j = 1:n
%         if rand() <= 0.15
%             Y(i,j) = NaN;
%         end
%     end
% end
% [mu_hat, Psi_hat] = EM_unknownTau_misY(Y, nu, 200);
% error_mu = norm((mu-mu_hat),'fro')
% error_Psi = norm((Psi-Psi_hat),'fro')
%% new test
% clear;clc;close all;
% p=10;
% n = 1000:1000:10000;
% missing_rate = 0.1;
% error_mu =zeros(size(n));
% error_Psi = zeros(size(n));
% try_times= 20;
% for iter = 1:try_times
%     for i = 1:size(n,2)
%         [Y, tau, nu, mu, Psi] = GenData(p, n(i), missing_rate);
%         [mu_hat, Psi_hat] = EM_unknownTau_misY(Y, nu, 50);
%         error_mu(i) = error_mu(i)+norm((mu-mu_hat),'fro');
%         error_Psi(i) = error_Psi(i)+norm((Psi-Psi_hat),'fro');
%     end
% end
% 
% figure;
% subplot(2,1,1);
% plot(n, error_mu./try_times, 'r', 'LineWidth', 2, 'DisplayName', 'error_mu');
% hold on;
% title('error of estimated \mu v.s. sample numbers n')
% subplot(2,1,2);
% plot(n, error_Psi./try_times, 'b', 'LineWidth', 2, 'DisplayName', 'error_Psi');
% hold on;
% title('error of estimated \Psi v.s. sample numbers n')
% func_collection = FuncCollection;
% Y_ob = func_collection.FindOb(Y);
% figure;
% iter_time = 100:100:1000

clear;clc;close all;
p=50;
n = 1000;
missing_rate = 0.0;
[Y, tau, nu, mu, Psi] = GenData(p, n, missing_rate);
[mu_hat, Psi_hat] = validationOfAlgorithm(Y, nu, 500);




