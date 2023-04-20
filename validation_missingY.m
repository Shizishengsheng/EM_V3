%% test for missing Y
clear;clc;close all;
p=10;
n = 1000:1000:10000;
missing_rate = 0.15;
error_mu =zeros(size(n));
% error_mu2 =zeros(size(n));
error_Psi = zeros(size(n));
% error_Psi2 = zeros(size(n));
try_times= 20;
for i = 1:try_times
    for j = 1:length(n)
        [Y, tau, nu, mu, Psi] = GenData(p, n(j), missing_rate);
        [mu_hat, Psi_hat] = EM_unknown_Tau_mis_Y(Y, nu, 150);
        error_mu(j) = error_mu(j) + norm((mu-mu_hat),'fro');
        error_Psi(j) = error_Psi(j)+ norm((Psi-Psi_hat),'fro');
%         [mu_hat, Psi_hat] = EM_unknownTau_misY(Y, nu, 100);
%         error_mu2(j) = error_mu2(j) + norm((mu-mu_hat),'fro');
%         error_Psi2(j) = error_Psi2(j)+ norm((Psi-Psi_hat),'fro');
    end
end
% for i =1:length(n)
%     tmp_mu =[];
%     tmp_Psi=[];
%     [Y, tau, nu, mu, Psi] = GenData(p, n(i), missing_rate);
%     for j = 1: try_times
%         [mu_hat, Psi_hat] = EM_unknown_Tau_mis_Y(Y, nu, 100);
%         tmp_mu(j) = norm((mu-mu_hat),'fro');
%         tmp_Psi(j) = norm((Psi-Psi_hat),'fro');
%     end
%     error_mu(i) = min(tmp_mu);
%     error_Psi(i) = min(tmp_Psi);
% end
figure;
subplot(2,1,1);
title('errors of \mu v.s. # of samples')
plot(n,error_mu/try_times, 'r', 'LineWidth', 2);
hold on;
% plot(1:length(n),error_mu2/try_times, 'b', 'LineWidth', 2);
hold on;
subplot(2,1,2);
title('errors of \Psi v.s. # of samples')
plot(n,error_Psi/try_times, 'b', 'LineWidth', 2);
hold on;
% plot(1:length(n),error_Psi2/try_times, 'b', 'LineWidth', 2);
hold on;