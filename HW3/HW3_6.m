clear, clc;
mkdir 'figure_jpg'
mkdir 'figure_fig'

P = 6;
f = [-0.25 -0.15 0.01 0.1 0.2 0.3].';
sigma_p = [1 3 1 1 1.5 0.5].';
sigma_w = 1;
M = 10;
N = 100;
n_MC = 5000;

load('MSA_HW3_Problem_6.mat');

L = N - M + 1;
f_est_LS = zeros(n_MC,P);
f_est_TLS = zeros(n_MC,P);
for i_MC = 1 : n_MC
    x = sum(sigma_p.*exp(1j*Psi(i_MC,:).').*exp(1j*2*pi*f*[0:N-1]),1).' + V(i_MC,:).';
    R_L = zeros(M,M);
    for i_L = 1 : L
        x_vec = x(i_L:i_L+M-1);
        R_L = R_L + x_vec*x_vec'/L;
    end
    [U,~,~] = svd(R_L);
    Us = U(:,1:P);
    J1 = [eye(M-1) zeros(M-1,1)];
    J2 = [zeros(M-1,1) eye(M-1) ];
    Us1 = J1*Us;
    Us2 = J2*Us;

    Phi_LS = inv(Us1'*Us1)*Us1'*Us2;
    f_est_LS(i_MC,:) = sort(angle(eig(Phi_LS))/2/pi);

    Q = [Us1 Us2];
    [U,~,~] = svd(Q'*Q);
    U12 = U(1:P,P+1:end);
    U22 = U(P+1:end,P+1:end);
    Phi_TLS = -U12*inv(U22);
    f_est_TLS(i_MC,:) = sort(angle(eig(Phi_TLS))/2/pi);
end

fprintf("LS:\n")
fprintf("Var: [%.2e %.2e %.2e %.2e %.2e %.2e]\n", var(f_est_LS,1))
fprintf("Mean: [%.2e %.2e %.2e %.2e %.2e %.2e]\n", mean(f_est_LS,1))
err = f_est_LS - f.';
fprintf("MSE: [%.2e %.2e %.2e %.2e %.2e %.2e]\n", mean(err.^2,1))
figure
subplot(231); hist(f_est_LS(:,1),1000); title("f_1 (LS)")
subplot(232); hist(f_est_LS(:,2),1000); title("f_2 (LS)")
subplot(233); hist(f_est_LS(:,3),1000); title("f_3 (LS)")
subplot(234); hist(f_est_LS(:,4),1000); title("f_4 (LS)")
subplot(235); hist(f_est_LS(:,5),1000); title("f_5 (LS)")
subplot(236); hist(f_est_LS(:,6),1000); title("f_6 (LS)")
jpg_name = sprintf('figure_jpg/HW3_6a.jpg');
fig_name = sprintf('figure_fig/HW3_6a.fig');
saveas(gcf, jpg_name)
saveas(gcf, fig_name)

fprintf("TLS:\n")
fprintf("Var: [%.2e %.2e %.2e %.2e %.2e %.2e]\n", var(f_est_TLS,1))
fprintf("Mean: [%.2e %.2e %.2e %.2e %.2e %.2e]\n", mean(f_est_TLS,1))
err = f_est_TLS - f.';
fprintf("MSE: [%.2e %.2e %.2e %.2e %.2e %.2e]\n", mean(err.^2,1))
figure
subplot(231); hist(f_est_TLS(:,1),1000); title("f_1 (TLS)")
subplot(232); hist(f_est_TLS(:,2),1000); title("f_2 (TLS)")
subplot(233); hist(f_est_TLS(:,3),1000); title("f_3 (TLS)")
subplot(234); hist(f_est_TLS(:,4),1000); title("f_4 (TLS)")
subplot(235); hist(f_est_TLS(:,5),1000); title("f_5 (TLS)")
subplot(236); hist(f_est_TLS(:,6),1000); title("f_6 (TLS)")
jpg_name = sprintf('figure_jpg/HW3_6b.jpg');
fig_name = sprintf('figure_fig/HW3_6b.fig');
saveas(gcf, jpg_name)
saveas(gcf, fig_name)

