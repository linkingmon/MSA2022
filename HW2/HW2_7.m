clear, clc;
mkdir 'figure_jpg'
mkdir 'figure_fig'

P = 5;
f = [-0.3 -0.2 0.05 0.1 0.4].';
sigma = [2 1 1.5 1.2 3].';
sigma_w = 1;
M = 15;
N = 200;
n_MC = 1000;

load('MSA_HW2_Problem_7.mat');

x = @(i_MC) (sum(sigma.*exp(1j*Psi(i_MC,:).').*exp(1j*2*pi*f*[0:N-1]),1) + sigma_w*V(i_MC,:)).';

%% Prolbem a
for i_MC = 1 : 3
    figure
    title_name = sprintf('Mag. of x of Monte Carlo sim. %d', i_MC);
    sgtitle(title_name)
    subplot(211)
    plot(real(x(i_MC)), 'linewidth', 1);
    ylabel('Real part Mag.', 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Sample', 'FontSize', 12, 'FontName', 'Arial');
    subplot(212)
    plot(imag(x(i_MC)), 'linewidth', 1);
    ylabel('Imag part Mag.', 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Sample', 'FontSize', 12, 'FontName', 'Arial');
    jpg_name = sprintf('figure_jpg/HW2_7a-%d.jpg', i_MC);
    fig_name = sprintf('figure_fig/HW2_7a-%d.fig', i_MC);
    saveas(gcf, jpg_name)
    saveas(gcf, fig_name)
end

%% Prolbem b
R_s = diag(sigma.^2);
V = exp(1j*2*pi*[0:M-1].'*f.');
R = V*R_s*V' + sigma_w^2*eye(M);
R_hat_mean = zeros(M,M,N-M+1);
for i_MC = 1 : n_MC
    x_MC_R = x(i_MC);
    for i_s =  0 : N-M
        R_hat_mean(:,:,i_s+1) = R_hat_mean(:,:,i_s+1) + x_MC_R(i_s+1:i_s+M)*x_MC_R(i_s+1:i_s+M)'/n_MC;
    end
end
for i_s = 1 : N-M+1
    err_mean(i_s) = norm(R_hat_mean(:,:,i_s) - R, 'fro')/M;
end

R_hat = zeros(M,M,N-M+1,n_MC);
for i_MC = 1 : n_MC
    clc;
    fprintf('Calculate Sample Coorelation Matrix (%d/%d)',i_MC,n_MC)
    x_MC_R = x(i_MC);
    for i_l = 1 : N-M+1
        for i_s =  0 : i_l-1
            R_hat(:,:,i_l,i_MC) = R_hat(:,:,i_l,i_MC) + x_MC_R(i_s+1:i_s+M)*x_MC_R(i_s+1:i_s+M)'/i_l;
        end
    end
end
figure
subplot(2,2,1)
plot(0:N-M, err_mean, 'linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Err Mag.', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Sample', 'FontSize', 12, 'FontName', 'Arial');
title('Err mag. of the est. mean', 'FontSize', 10)
for i_MC = 1 : 3
    for i_l = 1 : N-M+1
        err(i_l) = norm(R_hat(:,:,i_l,i_MC) - R, 'fro')/M;
    end
    subplot(2,2,i_MC+1)
    plot(1:N-M+1, err, 'linewidth', 2);
    set(gca, 'YScale', 'log')
    ylabel('Err Mag.', 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Sample', 'FontSize', 12, 'FontName', 'Arial');
    title_name = sprintf('Err. mag. of the samp. corr. in MC=%d',i_MC);
    title(title_name, 'FontSize', 10)
    saveas(gcf, 'figure_jpg/HW2_7b.jpg')
    saveas(gcf, 'figure_fig/HW2_7b.fig')
end

%% Prolbem c
eig_val = svd(R);
subplot(311)
plot(1:M, eig_val(1:M), 'r*');
ylabel('Eig. val Mag.', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Num. of eig. val', 'FontSize', 12, 'FontName', 'Arial');
title('Eig. val distrubution of the correlation matrix')
% title_name = sprintf('Err. mag. of the sample correlation matrix in MC=%d',i_MC);
% title(title_name, 'FontSize', 16)

L = floor((N-M+1)/2);
eig_val = svd(R_hat(:,:,L,1));
subplot(312)
plot(1:M, eig_val(1:M), 'r*');
ylabel('Eig. val Mag.', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Num. of eig. val', 'FontSize', 12, 'FontName', 'Arial');
title('Eig. val distrubution of the estimated mean')

L = N-M+1;
eig_val = svd(R_hat(:,:,L,1));
subplot(313)
plot(1:M, eig_val(1:M), 'r*');
ylabel('Eig. val Mag.', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Num. of eig. val', 'FontSize', 12, 'FontName', 'Arial');
title('Eig. val distrubution of the estimated mean')
saveas(gcf, 'figure_jpg/HW2_7c.jpg')
saveas(gcf, 'figure_fig/HW2_7c.fig')

%% Problem d
L = N-M+1;
f_grid = -0.5:0.001:0.5;
S_MVDR = zeros(numel(f_grid),n_MC);
v = @(f) (exp(1j*2*pi*f*[0:M-1].'));
for i_MC = 1 : n_MC
    clc;
    fprintf('Calculate MVDR (%d/%d)',i_MC,n_MC)
    for i_f = 1 : numel(f_grid)
        S_MVDR(i_f, i_MC) = M / real(v(f_grid(i_f))'*inv(R_hat(:,:,L,i_MC))*v(f_grid(i_f)));
    end
end
figure;
for i_MC = 1 : 3
    [pks,locs] = findpeaks(S_MVDR(:,i_MC));
    [sort_pks, idx] = sort(pks,'descend');
    sort_locs = locs(idx);
    subplot(2,2,i_MC)
    plot(f_grid, S_MVDR(:,i_MC), 'linewidth', 2);
    set(gca, 'YScale', 'log')
    hold on;
    plot(f_grid(sort_locs(1:P)),sort_pks(1:P), 'ro', 'MarkerSize', 8, 'LineWidth', 2.5);
    for i_P = 1 : P
        text(f_grid(sort_locs(i_P))+0.03,sort_pks(i_P),['(' num2str(f_grid(sort_locs(i_P))) ',' num2str(sort_pks(i_P)) ')'], 'FontSize', 7);
    end
    ylabel('MVDR Mag.', 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Freq.', 'FontSize', 12, 'FontName', 'Arial');
    title_name = sprintf('MVDR in MC=%d',i_MC);
    title(title_name, 'FontSize', 10)
end
[pks,locs] = findpeaks(mean(S_MVDR,2));
[sort_pks, idx] = sort(pks,'descend');
sort_locs = locs(idx);
subplot(2,2,4)
plot(f_grid, mean(S_MVDR,2), 'linewidth', 2);
set(gca, 'YScale', 'log')
hold on;
plot(f_grid(sort_locs(1:P)),sort_pks(1:P), 'ro', 'MarkerSize', 8, 'LineWidth', 2.5);
for i_P = 1 : P
    text(f_grid(sort_locs(i_P))+0.03,sort_pks(i_P),['(' num2str(f_grid(sort_locs(i_P))) ',' num2str(sort_pks(i_P)) ')'], 'FontSize', 7);
end
ylabel('MVDR Mag.', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Freq.', 'FontSize', 12, 'FontName', 'Arial');
title('Est. mean of MVDR', 'FontSize', 10)
saveas(gcf, 'figure_jpg/HW2_7d.jpg')
saveas(gcf, 'figure_fig/HW2_7d.fig')

%% Problem e
L = N-M+1;
f_grid = -0.5:0.001:0.5;
P_MUSIC = zeros(numel(f_grid),n_MC);
v = @(f) (exp(1j*2*pi*f*[0:M-1].'));
for i_MC = 1 : n_MC
    clc;
    fprintf('Calculate MUSIC (%d/%d)',i_MC,n_MC)
    [V,S,U] = svd(R_hat(:,:,L,i_MC));
    Un = V(:,P+1:M);
    for i_f = 1 : numel(f_grid)
        P_MUSIC(i_f, i_MC) = 1 / real(v(f_grid(i_f))'*Un*Un'*v(f_grid(i_f)));
    end
end
figure;
for i_MC = 1 : 3
    [pks,locs] = findpeaks(P_MUSIC(:,i_MC));
    [sort_pks, idx] = sort(pks,'descend');
    sort_locs = locs(idx);
    subplot(2,2,i_MC)
    plot(f_grid, P_MUSIC(:,i_MC), 'linewidth', 2);
    set(gca, 'YScale', 'log')
    hold on;
    plot(f_grid(sort_locs(1:P)),sort_pks(1:P), 'ro', 'MarkerSize', 8, 'LineWidth', 2.5);
    for i_P = 1 : P
        text(f_grid(sort_locs(i_P))+0.03,sort_pks(i_P),['(' num2str(f_grid(sort_locs(i_P))) ',' num2str(sort_pks(i_P)) ')'], 'FontSize', 7);
    end
    ylabel('MUSIC Mag.', 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Freq.', 'FontSize', 12, 'FontName', 'Arial');
    title_name = sprintf('MUSIC in MC=%d',i_MC);
    title(title_name, 'FontSize', 10)
end
[pks,locs] = findpeaks(mean(P_MUSIC,2));
[sort_pks, idx] = sort(pks,'descend');
sort_locs = locs(idx);
subplot(2,2,4)
plot(f_grid, mean(P_MUSIC,2), 'linewidth', 2);
set(gca, 'YScale', 'log')
hold on;
plot(f_grid(sort_locs(1:P)),sort_pks(1:P), 'ro', 'MarkerSize', 8, 'LineWidth', 2.5);
for i_P = 1 : P
    text(f_grid(sort_locs(i_P))+0.03,sort_pks(i_P),['(' num2str(f_grid(sort_locs(i_P))) ',' num2str(sort_pks(i_P)) ')'], 'FontSize', 7);
end
ylabel('MUSIC Mag.', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Freq.', 'FontSize', 12, 'FontName', 'Arial');
title('Est. mean of MUSIC', 'FontSize', 10)
saveas(gcf, 'figure_jpg/HW2_7e.jpg')
saveas(gcf, 'figure_fig/HW2_7e.fig')