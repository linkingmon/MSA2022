clear, clc;

load('MSA_HW1_Problem_7.mat'); % load V
[R, L] = size(V);
N = L/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = nan(R, L); % x(-L/2:L/2-1)

% x(n) = 0.5 x(n-1) - x(n-2) + 0.25 x(n-3) - 0.25 x(n-4) + v(n)
for i = 1 : L
    res = V(:,i) ;
    if(i > 1) res = res + 0.5*X(:,i-1); end
    if(i > 2) res = res - X(:,i-2); end
    if(i > 3) res = res + 0.25*X(:,i-3); end
    if(i > 4) res = res - 0.25*X(:,i-4); end
    X(:,i) = res;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Problem 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_grid = 0:0.001:1;
PSD_true = @(f) (1 ./ [3/2-3/2*cos(2*pi*f)+cos(4*pi*f)] ./ [5/4+cos(4*pi*f)] );

PSD_per = nan(R,N);
for r = 1 : R
    for k = 0 : N-1
        PSD_per(r,k+1) = abs(sum( X(r,[L/2:L/2+N-1]) .* exp(-j*2*pi*k/N*[0:N-1]) ))^2 / N;
    end
end

figure;
hold on;
h1 = plot([0:N-1]/N, PSD_per(1,:), 'linewidth', 2);
h2 = plot([0:N-1]/N, PSD_per(2,:), 'linewidth', 2);
h3 = plot([0:N-1]/N, PSD_per(3,:), 'linewidth', 2);
h4 = plot(f_grid, PSD_true(f_grid), 'linewidth', 2);
set(gca, 'YScale', 'log')
xlabel('Freq.', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Mag.', 'FontSize', 16, 'FontName', 'Arial');
legend([h1 h2 h3 h4], {'Periodogram in MC_1','Periodogram in MC_2','Periodogram in MC_3', 'True PSD'})
title('Monte Carlo for Periodogram', 'FontSize', 16)
grid on;
saveas(gcf, 'gen_fig/HW1_7a.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Problem 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSD_per_avg = mean(PSD_per,1);
W_B = @(f) diric(2*pi*f,N).^2*N;

PSD_per_true = nan(1,length(f_grid));
for i_f = 1 : length(f_grid)
    f = f_grid(i_f);
    PSD_true_conv_W_B = @(mu) (PSD_true(f-mu).*W_B(mu));
    PSD_per_true(i_f) = integral(PSD_true_conv_W_B,-0.5,0.5);
end 

figure;
hold on;
h3 = plot([0:N-1]/N, PSD_per_avg, 'linewidth', 2);
h1 = plot(f_grid, PSD_true(f_grid), 'linewidth', 2);
h2 = plot(f_grid, PSD_per_true, 'linewidth', 2);
set(gca, 'YScale', 'log')
xlabel('Freq.', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Mag.', 'FontSize', 16, 'FontName', 'Arial');
legend([h1 h2 h3], {'True PSD','True mean periodogram','Estimated mean periodogram'})
title('Estimated Periodogram', 'FontSize', 16)
grid on;
saveas(gcf, 'gen_fig/HW1_7b.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Problem 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSD_square = PSD_true(f_grid).^2;
MC_Var = sum((PSD_per-PSD_per_avg).^2,1)/(R-1);

figure;
hold on;
h2 = plot([0:N-1]/N, MC_Var, 'linewidth', 2);
h1 = plot(f_grid, PSD_square, 'linewidth', 2);
set(gca, 'YScale', 'log')
xlabel('Freq.', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Mag.', 'FontSize', 16, 'FontName', 'Arial');
legend([h1 h2], {'PSD square','Estimated periodogram var.'})
title('Estimated Periodogram', 'FontSize', 16)
grid on;
saveas(gcf, 'gen_fig/HW1_7c.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Problem 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_a = @(k) (1-abs(k)/N);
E_w = sum(w_a([-N+1:N-1]).^2);

r_x_hat = nan(R,2*N-1);
for r = 1 : R
    for k = -N+1 : 0
        r_x_hat(r,k+N) = mean(X(r,1:N).*conj(X(r,1-k:N-k)));
    end
    for k = 1 : N-1
        r_x_hat(r,k+N) = mean(X(r,L-N+1:L).*conj(X(r,L-N+1-k:L-k)));
    end
end

PSD_PS_true = PSD_true(f_grid).^2*E_w/N;

PSD_PS = nan(R,length(f_grid));
for r = 1 : R
    for i_f = 1 : length(f_grid)
        clc;
        fprintf("Calculating Freq. Grid = %4d/%d Under MC = %d/%d\n",i_f,length(f_grid),r,R)
        f = f_grid(i_f);
        sum_val = 0;
        for k = -N+1 : N-1
            sum_val = sum_val + r_x_hat(r,k+N) .* w_a(k) .* exp(-j*2*pi*f*k);
        end
        PSD_PS(r,i_f) = sum_val;
    end
end

PSD_PS = abs(PSD_PS);
PSD_PS_avg = mean(PSD_PS,1);
PS_Var = sum((PSD_PS-PSD_PS_avg).^2,1)/(R-1);

figure;
hold on;
h1 = plot(f_grid, PSD_PS_true, 'linewidth', 2);
h2 = plot(f_grid, PS_Var, 'linewidth', 2);
set(gca, 'YScale', 'log')
xlabel('Freq.', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Mag.', 'FontSize', 16, 'FontName', 'Arial');
legend([h1 h2], {'True PS var.','Estimated PS var.'})
title('Variance of Periodogram Smoothing', 'FontSize', 16)
grid on;
saveas(gcf, 'gen_fig/HW1_7d.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Problem 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 20;
seg_len = N/K;
PSD_square = PSD_true(f_grid).^2/K;
PSD_PA = nan(R,seg_len);
window = seg_len/2-abs([0:seg_len-1]-seg_len/2); window = window * sqrt(seg_len/sum(window.^2)); % Normalize window power to one
% window = ones(1,seg_len); window = window * sqrt(seg_len/sum(window.^2));
for r = 1 : R
    PSD_per = zeros(1,seg_len);
    for k = 1 : K
        idx_range = (k-1)*seg_len+L/2:k*seg_len+L/2-1;
        for l = 0 : seg_len-1 % f : [0,seg_len-1]/seg_len
            PSD_per(1,l+1) = PSD_per(1,l+1) + ...
                abs(sum( X(r,idx_range) .* window .* exp(-j*2*pi*l/seg_len*[0:seg_len-1]) ))^2 / seg_len;
        end
    end
    PSD_PA(r,:) = PSD_per / K;
end
PSD_PA_avg = mean(PSD_PA,1);
MC_Var = sum((PSD_PA-PSD_PA_avg).^2,1)/(R-1);

figure;
hold on;
h2 = plot([0:seg_len-1]/seg_len, MC_Var, 'linewidth', 2);
h1 = plot(f_grid, PSD_square, 'linewidth', 2);
set(gca, 'YScale', 'log')
xlabel('Freq.', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Mag.', 'FontSize', 16, 'FontName', 'Arial');
legend([h1 h2], {'True PA var.','Estimated PA var.'})
title('Variance of Periodogram Averaging', 'FontSize', 16)
grid on;
saveas(gcf, 'gen_fig/HW1_7e.jpg')
