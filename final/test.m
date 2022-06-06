clear, clc;

n_MC = 1000;                           % number of monte carlo test
N = 100;

sigma_err = @(sigma_hat,sigma_true) (max([ round( 10 * ( 1 - norm(sigma_hat(:) - sigma_true(:)) / norm(sigma_true) ) ), 0]));
f_err = @(f_hat,f_true) (max([ round( 10 * ( 1 - norm(f_hat(:) - f_true(:)) / norm(f_true) ) ), 0]));
sigma_w_err = @(sigma_w_hat,sigma_w_true) (max([ round( 5 * ( 1 - norm(sigma_w_hat(:) - sigma_w_true(:)) / norm(sigma_w_true) ) ), 0]));

sigma_score1 = zeros(n_MC,1);
sigma_score2 = zeros(n_MC,1);
sigma_score3 = zeros(n_MC,1);
sigma_score4 = zeros(n_MC,1);
f_score1 = zeros(n_MC,1);
f_score2 = zeros(n_MC,1);
f_score3 = zeros(n_MC,1);
f_score4 = zeros(n_MC,1);
sigma_w_score1 = zeros(n_MC,1);
sigma_w_score2 = zeros(n_MC,1);
sigma_w_score3 = zeros(n_MC,1);
sigma_w_score4 = zeros(n_MC,1);
tic;
% for i_MC = 1 : n_MC
parfor i_MC = 1 : n_MC
    disp(['Running MC = ', num2str(i_MC)])
    f_true = sort(rand([1 5]) - 0.5);    % [-0.5, 0.5]
    sigma_true = rand([1 5]) * 1 + 0.5;  % [50, 150]
    sigma_w_true = rand()*0.5 + 1;       % [0.5, 1.5]
    noise = (randn(1,N) + 1j*randn(1,N)) / sqrt(2);
    x_tilde = sum(sigma_true.'.*exp(1j*2*pi*f_true.'*[0:N-1]),1)*exp(1j*(pi/6)) + sigma_w_true*noise; % 1 x N
    
    [f_hat, sigma_hat, sigma_w_hat] = OMP(x_tilde*exp(-1j*pi/6));
    sigma_score1(i_MC) = sigma_err(sigma_hat, sigma_true);
    f_score1(i_MC) = f_err(f_hat, f_true);
    sigma_w_score1(i_MC) = sigma_w_err(sigma_w_hat, sigma_w_true);
    
    [f_hat, sigma_hat, sigma_w_hat] = CoSaMP(x_tilde*exp(-1j*pi/6));
    sigma_score2(i_MC) = sigma_err(sigma_hat, sigma_true);
    f_score2(i_MC) = f_err(f_hat, f_true);
    sigma_w_score2(i_MC) = sigma_w_err(sigma_w_hat, sigma_w_true);

    [f_hat, sigma_hat, sigma_w_hat] = IHT(x_tilde*exp(-1j*pi/6));
    sigma_score3(i_MC) = sigma_err(sigma_hat, sigma_true);
    f_score3(i_MC) = f_err(f_hat, f_true);
    sigma_w_score3(i_MC) = sigma_w_err(sigma_w_hat, sigma_w_true);

    [f_hat, sigma_hat, sigma_w_hat] = HTP(x_tilde*exp(-1j*pi/6));
    sigma_score4(i_MC) = sigma_err(sigma_hat, sigma_true);
    f_score4(i_MC) = f_err(f_hat, f_true);
    sigma_w_score4(i_MC) = sigma_w_err(sigma_w_hat, sigma_w_true);
end
disp(['Number of Monte Carlo Sim. = ', num2str(n_MC)])
fprintf("| Algorithm     | %6s | %6s | %6s | %6s |\n"," OMP ", "CoSaMP", " IHT ", " HTP ")
fprintf("| Sigma Score   | %.4f | %.4f | %.4f | %.4f | // Max score = 10\n",mean(sigma_score1), mean(sigma_score2), mean(sigma_score3), mean(sigma_score4))
fprintf("| F Score       | %.4f | %.4f | %.4f | %.4f | // Max score = 10\n",mean(f_score1), mean(f_score2), mean(f_score3), mean(f_score4))
fprintf("| Sigma_w Score | %.4f | %.4f | %.4f | %.4f | // Max score = 5\n",mean(sigma_w_score1), mean(sigma_w_score2), mean(sigma_w_score3), mean(sigma_w_score4))
toc
