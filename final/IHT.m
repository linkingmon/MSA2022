function [f_hat, sigma_hat, sigma_w_hat] = IHT(x_tilde)
    % x_tilde [1 x N]
    f_hat = [];        % 1 x P
    sigma_hat = [];    % 1 x P
    sigma_w_hat = [];  % 1 x 1

    % system parameters
    plot_on = true;
    plot_on = false;
    G = 200; % number of grids
    N = 100;  % number of samples
    s = 5;    % number of sparse level

    % initialize parameters
    DFTN = dftmtx(N);
    DFT_grid = linspace(0,1,N);
    DFT_grid(DFT_grid > 0.5) = DFT_grid(DFT_grid > 0.5) - 1;
    [~,DFT_seq] = sort(DFT_grid);
    f_grid = linspace(-0.5,0.5,G);     % 1 x G
    A = exp(j*2*pi*f_grid.*[0:N-1].'); % N x G

    if plot_on
        figure
        subplot(611)
        x_residual_freq = abs(DFTN*x_tilde.');
        plot(DFT_grid(DFT_seq),x_residual_freq(DFT_seq))
    end

    % run OMP
    f_sparse = zeros(G,1);
    sup_f = [];
    for i_P = 1 : 5
        f_new = f_sparse + A'*(x_tilde.'-A*f_sparse)/G;
        % f_new = f_sparse + pinv(A)*(x_tilde.'-A*f_sparse);
        [~, sup_f] = maxk(real(f_new),s); % sup_f: s x 1
        f_sparse = zeros(G,1);
        f_sparse(sup_f) = real(f_new(sup_f));
        f_hat = f_grid(sup_f);
        sigma_hat = real(f_new(sup_f));
        if plot_on
            subplot(6,1,i_P+1)
            x_residual = x_tilde.' - A(:,sup_f)*sigma_hat;
            x_residual_freq = abs(DFTN*x_residual);
            plot(DFT_grid(DFT_seq),x_residual_freq(DFT_seq))
        end
    end
    x_residual = x_tilde.' - A(:,sup_f)*sigma_hat;
    sigma_w_hat = sqrt(x_residual'*x_residual/N);
    % Output Reorder
    sigma_hat = sigma_hat.';
    [~, f_seq] = sort(f_hat);
    f_hat = f_hat(f_seq);
    sigma_hat = sigma_hat(f_seq);
end