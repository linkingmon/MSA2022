function [f_hat, sigma_hat, sigma_w_hat] = OMP(x_tilde)
    % x_tilde [1 x N]
    f_hat = [];        % 1 x P
    sigma_hat = [];    % 1 x P
    sigma_w_hat = [];  % 1 x 1

    % system parameters
    plot_on = true;
    plot_on = false;
    G = 1000; % number of grids
    N = 100;  % number of samples

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
    sup_f = [];
    f_sparse = zeros(G,1);
    for i_P = 1 : 5
        amp = A'*(x_tilde.'-A*f_sparse);
        [~, idx] = max(abs(amp));
        sup_f = [sup_f idx];
        f_hat = [f_hat f_grid(idx)];
        sigma_hat = real(pinv(A(:,sup_f)) * x_tilde.'); % i_P x 1
        f_sparse(sup_f) = sigma_hat;
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