clear, clc;
load("MSA_Final.mat")

[f_hat, sigma_hat, sigma_w_hat] = OMP(x_tilde*exp(-1j*pi/6));
save('MSA_Final_Results.mat')
[f_hat, sigma_hat, sigma_w_hat] = CoSaMP(x_tilde*exp(-1j*pi/6));
[f_hat, sigma_hat, sigma_w_hat] = IHT(x_tilde*exp(-1j*pi/6));
[f_hat, sigma_hat, sigma_w_hat] = HTP(x_tilde*exp(-1j*pi/6));
