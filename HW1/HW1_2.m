clear, clc;

H = @(f,k,N) diric(2*pi*(f-k/N),N);

N = 16;
f = 0:0.001:1;
figure;
hold on;
for k = 0 : N-1
    plot(f,abs(H(f,k,N)), 'linewidth', 2)
end
xlabel('Freq.', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Mag.', 'FontSize', 16, 'FontName', 'Arial');
title('Filter Bank for Periodogram', 'FontSize', 16)
grid on
saveas(gcf, 'gen_fig/HW1_2.jpg')
