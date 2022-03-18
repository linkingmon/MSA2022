clear, clc;
p1 = (1+sqrt(7)*1j) / 4;
p2 = (1-sqrt(7)*1j) / 4;
p3 = -1j / sqrt(2);
p4 = 1j / sqrt(2);

hold on;
xlim([-2 2])
ylim([-2 2])

max_r = max(abs(p1),abs(p3));
rectangle('Position',[-3 -3 6 6], 'LineWidth', 1.5,'FaceColor',[.5 .5 .5]);
rectangle('Position',[-1 -1 2 2]*max_r,'Curvature',[1,1], 'LineWidth', 1.5,'FaceColor',[1 1 1],'LineStyle','--');
plot(p1,'rx', 'MarkerSize', 10, 'LineWidth', 3)
plot(p2,'rx', 'MarkerSize', 10, 'LineWidth', 3)
plot(p3,'rx', 'MarkerSize', 10, 'LineWidth', 3)
plot(p4,'rx', 'MarkerSize', 10, 'LineWidth', 3)

pbaspect([1 1 1])
xlabel('Re', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Im', 'FontSize', 16, 'FontName', 'Arial');
saveas(gcf, 'gen_fig/HW1_6a.jpg')
