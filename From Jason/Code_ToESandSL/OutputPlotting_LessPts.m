close all; clear all; clf;
    
load general.dat
varsout(1,:) = general(:,2)*10^6;
varsout(2,:) = general(:,4)*10^6;
varsout(3,:) = general(:,13)/(450*10^-6);
varsout(4,:) = general(:,9)/(10^3);

vars_target = ones(length(general(:,1)),4);
vars_target(:,1) = 550;
vars_target(:,2) = 20;
vars_target(:,3) = 0.042;
vars_target(:,4) = 80;

time = general(:,1)/7;

figure(1)
subplot(2,2,1)
plot(time, varsout(1,:), 'LineWidth', 1.5, 'Color', 'k')
hold on
plot(time, vars_target(:,1), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
xlabel('Time (days)'); ylabel('Radius (um)')
axis([0 100 200 800])
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold',...
    'Xtick', 0:20:100, 'Ytick', 200:150:800)

subplot(2,2,2)
plot(time, varsout(2,:), 'LineWidth', 1.5, 'Color', 'k')
hold on
plot(time, vars_target(:,2), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
xlabel('Time (days)'); ylabel('Thickness (um)')
axis([0 100 0 300])
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold',...
    'Xtick', 0:20:100, 'Ytick', 0:75:300)

subplot(2,2,3)
plot(time, varsout(3,:), 'LineWidth', 1.5, 'Color', 'k')
hold on
plot(time, vars_target(:,3), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
xlabel('Time (days)'); ylabel('Compliance')
axis([0 100 0 0.1])
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold',...
    'Xtick', 0:20:100, 'Ytick', 0:0.02:0.1)

subplot(2,2,4)
semilogy(time, varsout(4,:), 'LineWidth', 1.5, 'Color', 'k')
hold on
plot(time, vars_target(:,4), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
xlabel('Time (days)'); ylabel('Stiffness (kPa)')
axis([0 100 1 1D4])
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold',...
    'Xtick', 0:20:100, 'Ytick', [1, 10, 100, 1000, 10000])
 