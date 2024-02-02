% clc; clear all; clf;

data = load('GnR_out_');
% load('orig.mat')
time = 0:1:length(data) - 1;
a = data(1:time(end) + 1,1);
h = data(1:time(end) + 1,2);
rhoR = data(1:time(end) + 1,3);
rhoR_p1 = data(1:time(end) + 1,4);
rhoR_p2 = data(1:time(end) + 1,5);
rhoR_i = data(1:time(end) + 1,6);
wss = data(1:time(end) + 1,7);
wss_h = data(1:time(end) + 1,8);
P = data(1:time(end) + 1,9);
P_h = data(1:time(end) + 1,10);
f = data(1:time(end) + 1,11);
f_h = data(1:time(end) + 1,12);
Q = data(1:time(end) + 1,13);
Q_h = data(1:time(end) + 1,14);
time = time;
% a = data(1:time(end) + 1,1)*1001;
% h = data(1:time(end) + 1,2)*1001;
% rhoR = data(1:time(end) + 1,3);
% rhoR_p1 = data(1:time(end) + 1,4);
% rhoR_p2 = data(1:time(end) + 1,5);
% rhoR_i = data(1:time(end) + 1,6);
% wss = data(1:time(end) + 1,7);
% wss_h = data(1:time(end) + 1,8);
% P = data(1:time(end) + 1,9);
% P_h = data(1:time(end) + 1,10);
% f = data(1:time(end) + 1,11);
% f_h = data(1:time(end) + 1,12);
% Q = data(1:time(end) + 1,13);
% Q_h = data(1:time(end) + 1,14);

figure(1)
subplot(3,2,1)
hold on
plot(time / 7, a / a(1))
% plot(time(end), a_e/a(1), 's')
xlabel('Time (weeks)'); ylabel('radius')
% axis([ 0 180 0.95 1.2 ])

subplot(3,2,2)
hold on
plot(time / 7, h / h(1))
% plot(time(end), h_e/h(1), 's')
xlabel('Time (weeks)'); ylabel('thickness')
% axis([ -50 560 0.8 2.0])

% subplot(4,3,3)
% hold on
% plot(time, rhoR)
% % plot(time(end), h_e/h(1), 's')
% xlabel('Time (days)'); ylabel('rhoR (-)')
% % axis([ -50 560 0.8 2.0])

subplot(3,2,3)
hold on
plot(time, rhoR_p1)
% plot(time(end), h_e/h(1), 's')
xlabel('Time (days)'); ylabel('rho_p1 (-)')
% axis([ -50 560 0.8 2.0])

subplot(3,2,4)
hold on
plot(time, rhoR_p2)
% plot(time(end), h_e/h(1), 's')
xlabel('Time (days)'); ylabel('rho_p2 (-)')
% axis([ -50 560 0.8 2.0])

% subplot(4,3,6)
% hold on
% plot(time, rhoR_i)
% % plot(time(end), h_e/h(1), 's')
% xlabel('Time (days)'); ylabel('rho_i (-)')
% % axis([ -50 560 0.8 2.0])

subplot(3,2,5)
hold on
plot(time, wss)
plot(time,wss_h,'k--')
% plot(time(end), h_e/h(1), 's')
xlabel('Time (days)'); ylabel('wss')
% legend('WSS','Reference WSS')
% axis([ -50 560 0.8 2.0])

subplot(3,2,6)
hold on
plot(time, rhoR_i)
% plot(time(end), h_e/h(1), 's')
xlabel('Time (days)'); ylabel('infl mat')
% axis([ -50 560 0.8 2.0])
% 
% subplot(4,3,9)
% hold on
% plot(time,f)
% % plot(time,f_h)
% % plot(time(end), h_e/h(1), 's')
% xlabel('Time (days)'); ylabel('F')
% % legend('F','Reference F')
% % axis([ -50 560 0.8 2.0])
% 
% subplot(4,3,10)
% hold on
% plot(time, Q)
% plot(time, Q_h)
% % plot(time(end), h_e/h(1), 's')
% xlabel('Time (days)'); ylabel('Q')
% % axis([ -50 560 0.8 2.0])

% subplot(3,3,3)
% hold on
% plot(time,ups_c/ups_c(1))
% plot(time(end), ups_c_e/ups_c(1), 's')
% xlabel('Time (days)'); ylabel('Ups_c (-)')
% % axis([ -50 560 0.95 1.3])
% 
%subplot(3,2,6)
%hold on
%plot(time,rhoR_m)
% plot(time(end), rho_m_e/rho_m(1), 's')
%xlabel('Time (days)'); ylabel('rhoR_m')
% % axis([-50 560 0.75 3.5])
% 
% subplot(3,3,5)
% hold on
% plot(time,rho_c/rho_c(1))
% plot(time(end), rho_c_e/rho_c(1), 's')
% xlabel('Time (days)'); ylabel('rhoR_c (-)')
% % axis([ -50 560 0.75 3.5])
% 
% subplot(3,3,6)
% hold on
% plot(time,perturb/perturb(1))
% plot(time(end), perturb_e/perturb(1), 's')
% xlabel('Time (days)'); ylabel('Perturb (-)')
% % axis([ -50 560 0.95 1.6])
