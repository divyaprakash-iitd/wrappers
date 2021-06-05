clear; clc; close all;

N = flip(logspace(1,6,6));

GS = [69.182 6.934 0.782 0.076 0.007 0.001];
AMGX = [9.942 4.085 3.586 3.610 3.474 3.516];

plot(N,GS,'-o','DisplayName','GS')
hold on
plot(N,AMGX,'r-*','DisplayName','AMGX')
set(gca,'xscale','log')
xlabel('N')
ylabel('Time [s]')
legend('location','best')