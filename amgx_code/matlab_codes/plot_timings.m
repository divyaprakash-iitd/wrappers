clear; clc; close all;

gs_dense = 93.744;
gs_sparse = [0.00195 0.01087 2.44072 32.62450];
gs_amgx = [0.00321104 0.04456 0.761726 6.70354];

N = [16*11 31*21 152*102 302*202];

hold on
plot(N(1),gs_dense,'-^','DisplayName','SOR-Dense')
plot(N,gs_sparse,'-o','DisplayName','SOR-Sparse')
plot(N,gs_amgx,'-*','DisplayName','AMGX-CRS')
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('location','best')
xlabel('No. of Grid Points')
ylabel('Execution Time [s]')
grid on
set(gcf,'color','w')
export_fig plot_timings.pdf