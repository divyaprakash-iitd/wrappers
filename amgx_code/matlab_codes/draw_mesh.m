clear; clc; close all;

delta   = 0.25;

Lx = 1.5;
Ly = 1;

Ni = Lx/delta + 1; 
Nj = Ly/delta + 1;

xp = linspace(0,Lx,Ni);
yp = linspace(0,Ly,Nj);

[x, y] = meshgrid(xp,yp);

z = -0.1 + 0*x;

hold on
mesh(x,y,z)
plot(x,y,'ko')
axis equal
axis off
view(2)
set(gcf,'color','w')
export_fig mesh_ht.pdf

