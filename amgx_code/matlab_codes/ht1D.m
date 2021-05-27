clear; clc; close all;

% Boundary conditions
TL = 300;
TR = 100;

n = 1e3; % Number of grid points, including boundary nodes
N = n-2; % Number of unknowns

A1 = diag(ones(1,N-1),1);
A2 = diag(-2*ones(1,N));
A3 = diag(ones(1,N-1),-1);

% Coefficient Matrix
A = A1 + A2 + A3;

b = zeros(N,1);
b(1) = -TL; b(end) = -TR;

T = A\b;







