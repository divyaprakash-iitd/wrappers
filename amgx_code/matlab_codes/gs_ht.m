clear; clc; close all;

N = 10;          % Number of unknowns
T = zeros(1,N);  % Initial guess
B = zeros(1,N);  % RHS vector
R = ones(1,N);   % Residual vector

% Boundary Conditions
TL = 300;
TR = 100;

NN      = 3;
LEFT    = 1;
CENTER  = 2;
RIGHT   = 3;

A = zeros(N,NN);

A(2:N-1,LEFT)   = 1;
A(2:N-1,RIGHT)  = 1;
A(:,CENTER)     = -2;

A(1,LEFT) = 0; A(1,RIGHT) = 1;
A(N,LEFT) = 1; A(N,RIGHT) = 0;

B(1) = -TL; B(N) = -TR;

MAXITER = 1000;
ITER = 0;
ERR = 1;
TOLERANCE = 1E-10;

figure(1)
hold on
xlabel('Iterations')
ylabel('||R||')
while and((ITER < MAXITER),(ERR > TOLERANCE))
    for i = 1:N
        if (i == 1)
            TTemp = B(i) - A(i,LEFT)*TL - A(i,RIGHT)*T(i+1);
        elseif (i == N)
            TTemp = B(i) - A(i,LEFT)*T(i-1) - A(i,RIGHT)*TR;
        else
            TTemp = B(i) - A(i,LEFT)*T(i-1) - A(i,RIGHT)*T(i+1);
        end
        T(i) = 1/A(i,CENTER)*TTemp;
    end
    
    % Calculate Residual: R = AT - B
    for i = 1:N
        if (i == 1)
            R(i) = A(i,LEFT)*TL + A(i,CENTER)*T(i) + A(i,RIGHT)*T(i+1) - B(i);
        elseif (i == N)
            R(i) = A(i,LEFT)*T(i-1) + A(i,CENTER)*T(i) + A(i,RIGHT)*TR - B(i);
        else
            R(i) = A(i,LEFT)*T(i-1) + A(i,CENTER)*T(i) + A(i,RIGHT)*T(i+1) - B(i);
        end
    end
    ITER = ITER + 1;
    ERR = norm(R(:));
    plot(ITER,ERR,'rx')
end
set(gca,'yscale','log')