clear; clc; close all;

n = 7;
N = n-2; % Number of unknowns

% Number of nnz
nnz = 2*2 + 3*(N-2);

val = zeros(1,nnz);
col_ind = zeros(1,nnz);

id = 1;
for row = 1:N    
    if (row == 1)
        val(id) = -2; col_ind(id) = 0; id = id + 1;
        val(id) =  1; col_ind(id) = 1; id = id + 1;
    elseif (row == N)
        val(id) =  1; col_ind(id) = row - 3 + 1; id = id + 1;
        val(id) = -2; col_ind(id) = row - 3 + 2; id = id + 1;
    else
        val(id) =  1; col_ind(id) = row - 3 + 1; id = id + 1;
        val(id) = -2; col_ind(id) = row - 3 + 2; id = id + 1;
        val(id) =  1; col_ind(id) = row - 3 + 3; id = id + 1;
    end
end

row_ptr = [0, 2:3:nnz-2,nnz];
