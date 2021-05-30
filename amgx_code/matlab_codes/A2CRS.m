% Generates the coefficient matrix of the laplacian operator in the
% pressure poisson equation
%
% Author: Divyaprakash (2020AMZ8461), PhD Candidate
%         Applied Mechanics Department, IIT Delhi
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 06 April 2021

clear; clc; close all;

global Ni Nj Nk
Ni = 3; Nj = 3; Nk = 3;
dx = 0.1; dz = 0.1;

% Stretched grid in y-direction (% Placeholder values)
deltayp = 0.1*ones(1,Nj);
deltayv = 0.1*ones(1,Nj);

invdeltax2 = 1/dx/dx;
invdeltaz2 = 1/dz/dz;

N = Ni*Nj*Nk; % Matrix size
row = 1;

% Coefficient Matrix
C = zeros(N,N);

% Write data to a file
fileID = fopen('mesh.txt','w');

% CRS Format data
%nz = 6*Ni*Nj*Nk + 7*(Ni-2)*(Nj-2)*(Nk-2); % Number of non-zero elements (7 per row)
nz = 171;
val = zeros(1,nz); % Non-zero values of the coefficient matrix
col_ind = zeros(1,nz); % Column-wise locations of each element
row_ptr = zeros(1,N+1); % A value for each row
row_ptr(1,N+1) = nz +1; 
NNZ = 1; % Initialze the NNZ count

for k =1:Nk
	for j = 1:Nj
		for i = 1:Ni 
			% x-direction
			C(row, idx(i-1,j,k)) = invdeltax2;				 	% Left Node
			C(row, idx(i+1,j,k)) = invdeltax2; 					% Right Node

            % CRS Data
            
			% z-direction
			C(row, idx(i,j,k-1)) = invdeltaz2;				 	% Back Node
			C(row, idx(i,j,k+1)) = invdeltaz2;				 	% Front Node
            
            % y-direction
            if j == 1
                C(row, idx(i,j+1,k)) = 1/deltayv(j+1)/deltayp(j); % Top Node
                
                % Center Node  
                C(row, idx(i,j,k))   = -2.1*invdeltax2 - 2.1*invdeltaz2 -1/deltayv(j+1)/deltayp(j);
            elseif j == Nj
                C(row, idx(i,j-1,k)) = 1/deltayv(j)/deltayp(j); 	% Bottom Node

                % Center Node
                C(row, idx(i,j,k))   = -2.1*invdeltax2 - 2.1*invdeltaz2 -1/deltayv(j)/deltayp(j);
            else
                C(row, idx(i,j-1,k)) = 1/deltayv(j)/deltayp(j); 	% Bottom Node
                C(row, idx(i,j+1,k)) = 1/deltayv(j+1)/deltayp(j); 	% Top Node
                
                % Center Node
                C(row, idx(i,j,k))   = -2.1*invdeltax2 - 2.1*invdeltaz2 -1/deltayp(j)*(1/deltayv(j+1) + 1/deltayv(j));
            end
            
            fprintf(fileID,'%d, %d, %d: %d\n',i,j,k,idx(i,j,k));
            
            % Inspect the row for non-zero values
            chk = 0; % Initialize before each row
            for m = 1:N
                if (abs(C(row,m)) > 1e-14)
                    if chk == 0
                        row_ptr(row) = NNZ;
                        chk = 1;
                    end
                    
                    val(NNZ) = C(row,m);
                    col_ind(NNZ) = m;
                    NNZ = NNZ + 1;
                end
            end
            
			row = row + 1;
		end
	end
end

%spy(C)

function id = idx(i,j,k)
% Transforms 3D location to 1D

    global Ni Nj Nk
    
    % Periodicity in x
    i = (i < 1)*Ni + (i>Ni)*1 + (i<=Ni && i>=1)*i;
    
    % Periodicity in z
    k = (k < 1)*Nk + (k>Nk)*1 + (k<=Nk && k>=1)*k;
    
    id = (k-1)*(Ni*Nj) + (j-1)*Ni + i;
end
