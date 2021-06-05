clear; clc; close all;

eGen    = 5000;
k       = 2;
hT      = 70;
hB      = 10;
Tinf    = 25;
delta   = 0.05 ;
qf      = 50;
Tright  = 45;

Lx = 1.5;
Ly = 1;

Ni = Lx/delta + 1; 
Nj = Ly/delta + 1;

N = Ni*(Nj-1); % Because of right boundary

A = zeros(N,N);
b = zeros(N,1);

id = @(i,j) (j-1)*Ni+i;

row = 1;
% Internal nodes
for j = 1:Nj-1
    for i = 1:Ni
 
        LEFT    = j-1;
        RIGHT   = j+1;
        BOTTOM  = i+1;
        TOP     = i-1;

        %% Internal Nodes
        if ((LEFT ~=0) && (RIGHT~=Ni) && (TOP~=Nj+1) && (BOTTOM~=0))        
            A(row,id(i,LEFT))   =  1;
            A(row,id(i,RIGHT))  =  1; 
            A(row,id(BOTTOM,j)) =  1;  
            A(row,id(TOP,j))    =  1; 
            A(row,id(i,j))      = -4;   % Center
            % RHS vector
            b(row) = -eGen*delta^2/k;
        end
         
        %% Boundary Nodes
        % Left
        if (LEFT == 0 && (TOP ~=0) && (BOTTOM ~= Ni+1))
            A(row,id(i,RIGHT))  =  2;  
            A(row,id(BOTTOM,j)) =  1;  
            A(row,id(TOP,j))    =  1; 
            A(row,id(i,j))      = -4;   % Center
            b(row) = -eGen*delta^2/k - 2*qf*delta/k;
        end
         
        % Right
        if (RIGHT == Nj && (TOP ~= 0) && (BOTTOM ~= Ni+1))
            A(row,id(i,LEFT))   =  1;
            A(row,id(BOTTOM,j)) =  1;  
            A(row,id(TOP,j))    =  1; 
            A(row,id(i,j))      = -4;   % Center
            b(row) = -eGen*delta^2/k -Tright;
        end
        
         
        % Top
        if ((TOP == 0) && LEFT ~=0 && RIGHT ~=Nj)
            A(row,id(i,LEFT))   =  1;
            A(row,id(i,RIGHT))  =  1;  
            A(row,id(BOTTOM,j)) =  2;
            A(row,id(i,j)) = -4 - 2*delta*hT/k;
            b(row) = -eGen*delta^2/k - 2*delta*hT/k*Tinf;
        end
         
        % Bottom
        if (BOTTOM == Ni && LEFT ~=0 && RIGHT ~=Nj)
            A(row,id(i,LEFT))   =  1;
            A(row,id(i,RIGHT))  =  1;   
            A(row,id(TOP,j))    =  2; 
            A(row,id(i,j)) = -4 - 2*delta*hB/k;
            b(row) = -eGen*delta^2/k - 2*delta*hB/k*Tinf;
        end
         
        %% Corner Nodes
        if ((TOP == 0) && LEFT == 0) % Top-Left
            A(row,id(i,RIGHT))  =  2;  
            A(row,id(BOTTOM,j)) =  2;   
            A(row,id(i,j))      = -4 - 2*delta*hT/k;   % Center

            % RHS vector
            b(row) = -eGen*delta^2/k - 2*qf*delta/k - 2*delta*hT/k*Tinf;
        end
         
        if ((TOP == 0) && RIGHT == Nj) % Top-Right
            A(row,id(i,LEFT))   =  1;
            A(row,id(BOTTOM,j)) =  2;  
            A(row,id(i,j))      = -4-2*delta*hT/k;   % Center
            
            b(row) = -eGen*delta^2/k - 2*delta*hT/k*Tinf - Tright;
        end
         
        if (BOTTOM == Ni && LEFT == 0) % Bottom-Left
            A(row,id(i,RIGHT))  =  2;   
            A(row,id(TOP,j))    =  2; 
            A(row,id(i,j))      = -4 - 2*delta*hB/k;   % Center
            
            % RHS vector
            b(row) = -eGen*delta^2/k - 2*qf*delta/k - 2*delta*hB/k*Tinf;
        end
              
        if (BOTTOM == Ni && RIGHT == Nj) % Bottom-Right
            A(row,id(i,LEFT)) =  1;
            A(row,id(TOP,j))  =  2;  
            A(row,id(i,j))    = -4 - 2*delta*hB/k;   % Center
            
            b(row) = -eGen*delta^2/k - 2*delta*hB/k*Tinf - Tright;
        end
               
        row = row + 1; 
    end
end
%size(A)
T = A\b;
T = reshape(T,Ni,Nj);
contourf(T',20,'edgecolor','none');
colormap(jet)
colorbar
axis equal