clear; clc; close all;

eGen    = 5000;
k       = 2;
hT      = 70;
hB      = 10;
Tinf    = 25;
delta   = 0.05;
qf      = 500;
Tright  = 45;

Lx = 1.5;
Ly = 1;

Ni = ceil(Lx/delta) + 1; 
Nj = ceil(Ly/delta) + 1;

N = (Ni-1)*(Nj); % Because of right boundary

A = sparse(zeros(N,N));
b = sparse(zeros(N,1));

id = @(i,j) (j-1)*(Ni-1)+i;

row = 1;
% Internal nodes
for j = 1:Nj
    for i = 1:Ni-1
        
        LEFT    = i-1;
        RIGHT   = i+1;
        BOTTOM  = j-1;
        TOP     = j+1;

        %% Internal Nodes
        if (LEFT ~=0 && RIGHT~=Ni && TOP~=Nj+1 && BOTTOM~=0)        
            A(row,id(LEFT,j))   =  1;
            A(row,id(RIGHT,j))  =  1;  
            A(row,id(i,BOTTOM)) =  1;  
            A(row,id(i,TOP))    =  1; 
            A(row,id(i,j))      = -4;   % Center

            % RHS vector
            b(row) = -eGen*delta^2/k;
        end
        
        %% Boundary Nodes
        % Left
        if (LEFT == 0 && TOP ~=Nj+1 && BOTTOM ~= 0)
            A(row,id(RIGHT,j))  =  2;  
            A(row,id(i,BOTTOM)) =  1;  
            A(row,id(i,TOP))    =  1; 
            A(row,id(i,j))      = -4;   % Center
            b(row) = -eGen*delta^2/k - 2*qf*delta/k;
        end
        
        % Right
        if (RIGHT == Ni && TOP ~= Nj+1 && BOTTOM ~= 0)
            A(row,id(LEFT,j))   =  1;
            A(row,id(i,BOTTOM)) =  1;  
            A(row,id(i,TOP))    =  1; 
            A(row,id(i,j))      = -4;   % Center
            b(row) = -eGen*delta^2/k -Tright;
        end
        
        
        % Top
        if (TOP == Nj + 1 && LEFT ~=0 && RIGHT ~=Ni)
            A(row,id(LEFT,j))   =  1;
            A(row,id(RIGHT,j))  =  1;  
            A(row,id(i,BOTTOM)) =  2;
            A(row,id(i,j)) = -4 - 2*delta*hT/k;
            b(row) = -eGen*delta^2/k - 2*delta*hT/k*Tinf;
        end
        
        % Bottom
        if (BOTTOM == 0 && LEFT ~=0 && RIGHT ~=Ni)
            A(row,id(LEFT,j))   =  1;
            A(row,id(RIGHT,j))  =  1;   
            A(row,id(i,TOP))    =  2; 
            A(row,id(i,j)) = -4 - 2*delta*hB/k;
            b(row) = -eGen*delta^2/k - 2*delta*hB/k*Tinf;
        end
        
        %% Corner Nodes
        if (TOP == Nj+1 && LEFT == 0) % Top-Left
            A(row,id(RIGHT,j))  =  2;  
            A(row,id(i,BOTTOM)) =  2;   
            A(row,id(i,j))      = -4 - 2*delta*hT/k;   % Center

            % RHS vector
            b(row) = -eGen*delta^2/k - 2*qf*delta/k - 2*delta*hT/k*Tinf;
        end
        
        if (TOP == Nj+1 && RIGHT == Ni) % Top-Right
            A(row,id(LEFT,j))   =  1;
            A(row,id(i,BOTTOM)) =  2;  
            A(row,id(i,j))      = -4-2*delta*hT/k;   % Center
            
            b(row) = -eGen*delta^2/k - 2*delta*hT/k*Tinf - Tright;
        end
        
        if (BOTTOM == 0 && LEFT == 0) % Bottom-Left
            A(row,id(RIGHT,j))  =  2;   
            A(row,id(i,TOP))    =  2; 
            A(row,id(i,j))      = -4 - 2*delta*hB/k;   % Center
            
            % RHS vector
            b(row) = -eGen*delta^2/k - 2*qf*delta/k - 2*delta*hB/k*Tinf;
        end
             
        if (BOTTOM == 0 && RIGHT == Ni) % Bottom-Right
            A(row,id(LEFT,j)) =  1;
            A(row,id(i,TOP))  =  2;  
            A(row,id(i,j))    = -4 - 2*delta*hB/k;   % Center
            
            b(row) = -eGen*delta^2/k - 2*delta*hB/k*Tinf - Tright;
        end
               
        row = row + 1; 
    end
end

T = A\b;
T = reshape(T,Ni-1,Nj);
T = T';
T(:,end+1) = 45;
contourf(T,40,'edgecolor','none');
colormap(jet)
colorbar
%contour(T,40,'showText','off');
axis equal