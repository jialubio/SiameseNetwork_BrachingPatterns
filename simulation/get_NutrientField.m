% Get nutrient field by solving 2d diffusion equation by FED
% N0 is the defined geometry of nutrient field of size(nx,ny)
function N = get_NutrientField(xx,yy,Nmax, N, nt)
[nx,ny] = size(N);
dx = xx(1,2) - xx(1,1);
dy = dx;
Nmax = max(max(N));
L = max(max(xx))*2;
dt = 0.5;

Nn = zeros(nx, ny);
% Dirichlet boundary condition
NW = 0;                            %x=0 Dirichlet B.C
NE = 0;                            %x=L Dirichlet B.C
NS = 0;                            %y=0 Dirichlet B.C
NN = 0;                            %y=L Dirichlet B.C

vis = 5;  % diffusion coefficient

% Assign boundary condition
%B.C vector
bc = zeros(nx-2,ny-2);
bc(1,:) = NW/dx^2; bc(nx-2,:) = NE/dx^2;  %Dirichlet B.Cs
bc(:,1) = NS/dy^2; bc(:,ny-2) = NN/dy^2;  %Dirichlet B.Cs
%B.Cs at the corners:
bc(1,1)    = NW/dx^2+NS/dy^2; bc(nx-2,1)    = NE/dx^2+NS/dy^2;
bc(1,ny-2) = NW/dx^2+NN/dy^2; bc(nx-2,ny-2) = NE/dx^2+NN/dy^2;
bc = vis * dt * bc;

%Calculating the coefficient matrix for the implicit scheme
Ex = sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax = Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
Ey = sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay = Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
A  = kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D  = speye((nx-2)*(ny-2))-vis*dt*A;

for it = 0:nt,
    Nn = N;    
    
    %NNcomment as necessary
    %Implicit method:
    U = Nn;
    U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[]; % BC
    U = reshape(U+bc,[],1);
    U = D\U;
    U = reshape(U,nx-2,ny-2);
    N(2:nx-1,2:ny-1) = U;
    %BoNNdary conditions
    %Dirichlet:
    N(1,:)  = NW;
    N(nx,:) = NE;
    N(:,1)  = NS;
    N(:,ny) = NN;
end
N = N/max(max(N));
end