function [dF_dx,dF_dy,F]=...
  BONE_diff_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre)

F        = BONE_multipoles_laplace_doubly_periodic(...
             X,Y,M+1,AB,centre);
AB2      = [AB(2),AB(1)];
centre2  = [centre(2),-centre(1)];
Fiz      = BONE_multipoles_laplace_doubly_periodic(...
             Y,-X,2,AB2,centre2);
%%
Z     = X-centre(1)+1i*(Y-centre(2));
np    = length(Z);
dF_dx = zeros(np,M);
dF_dy = dF_dx;
%%
A     = AB(1);
B     = AB(2);
alp   = pi/(2*A);
bet   = pi/(2*B);
Zm    = -1i*Z;
%%
for j=2:2:M
  dF_dx(:,j)   = -j*F(:,j+1);
end
for j=3:2:M
  dF_dx(:,j)   = -j*F(:,j+1)+alp^2*(j-1)*F(:,j-1);
end
dF_dy = 1i*dF_dx;
%%
dF_dx(:,1)  = -real(F(:,2))-1i*real(-1i*Fiz(:,2));
dF_dy(:,1)  = -real(1i*F(:,2))-1i*real(Fiz(:,2));
%%
F(:,M+1) = [];
