function [dnW,dsW,W] = BONE_MP_kernelmatV2_dbl_periodic(...
              geom_stuff,MP_stuff,M)

soL   = geom_stuff{1};
xx    = geom_stuff{2};
yy    = geom_stuff{3};
LC    = geom_stuff{4};
theta = geom_stuff{5};
%%
centre   = MP_stuff{1};
AB       = MP_stuff{2};
%%
[W_x,W_y]   = BONE_diff_multipoles_laplace_doubly_periodic(...
             xx,yy,M,AB,centre);
dnW   = diag(sin(theta))*W_x - diag(cos(theta))*W_y;
dsW   = diag(cos(theta))*W_x + diag(sin(theta))*W_y;
dnW   = [real(dnW),imag(dnW)];
dsW   = [real(dsW),imag(dsW)];
%%
if nargout==3
  W   = BONE_multipoles_laplace_doubly_periodic(...
             xx,yy,M,AB,centre);
  W   = [real(W),imag(W)];
end
