function [M_dnW,M_W,xtra]  = BONE_MP_kernelmat_dbl_periodic(...
              geom_stuff,ip_stuff,MP_stuff,M)

soL   = geom_stuff{1};
xx    = geom_stuff{2};
yy    = geom_stuff{3};
LC    = geom_stuff{4};
Xtra  = geom_stuff{5};
theta = Xtra{1};
%%
ip_d1chi = ip_stuff{1}{1}(2:end,:);
centre   = MP_stuff{1};
AB       = MP_stuff{2};
%%
W     = BONE_multipoles_laplace_doubly_periodic(...
             xx,yy,M,AB,centre);
M_W   = ip_d1chi*[real(W),imag(W)];
%%
[W_x,W_y]   = BONE_diff_multipoles_laplace_doubly_periodic(...
             xx,yy,M,AB,centre);
dnW   = diag(sin(theta))*W_x - diag(cos(theta))*W_y;
M_dnW = ip_d1chi*[real(dnW),imag(dnW)];
%%
if nargout==3
  dsW    = diag(cos(theta))*W_x + diag(sin(theta))*W_y;
  xtra   = {[real(W),imag(W)],[real(dnW),imag(dnW)],...
             [real(dsW),imag(dsW)]};
end
