function Y=...
  GRN_laplace_hexag_periodic(dX,dY,a,b,th_vec,USE_SING)
%% periodic Green's function for Laplace's eqn,
%%  using hexagonal cells:
%%  _
%% / \   of height 2a and width 2b
%% \_/

if isempty(a)
  a=b*sqrt(3)/2;
end

Y = GRN_laplace_doubly_periodic(dX,dY,...
      2*a,3*b,th_vec,USE_SING) + ...
  + GRN_laplace_doubly_periodic(dX+1.5*a,dY+b,...
      2*a,3*b,th_vec,USE_SING);