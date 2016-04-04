function Y=...
  GRN_laplace_singly_periodic(GF_args)
%% periodic Green's function for Laplace's eqn,
%%  with period of B in the y direction;
%% CALL: Y=GRN_laplace_singly_periodic(GF_args);
%%  where GF_args={dX,dY,B,USE_ABS,normal_deriv,USE_SING};
%%  - normal_deriv==[] => just get G;
%%      ELSE
%%    normal_derive==th_vec => get \pa_n.G;
%%  - USE_SING==0 => subtract off log(R)/2/pi;
%%      ELSE
%%    USE_SING==1 => retain log(R)/2/pi singularity;
%%  - USE_ABS==0 => leave out |dX|/4/B
%%      ELSE
%%  - USE_ABS=1 => retain |dX|/4/B;

dX             = GF_args{1};
dY             = GF_args{2};
B              = GF_args{3};
USE_ABS        = GF_args{4};
normal_deriv   = GF_args{5};
USE_SING       = GF_args{6};
%%
k           = pi/B;
[R,varTH]   = GEN_polar_coords(dX,dY);
jnz         = find(R);

if isempty(normal_deriv)
  ndiff  = 0;
else
  ndiff     = normal_deriv{1};
  theta_vec = normal_deriv{2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndiff==0
  %% calculate G:
  Y      = 0*dX+log(k);
  dummy  = 1-exp( k*(1i*dY-abs(dX)) );
  if USE_SING==0
    %% allow for log singularity:
    Y(jnz)  = log(abs(dummy(jnz)./R(jnz)));
  else
    %% don't need to worry about singularity:
    Y = log(abs(dummy));
  end
  Y   = Y/2/pi;
  if USE_ABS
    Y = Y+abs(dX)/4/B;
  end
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if make it this far,
%% calculate normal derivative of G:
Sign  = sign(dX);
Exp   = exp( k*(1i*dY-abs(dX)) );
Y_x   = 0*dX;
Y_y   = Y_x;
%%
if USE_SING%% RETAIN SINGULARITY:
  Y_x = k*Sign.*Exp./(1-Exp);
  Y_y = -1i*k*Exp./(1-Exp);
  %%
  Y_x = real(Y_x)/2/pi;
  Y_y = real(Y_y)/2/pi;
  %%
  if USE_ABS
    Y_x  = Y_x+Sign/4/B;
  end
else%% SUBTRACT SINGULARITY:
  subx   = dX(jnz)./(dX(jnz).^2+dY(jnz).^2);
  suby   = dY(jnz)./(dX(jnz).^2+dY(jnz).^2);
  dummy  = Exp(jnz)./(1-Exp(jnz));
  %%
  Y_x(jnz)  = k*Sign(jnz).*dummy-subx;
  Y_y(jnz)  = -1i*k*dummy-suby;
  %%
  Y_x = real(Y_x)/2/pi;
  Y_y =real(Y_y)/2/pi;
  %%
  if USE_ABS
    Y_x  = Y_x+Sign/4/B;
  end
end
%%
n1 = diag(sin(theta_vec));
n2 = -diag(cos(theta_vec));
if ndiff==1%% calc \pa_n:
  Y   = n1*Y_x+n2*Y_y;
else%% calc \pa_n':
  Y_x0   = -Y_x;
  Y_y0   = -Y_y;
  Y      = Y_x0*n1+Y_y0*n2;
end
