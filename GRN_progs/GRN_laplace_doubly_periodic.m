function Y = GRN_laplace_doubly_periodic(GF_args)
%% CALL: Y = GRN_laplace_doubly_periodic(GF_args)
%%  where GF_args={dX,dY,A,B,normal_deriv,USE_SING};
%% function calc's periodic Green's function for Laplace's eqn,
%%  - using rectangular cells of height 2B and width 2A;
%%  - normal_deriv==[] => just get G;
%%    ELSE
%%    normal_derive==th_vec => get \pa_n.G;
%%  - USE_SING==0 => subtract off log(R)/2/pi;
%%    ELSE
%%    USE_SING==1 => retain log(R)/2/pi singularity;

dX             = GF_args{1};
dY             = GF_args{2};
A              = GF_args{3};
B              = GF_args{4};
normal_deriv   = GF_args{5};
USE_SING       = GF_args{6};
%%
k     = pi/B;
eps   = 1e-10;
Nsum  = ceil( -log(eps)/k/(2*A) );

Y           = 0*dX;
Y0          = Y;
USE_SING0   = 1;%% keep singularity
                %% when adding the terms from the translated sources:
%%
USE_ABS  = 0;
for n=Nsum:-1:1
  GF0_args     = {dX+2*n*A,dY,B,USE_ABS,normal_deriv,USE_SING0};
  Y0           = GRN_laplace_singly_periodic(GF0_args);
  GF0_args{1}  = dX-2*n*A;
  Y0           = Y0+GRN_laplace_singly_periodic(GF0_args);
  Y            = Y+Y0;
end
%%
USE_ABS  = 1;
GF0_args = {dX,dY,B,USE_ABS,normal_deriv,USE_SING};
Y        = Y+GRN_laplace_singly_periodic(GF0_args);
