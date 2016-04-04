function Y = GRN_helmholtz_doubly_periodic(GF_args)
%% CALL: Y = GRN_laplace_doubly_periodic(GF_args)
%%  where GF_args={dX,dY,A,B,normal_deriv,USE_SING};
%% function calc's periodic Green's function for Helmholtz's eqn,
%%  - using rectangular cells of height 2B and width 2A;
%%  - normal_deriv==[] => just get G;
%%    ELSE
%%    normal_derive==th_vec => get \pa_n.G;
%%  - USE_SING==0 => subtract off log(R)/2/pi;
%%    ELSE
%%    USE_SING==1 => retain log(R)/2/pi singularity;

DOTEST   = 1;

if nargin==0%% test inputs
   A              = 1;
   B              = 2;
   %%
   xx             = linspace(-A,A,30)';
   yy             = B/3+0*xx;
   GF_args        = {xx,yy};
   dX             = A/pi*angle(pi/A*GF_args{1});%%need dX\in[-A,A]
   dY             = B/pi*angle(pi/B*GF_args{2});%%need dY\in[-B,B]
   %%
   k              = 1.2;
   qq             = [0 0];
   normal_deriv   = [];
   USE_SING       = 1;
else
   A              = GF_args{3};
   B              = GF_args{4};
   %%
   dX             = A/pi*angle(pi/A*GF_args{1});%%need dX\in[-A,A]
   dY             = B/pi*angle(pi/B*GF_args{2});%%need dY\in[-B,B]
   %%
   k              = GF_args{5};
   qq             = GF_args{6};
   normal_deriv   = GF_args{7};
   USE_SING       = GF_args{8};
end

if ~isempty(normal_deriv)
   if normal_deriv{1}~=0
      %%haven't programmed normal derivative yet;
      disp('Have not programmed the normal derivative of this GF yet.')
      disp('Exiting GRN_helmholtz_doubly_periodic.m.')
      return;
   end
end

%%
a           = 2*A;
b           = 2*B;
q1          = qq(1);
q2          = qq(2);
qna_shift   = q1*(GF_args{1}-dX);
qnb_shift   = q2*(GF_args{2}-dY);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get singular bit of function;
%% (default periodicity is in y direction,
%%  decays exponentially as |x|\ar0);
fn_Gtilde   = @GRN_laplace_singly_periodic;
USE_ABS     = 1;
%%
R        = abs(dX+1i*dY);
Y        = 0*R;
jzero    = find(R==0);
if ~isempty(jzero)
   %% get limit as x,y\ar0, after taking off 1/2/pi*log(R);
   N0       = 20;%% sum decays like n^(-9);
   Y(jzero) = G0(k,q1,a,b,N0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_xdecay = find((R>0)&(abs(dX)>=abs(dY)));
if ~isempty(j_xdecay)
   x           = dX(j_xdecay);
   y           = dY(j_xdecay);
   N           = 15;%% sum decays exponentially;
   Y(j_xdecay) = G_xdecay(x,y,k,q1,a,b,N);
   %%
   Gtilde_args = {x,y,B,USE_ABS,normal_deriv,USE_SING};
   Gtilde      = feval(fn_Gtilde,Gtilde_args);
   Y(j_xdecay) = exp(1i*q2*y).*(Y(j_xdecay)+Gtilde);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_ydecay = find((R>0)&(abs(dX)<abs(dY)));
if ~isempty(j_ydecay)
   x           = dX(j_ydecay);
   y           = dY(j_ydecay);
   N           = 15;%% sum decays exponentially;
   Y(j_ydecay) = G_xdecay(y,x,k,q2,b,a,N);
   %%
   Gtilde_args = {y,x,A,USE_ABS,normal_deriv,USE_SING};
   Gtilde      = feval(fn_Gtilde,Gtilde_args);
   Y(j_ydecay) = exp(1i*q1*x).*(Y(j_ydecay)+Gtilde);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out  = G_xdecay(x,y,k,q1,a,b,N)

z0    = exp(1i*a*(k+q1));
z0b   = exp(1i*a*(k-q1));
out   = 1/2i/k/b*( exp(1i*k*x)-1i*k*abs(x)+...
         +z0*exp(-1i*k*x)/(1-z0)+z0b*exp(1i*k*x)/(1-z0b) );

%% S_0 sum;
nn       = (1:N)';
[NN,X]   = meshgrid(nn,x);
[NN,Y]   = meshgrid(nn,y);
kn       = 2*pi*nn/b;
Kn       = k*lam_fxn(kn/k);
%%
bet   = k*b/2/pi;
S0n   = -1/b./Kn+nn.^(-1)/2/pi+...
         +bet^2/4/pi*( nn.^(-3)+3/4*bet^2*nn.^(-5)+5/8*bet^4*nn.^(-7) );
out   = out + sum(S0n)- bet^2/4/pi*( zeta(3)+3/4*bet^2*zeta(5)+5/8*bet^4*zeta(7) );

%% exponentially decaying part of sum;
kn    = 2*pi*NN/b;
Kn    = k*lam_fxn(kn/k);
%%
zn    = exp(a*(1i*q1-Kn));
znb   = exp(-a*(1i*q1+Kn));
Zn    = exp(a*(1i*q1-Kn)+Kn.*X);
Znb   = exp(-a*(1i*q1+Kn)-Kn.*X);
%%
Sn    = exp(-Kn*abs(X))./Kn-exp(-kn*abs(X))./kn+...
         +(Zn./Kn)./(1-zn) + (Znb./Kn)./(1-znb);
         {out,Sn.*cos(kn.*Y)}
out   = out - 1/b*sum( Sn.*cos(kn.*Y),2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y  = G0(k,q1,a,b,N)

z0    = exp(1i*a*(k+q1));
z0b   = exp(1i*a*(k-q1));
Y     = 1/2/pi*log(2*pi/b)+1/2i/k/b*( 1+z0/(1-z0)+z0b/(1-z0b) );

%% S_0 sum;
nn    = (1:N)';
kn    = 2*pi*nn/b;
Kn    = k*lam_fxn(kn/k);
%%
bet   = k*b/2/pi;
S0n   = -1/b./Kn+nn.^(-1)/2/pi+...
         +bet^2/4/pi*( nn.^(-3)+3/4*bet^2*nn.^(-5)+5/8*bet^4*nn.^(-7) );
Y     = Y + sum(S0n)- bet^2/4/pi*( zeta(3)+3/4*bet^2*zeta(5)+5/8*bet^4*zeta(7) );

%% exponentially decaying part of sum;
zn    = exp(a*(1i*q1-Kn));
znb   = exp(-a*(1i*q1+Kn));
Y     = Y - 1/b*sum( (zn./Kn)./(1-zn) + (znb./Kn)./(1-znb) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = lam_fxn(t)

y  = conj(sqrt(t.^2-1));
%% we need sqrt(t.^2-1)=-1i*sqrt(1-t.^2)
%% - this is the opposite of matlab's sqrt fxn so take complex conjugate; 
