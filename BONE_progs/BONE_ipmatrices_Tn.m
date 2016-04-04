function [svec_on_LC,ip_stuff,xtra_ints,Mlog,w_quad]=...
   BONE_ipmatrices_Tn(Nint,Nterms)

%  addpath ~/WORK/programs/matlab/GEN_progs/OP_progs;
Nint  = max(Nterms+1,Nint);
hn0   = [];

%% CALC INNER PRODUCT MATRIX FOR \psi_n''(t)=(1-t^2)^{-.5}T_n(t)/h_n;
%%  T_n are chebyshev polynomials,
%%   h_n normalise so that \int_{-1}^1\psi_m''(t)T_n(t)dt=\delta_{mn};
[svec_on_LC,w_quad]  = OP_numint_chebyshev(Nint);
[ipT,hnT,Tn_vals]    = OP_inprod_chebyshev(svec_on_LC,w_quad,Nint-1);

%%degree of approx is >= Nunc-1 = Nterms
Nunc        = Nterms+1;
ip_d2psi    = ipT(1:Nunc,:);
hn2         = hnT(1:Nunc);
d2psi_vals  = Tn_vals(:,1:Nunc)*diag(1./hn2);

%% CALC INNER PRODUCT MATRIX FOR \psi_n'(t);
N_U               = max(2,Nterms-1);
alpU              = 1;
wU                = (1-svec_on_LC.^2).*w_quad;
[ipU,hn1,Un_vals] = OP_inprod_gegenbauer(svec_on_LC,wU,alpU,N_U);
nvec              = (1:Nterms)';

if Nterms>=1%% NB HAVE REMOVED THE (1-t^2)^.5 FACTOR:
  ip_d1psi(2:Nunc,:)    = diag(-1./nvec)*ipU(1:Nterms,:);
  d1psi_vals(:,2:Nunc)  = ...
             Un_vals(:,1:Nterms)*diag(-1./(nvec.*hn1(1:Nterms)));
end
%% Now need to calc inner prod for \psi_0'(t)=1-acos(t)/pi:
%  [ip1,ip_t,ip_acos,ip_t_times_acos]=xtra_integrals_chebyshev(ipT);
[ip1,ip_t,ip_acos,ip_t_times_acos]  = OP_diffweights_chebyshev(ipT);
%  ip1t=[ip1(1:8);ip_t(1:8)+ip1(1:8)],'hi'
ip_d1psi(1,:)     = ip1-ip_acos;
d1psi_vals(:,1)   = 1-acos(svec_on_LC)./pi;

%% CALC INNER PRODUCT MATRIX FOR \psi_n(t);
alpC  = 2;
wC    = (1-svec_on_LC.^2).*wU;

if Nterms>=2%% NB HAVE NOT INCLUDED THE FACTOR OF (1-t^2)^1.5
  nvec(1)            = [];%% nvec now (2:Nterms)
  [ipC,hn0,Cn_vals]  = OP_inprod_gegenbauer(svec_on_LC,wC,alpC,Nterms-2);
  ip_psi(3:Nunc,:)   = diag(.5./nvec)*ipC;
  psi_vals(:,3:Nunc) = ...
               Cn_vals*diag(.5./(nvec.*hn0));
end
%% ip_sqrt*ff = \int_{-1}^1[f(t)*sqrt(1-t^2)/pi]dt
%%       = \int_{-1}^1[f(t)*\cU_0(t)/2]dt;
%% \psi_0(t) = t-t*acos(t)/pi+\cU_0(t)/2
ip_sqrt        = ipU(1,:)/2;
ip_psi(1,:)    = ip_t-ip_t_times_acos+ip_sqrt;
psi_vals(:,1)  = svec_on_LC.*(1-acos(svec_on_LC)/pi)+...
                sqrt(1-svec_on_LC.^2)/pi;

%% ip_t_times_sqrt*ff = \int_{-1}^1[f(t)*t*sqrt(1-t^2)/pi]dt
%%       = \int_{-1}^1[f(t)*\cU_1(t)/4]dt;
%% \psi_1(t)= -1 +acos(t)/pi - \cU_1(t)/4;
if Nterms>=1
  ip_t_times_sqrt = ipU(2,:)/4;
  ip_psi(2,:)     = -ip1+ip_acos-ip_t_times_sqrt;
  psi_vals(:,2)   = -1+acos(svec_on_LC)/pi+...
                -svec_on_LC.*sqrt(1-svec_on_LC.^2)/pi;
end

%% OUTPUTS:
ip_stuff = { {ip_d2psi,ip_d1psi,ip_psi},{hn2,hn1,hn0},...
             {d2psi_vals,d1psi_vals,psi_vals} };
%% need 1 & t as part of \bfu,w,w_n,w_s;
xtra_ints   = {ip1,ip_t};

%% L(t,t0)=(1/2/pi)*\log|t-t0|=...
%%          -log(2)/2/pi-\sum_{n=1}^\infty[T_n(t)T_n(t0)]/(n*pi)
%%  => [Mlog]_{mn}=\iint T_m(t)*L(t,t0)*T_n(t0)*dt0*dt
%%                = {-log(2) m=n=0
%%                  {-1/n*pi m=n>0
%%                  {0       m~=n
Mlog  = diag([-log(2)/2/pi;(-1/pi)./(1:Nterms).']);

if 0%% do some tests:
  tt     = svec_on_LC;
  faze   = 2*pi*rand(1);
  ff     = tt.*exp(i*(pi*tt+faze));
  %%
  if 0%% check T_n results:
    an   = ip_d2psi*ff;
    %[an,hn,Tn_vals]=OP_inprod_chebyshev(tt,w_quad,Nterms,ff);
    f_ap = Tn_vals(:,1:Nterms+1)*an;
    subplot(1,2,1), plot(tt,real([ff,f_ap]));
    subplot(1,2,2), plot(tt,real([ff,f_ap]/i));
  elseif 0 & Nterms>=1%% check U_n results:
    an   = ipU*ff;
    %[an,hn,Un_vals]=OP_inprod_gegenbauer(tt,wU,alpU,N_U,ff);
    f_ap = Un_vals*an;
    subplot(1,2,1), plot(tt,real([ff,f_ap]));
    subplot(1,2,2), plot(tt,real([ff,f_ap]/i));
  elseif 1 & Nterms>=2%% check C^(2)_n results:
    an               = ipC*ff;
    [an,hn,Cn_vals]  = OP_inprod_gegenbauer(tt,wC,alpC,Nterms-2,ff);
    f_ap             = Cn_vals*an;
    subplot(1,2,1), plot(tt,real([ff,f_ap]));
    subplot(1,2,2), plot(tt,real([ff,f_ap]/i));
  elseif 0%% look at \psi_0,\psi_1,\psi_0';
          %% formulae checked with maple (xtra_ints_cheby.mw)
    fac  = 1;
    j0   = find( abs(tt)==min(abs(tt)) );tt0=tt(j0)
    subplot(1,3,1), plot(tt,real(psi_vals(:,1)/fac));
    psi00   = psi_vals(j0,1)
    hold on,plot(tt,tt.*(1-acos(tt)/pi)+sqrt(1-tt.^2)/pi,'--r');
    %%
    subplot(1,3,2), plot(tt,real(psi_vals(:,2)/fac));
    psi10   = psi_vals(j0,2)
    hold on,plot(tt,-1+acos(tt)/pi-tt.*sqrt(1-tt.^2)/pi,'--r');
    %%
    subplot(1,3,3), plot(tt,real(d1psi_vals(:,1)/fac));
    d1psi00 = d1psi_vals(j0,1)
    hold on,plot(tt,1-acos(tt)/pi,'--r');
  end

  if 0%% test some of the extra integrals:
    tol  = 1e-12;
    f_in = inline('tt.*exp(i*pi*tt)');
    tst1 = [exp(i*faze)*quad(f_in,-1,1,tol),ip1*ff]
    %%
    f_in    = inline('tt.^2.*exp(i*pi*tt)');
    tst_t   = [exp(i*faze)*quad(f_in,-1,1,tol),ip_t*ff]
    %%
    f_in    = inline('tt.*exp(i*pi*tt).*acos(tt)/pi');
    tst_ac  = [exp(i*faze)*quad(f_in,-1,1,tol),ip_acos*ff]
    %%
    f_in    = inline('tt.^2.*exp(i*pi*tt).*acos(tt)/pi');
    tst_tac = [exp(i*faze)*quad(f_in,-1,1,tol),ip_t_times_acos*ff]
    %%
    f_in       = inline('tt.*exp(i*pi*tt).*(1-acos(tt)/pi)');
    tst_dpsi0  = [exp(i*faze)*quad(f_in,-1,1,tol),ip_d1psi(1,:)*ff]
    %%
    f_in = inline(...
      'tt.*exp(i*pi*tt).*(tt.*(1-acos(tt)/pi)+sqrt(1-tt.^2)/pi)');
    tst_psi0   = [exp(i*faze)*quad(f_in,-1,1,tol),ip_psi(1,:)*ff]
    %%
    f_in = inline(...
      'tt.*exp(i*pi*tt).*(-1+acos(tt)/pi-tt.*sqrt(1-tt.^2)/pi)');
    tst_psi1   = [exp(i*faze)*quad(f_in,-1,1,tol),ip_psi(2,:)*ff]
  end
%    rmpath ~/WORK/programs/matlab/GEN_progs/OP_progs;
end

function [ip_1,ip_t,ip_acos,ip_t_acos]=xtra_integrals_chebyshev(ipT)

N        = size(ipT,1)-1;
nn       = (0:N)';
nn_ev    = (0:2:N)';
nn_odd   = (1:2:N)';

%% ip_1*ff=\int_{1}^{1}f(t)dt=\sum_{n=0}^Nf_n\int_{-1}^{1}T_n(t)dt
int_n          = zeros(N+1,1);
int_n(nn_ev+1) = 2./(1-nn_ev.^2);
ip_1           = int_n'*ipT;

%% ip_t*ff=\int_{1}^{1}[f(t)*t]dt=\sum_{n=0}^Nf_n\int_{-1}^{1}[T_n(t)*t]dt
int_n             = zeros(N+1,1);
int_n(nn_odd+1)   = 2./(4-nn_odd.^2);
ip_t              = int_n'*ipT;

%% ip_acos*ff=\int_{1}^{1}[f(t)*acos(t)/pi]dt
%% = \sum_{n=0}^Nf_n\int_{-1}^{1}[T_n(t)*acos(t)/pi]dt
int_n          = zeros(N+1,1)-1/4;
nn_            = nn;
nn_(2)         = [];
int_n(nn_+1)   = (-1).^nn_./(1-nn_.^2);
ip_acos        = int_n'*ipT;

%% ip4*ff=\int_{1}^{1}[f(t)*t*acos(t)/pi]dt
%% = \sum_{n=0}^Nf_n\int_{-1}^{1}[T_n(t)*t*acos(t)]dt
int_n          = zeros(N+1,1)-1/16;
nn_            = nn;
nn_(3)         = [];
int_n(nn_+1)   = -(-1).^nn_./(4-nn_.^2);
ip_t_acos      = int_n'*ipT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_plot_U(soL,FF,LC)

Nterms   = length(FF)-1;
nvec     = 0:Nterms;
DD       = diag((nvec+1)/LC);
ff       =   OP_interp_gegenbauer(soL,1,DD*FF);
plot(soL,real(ff/1i),'--m');
