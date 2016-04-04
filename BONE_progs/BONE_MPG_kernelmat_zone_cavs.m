function [M_dnW,M_dsW,FF,MK_inner,xtra]=...
  BONE_MPG_kernelmat_zone_cavs(...
    geom_stuff_outer,Irr_vars,NN,ipCS)

USE_FS   = 1;%% USE_FS=1=> USE FUNDAMENTAL SOLUTION AS GF,
         %% USE_FS=0=> USE DOUBLY PERIODIC GF;

Nterms   = NN(1);
Nint     = 50;
if length(NN)==2
  Nint   = NN(2);
end
TST1  = 0;
nn    = 3;
TST2  = 0;

%  ipCS=IP_stuff{1};
N_MP     = size(ipCS,1)/2;
kn       = pi*(1:N_MP)';
DD       = diag(kn);
MD_otr   = [0*DD,DD;-DD,0*DD];
DD       = diag(-.5./kn);
Mlog_otr = [DD,0*DD;0*DD,DD];
% ip_dCS = MD*ip_CS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GEOMETRY OF BUFFER ZONE:
%% geom_stuff_outer=...
%%   {soL,LC,xvec,yvec,LC_otr,theta_otr,A,B};
soL_otr     = geom_stuff_outer{1};
xvec_otr    = geom_stuff_outer{2};
yvec_otr    = geom_stuff_outer{3};
LC_otr      = geom_stuff_outer{4};
theta_otr   = geom_stuff_outer{5};
Np_otr      = length(soL_otr);
wq_otr      = 2/Np_otr+0*soL_otr;
A           = geom_stuff_outer{7};
B           = geom_stuff_outer{8};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GEOMETRY OF CAVITIES:
Nunc           = 2*Nterms;
[soL,ip_stuff] = ...
   BONE_ipmatrices_cos(Nint,Nterms);
ip_d1chi = ip_stuff{1}{1}(2:end,:);
ip_chi   = ip_stuff{1}{2}(2:end,:);
Np_in    = size(ip_d1chi,2);
wq_in    = 2/Np_in+0*(1:Np_in)';
kn       = pi*(1:Nterms)';
DD       = diag(-.5./kn);
Mlog_in  = [DD,0*DD;0*DD,DD];
%%
Nirregs           = size(Irr_vars,1);
Nunc              = 2*Nterms;
Ntot              = Nirregs*Nunc;
bn_coeffs         = zeros(Ntot,1);
geom_stuff_inner  = cell(Nirregs,6);

for j=1:Nirregs
  irr_vars     = Irr_vars(j,:);
  centre(j,:)  = irr_vars{3}{3};
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,area_j]   = ...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  Area_cavs(j,1)  = abs(area_j);
  if area_j>0
    %% NB normal should point INTO the cavity:
    xyvecs        = fliplr(xyvecs);
    ds_dt         = flipud(ds_dt);
    th_vec        = flipud(th_vec)+pi;
    dtheta_ds     = -flipud(dtheta_ds);
    d2s_dt2       = -flipud(d2s_dt2);
    d2theta_ds2   = flipud(d2theta_ds2);
    d2xy_ds2      = flipud(d2xy_ds2);
    tvec          = flipud(tvec);
  end
  Xvec{j}      = xyvecs(1,:)';
  Yvec{j}      = xyvecs(2,:)';
  LC_irrs(j)   = LC;
  Xtra(j,:)    = {th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ                    = (1:Nunc)+(j-1)*Nunc;
  bn_coeffs(JJ)         = ip_d1chi*xyvecs'*[1;1i];
  geom_stuff_inner(j,:) = ...
    {soL,Xvec{j},Yvec{j},LC_irrs(j),th_vec,d2xy_ds2};
end
b1_coeffs   = imag(-bn_coeffs);%% -\xi_2: -\pa_s->\sin(\th)=n_1;
b2_coeffs   = real(bn_coeffs);%% \xi_1: -\pa_s->-\cos(\th)=n_2;
Area        = sum(Area_cavs);
xtra        = {b1_coeffs,b2_coeffs,Area};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALLOCATE MEMORY FOR MATRICES:
M_dnW    = zeros(2*N_MP+Ntot,Np_otr);
M_dsW    = M_dnW;
FF       = [M_dnW(1:2*N_MP,1);bn_coeffs/2];%bn_coeffs/2
MK_inner = zeros(2*N_MP+Ntot,Ntot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALC #1:
%% 1st row of M_dsW & MK_inner
%%  (eqns for boundary of buffer zone,
%%   from \int [G\pa_s'W]ds'):
ip_left     = (MD_otr/LC_otr)*ipCS;
ip_right    = diag(LC_otr*wq_otr);
Mlog        = MD_otr*Mlog_otr*ipCS;
IP_stuff    = {ip_left,ip_right,Mlog};
GEOM_stuff  = {geom_stuff_outer,geom_stuff_outer};
JJ          = 1:2*N_MP;
JJ0         = 1:Np_otr;

if USE_FS
  M_dsW(JJ,JJ0)   = ...
     BONE_kernelmat_log(GEOM_stuff,IP_stuff,1);
  %%
  for j=1:Nirregs
    IP_stuff(2:3)    = {ip_d1chi,[]};
    GEOM_stuff{2}    = geom_stuff_inner(j,:);
    JJ0              = (1:Nunc)+(j-1)*Nunc;
    MK_inner(JJ,JJ0) = ...
        BONE_kernelmat_log(GEOM_stuff,IP_stuff,1);
  end
else
  GF              = {@GRN_laplace_doubly_periodic,A,B};
  M_dsW(JJ,JJ0)   = ...
     BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
%       M_dsW(JJ,JJ0),

  %%
  for j=1:Nirregs
    IP_stuff(2:3)    = {ip_d1chi,[]};
    GEOM_stuff{2}    = geom_stuff_inner(j,:);
    JJ0              = (1:Nunc)+(j-1)*Nunc;
    MK_inner(JJ,JJ0) = ...
        BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALC #2:
%% rest of M_dsW & MK_inner
%%  (eqns for on the irregularities,
%%   from \int [G\pa_s'W]ds'):
ip_right       = diag(LC_otr*wq_otr);
Mlog           = [];
GEOM_stuff{2}  = geom_stuff_outer;
JJ0            = 1:Np_otr;
if USE_FS
  for j=1:Nirregs
    IP_stuff      = {ip_d1chi,ip_right,Mlog};
    GEOM_stuff{1} = geom_stuff_inner(j,:);
    JJ            = (1:Nunc)+(j-1)*Nunc+(2*N_MP);
    M_dsW(JJ,JJ0) = ...
        BONE_kernelmat_log(GEOM_stuff,IP_stuff,1);
    for r=1:Nirregs
      if j==r
        IP_stuff  = {ip_d1chi,ip_d1chi,Mlog_in};
      else
        IP_stuff  = {ip_d1chi,ip_d1chi,[]};
      end
      GEOM_stuff{2}     = geom_stuff_inner(r,:);
      JJ0               = (1:Nunc)+(r-1)*Nunc;
      MK_inner(JJ,JJ0)  = ...
        BONE_kernelmat_log(GEOM_stuff,IP_stuff,1);
    end
  end
else
  for j=1:Nirregs
    IP_stuff      = {ip_d1chi,ip_right,Mlog};
    GEOM_stuff{1} = geom_stuff_inner(j,:);
    JJ            = (1:Nunc)+(j-1)*Nunc+(2*N_MP);
    M_dsW(JJ,JJ0) = ...
        BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
    for r=1:Nirregs
      if j==r
        IP_stuff  = {ip_d1chi,ip_d1chi,Mlog_in};
      else
        IP_stuff  = {ip_d1chi,ip_d1chi,[]};
      end
      GEOM_stuff{2}     = geom_stuff_inner(r,:);
      JJ0               = (1:Nunc)+(r-1)*Nunc;
      MK_inner(JJ,JJ0)  = ...
        BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALC #3:
%% 1st row of M_dnW & FF
%%  (eqns for boundary of buffer zone,
%%   from -\int[\pa_nG\pa_n'W]ds'):

ip_left     = ipCS;
ip_right    = diag(LC_otr*wq_otr);
IP_stuff    = {ip_left,ip_right,[]};
GEOM_stuff  = {geom_stuff_outer,geom_stuff_outer};
JJ          = 1:2*N_MP;
JJ0         = 1:Np_otr;

if USE_FS
  M_dnW(JJ,JJ0)   = -.5*ipCS+...
      -BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff);
  %%
  for j=1:Nirregs
    theta0     = geom_stuff_inner{j,5};
    ip_right   = (LC_irrs(j)*wq_in.*exp(1i*theta0))';
      %% NB need to take transpose & complex conjugate
    if TST1
      zz       = geom_stuff_inner{j,2}+1i*geom_stuff_inner{j,3};
      g_s      = nn*zz.^(nn-1).*exp(1i*theta0);
      g_n      = g_s/1i;
      ip_right = (LC_irrs(j)*wq_in.*g_n)';
    end

    IP_stuff      = {ip_left,ip_right,[]};
    GEOM_stuff{2} = geom_stuff_inner(j,:);
    FF(JJ)        = FF(JJ)+...
        +BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff);
  end
else
  M_dnW(JJ,JJ0)   = -.5*ipCS+...
      -BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
  %%
  for j=1:Nirregs
    theta0     = geom_stuff_inner{j,5};
    ip_right   = (LC_irrs(j)*wq_in.*exp(1i*theta0))';
      %% NB need to take transpose & complex conjugate
    if TST1
      zz       = geom_stuff_inner{j,2}+1i*geom_stuff_inner{j,3};
      g_s      = nn*zz.^(nn-1).*exp(1i*theta0);
      g_n      = g_s/1i;
      ip_right = (LC_irrs(j)*wq_in.*g_n)';
    end

    IP_stuff      = {ip_left,ip_right,[]};
    GEOM_stuff{2} = geom_stuff_inner(j,:);%FF(JJ)
    FF(JJ)        = FF(JJ)+...
        +BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALC #4:
%% rest of M_dnW & FF
%%  (eqns for on the irregularities,
%%   from -\int[\pa_nG\pa_n'W]ds'):
if USE_FS
  for j=1:Nirregs
    ip_left       = -LC_irrs(j)*ip_chi;%% NB minus sign
    ip_right      = diag(LC_otr*wq_otr);
    IP_stuff      = {ip_left,ip_right,[]};
    GEOM_stuff    = {geom_stuff_inner(j,:),geom_stuff_outer};
    JJ            = 2*N_MP+(1:Nunc)+(j-1)*Nunc;
    JJ0           = 1:Np_otr;
    M_dnW(JJ,JJ0) = ...
      - BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff);
    %%
    for r=1:Nirregs
      theta0   = geom_stuff_inner{r,5};
      ip_right = (LC_irrs(r)*wq_in.*exp(1i*theta0))';
        %% NB need to take complex conjugate
      IP_stuff       = {ip_left,ip_right,[]};
      GEOM_stuff{2}  = geom_stuff_inner(r,:);
      FF(JJ)   = FF(JJ)+...
          + BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff);%FF(JJ)
    end
  end
  if 1
    Mcav_Wn = M_dnW(JJ,:);
    Mcav_Ws = M_dsW(JJ,:);
    Mcav_U  = MK_inner(JJ,:);
    FFcav   = FF(JJ);
    save MPGcav Mcav_Wn Mcav_Ws Mcav_U FFcav
  end
else
  for j=1:Nirregs
    ip_left       = -LC_irrs(j)*ip_chi;%% NB minus sign
    ip_right      = diag(LC_otr*wq_otr);
    IP_stuff      ={ip_left,ip_right,[]};
    GEOM_stuff    = {geom_stuff_inner(j,:),geom_stuff_outer};
    JJ            = 2*N_MP+(1:Nunc)+(j-1)*Nunc;
    JJ0           = 1:Np_otr;
    M_dnW(JJ,JJ0) = ...
      - BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
    %%
    for r=1:Nirregs
      theta0   = geom_stuff_inner{r,5};
      ip_right = (LC_irrs(r)*wq_in.*exp(1i*theta0))';
        %% NB need to take complex conjugate
      IP_stuff       = {ip_left,ip_right,[]};
      GEOM_stuff{2}  = geom_stuff_inner(r,:);
      FF(JJ)         = FF(JJ)+...
          + BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
  load UUmpg;
  Umpg   = U2mpg-1i*U1mpg;
  %%
  J_inner_eqns = (2*N_MP+1):length(FF);
  Jie          = J_inner_eqns;
  %%
  tst_ie = [FF(Jie),-b00_theory*b1_coeffs+MK_inner(Jie,:)*Umpg];
  tstMK  = MK_inner(Jie,(1:10))
  tstFF  = FF(Jie);
  pause
end



if TST1|TST2%%test outer matrices
%    Jo=(1:2*N_MP)';
  Jo  = 2*N_MP+(1:Nunc);
  M1  = M_dsW(Jo,:);
  M2  = M_dnW(Jo,:);
  %%
  zz  = xvec_otr+1i*yvec_otr;
%    nn=2;
  if TST1
    f_s  = nn*zz.^(nn-1).*exp(1i*theta_otr);
    f_n  = f_s/1i;
  else
    load MPcav
    centre  = [0 0];
    AB      = .5*[1 1];
    %%
    [W_x,W_y]  = BONE_diff_multipoles_laplace_doubly_periodic(...
       xvec_otr,yvec_otr,N_MP,AB,centre);
    W_x  = [real(W_x),imag(W_x)]*mp_coeffs;
    W_y  = [real(W_y),imag(W_y)]*mp_coeffs;
    f_n  = sin(theta_otr).*W_x-...
          cos(theta_otr).*W_y;
    f_s  = cos(theta_otr).*W_x+...
          sin(theta_otr).*W_y;
    %%
    [W_x,W_y]  = BONE_diff_multipoles_laplace_doubly_periodic(...
       Xvec{1},Yvec{1},N_MP,AB,centre);
    W_x  = [real(W_x),imag(W_x)]*mp_coeffs;
    W_y  = [real(W_y),imag(W_y)]*mp_coeffs;
    g_n  = sin(th_vec).*W_x-...
          cos(th_vec).*W_y;
    g_s  = cos(th_vec).*W_x+...
          sin(th_vec).*W_y;
%      plot(soL,[g_n,exp(1i*th_vec)]);
  end
%    [M1*f_s,-M2*f_n]
  M3  = MK_inner(Jo,:);
  [FF(Jo),M1*f_s+M2*f_n+M3*LC_irrs(1)*ip_d1chi*g_s]
  %%
  r0  = .3;
%    for j=1:length(xvec_otr)
%      iGnfn0(j,1)=...
%        quad(@(theta0)test_ig_Gnfn0(theta0,r0,zz(j),nn),0,2*pi);
%      [Y,Gn(:,j),fn0]=test_ig_Gnfn0(-pi*soL,r0,zz(j),nn);
%      iGnfn0(j,2)=sum(pi*wq_in.*Y);
%    end
%    [FF(Jo),ipCS*iGnfn0]
%    plot(soL_otr,iGnfn0)
%    plot(soL,[fn0,g_n])
%    plot(soL_otr,Gn(:,10),'--r');z01=r0*exp(1i*pi*soL(1))
  return;
end

function [Y,Gn,fn0]  = test_ig_Gnfn0(theta0,r0,z,nn)

%f=(x0+1i*y0)^n=r0^n*e^(1i*n*th0)
fn0   = -nn*r0^(nn-1)*exp(1i*nn*theta0);
ds0   = r0;
%  dx = real(z-r0*exp(theta0));
%  dy = imag(z-r0*exp(theta0));
R  = abs(z-r0*exp(1i*theta0));
r  = abs(z);
Gn = (1/2/pi./R.^2).*(r-r0/r*real(z*exp(-1i*theta0)));
  %%z=r*exp(1i*theta)=>z_n=+z_r=exp(1i*theta)
Y  = Gn.*(fn0*ds0);
