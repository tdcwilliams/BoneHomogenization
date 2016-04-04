function [MU_dnW,MU_dsW,MU_inner,FF,...
           MV_dnW,MV_dsW,MV_inner,xtra]=...
  BONE_MPG_kernelmat_zone_fibres(...
    geom_stuff_outer,Irr_vars,NN,ipCS)

Nterms   = NN(1);
Nint     = 150;
if length(NN)==2
  Nint   = NN(2);
end
TST1     = 0;
nn       = 3;
TST2     = 0;
USE_FS   = 1;

%  ipCS=IP_stuff{1};
N_MP     = size(ipCS,1)/2;
kn       = pi*(1:N_MP)';
DD       = diag(kn);
MD_otr   = [0*DD,DD;-DD,0*DD];%%differentiates F.series
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
%% CHOOSE GF:
if USE_FS
  GF  = [];%%->log/2/pi
else
  GF  = {@GRN_laplace_doubly_periodic,A,B};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GEOMETRY OF FIBRES:
Nunc           = 2*Nterms;
[soL,ip_stuff] = ...
   BONE_ipmatrices_cos(Nint,Nterms);
ip_d1chi       =  ip_stuff{1}{1}(2:end,:);
ip_chi         =  ip_stuff{1}{2}(2:end,:);
Np_in          = size(ip_d1chi,2);
wq_in          = 2/Np_in+0*(1:Np_in)';
kn             = pi*(1:Nterms)';
DD             = diag(-.5./kn);
Mlog_in        = [DD,0*DD;0*DD,DD];
%%
TTf2c    = 2*[0*DD,DD;-DD,0*DD];%%fib-cav map
DD       = diag(kn);
MD_inr   = [0*DD,DD;-DD,0*DD];%% differentiates F.series
TTc2f    = MD_inr;%%cav->fib map
%%
Nirregs     = size(Irr_vars,1);
Nunc        = 2*Nterms;
Ntot        = Nirregs*Nunc;
bn_coeffs   = zeros(Ntot,1);
FF          = bn_coeffs;
%%
geom_stuff_inner  = cell(Nirregs,6);

for j=1:Nirregs
  irr_vars     = Irr_vars(j,1:3);
  centre(j,:)  = irr_vars{3}{3};
  %%
  mdBC         =  Irr_vars{j,4};
  m_irrs(j,1)  = mdBC{1};%% internal shear modulus (rel to host);
  d_irrs(j,1)  = mdBC{2};%% internal density (rel to host);
  BC(j,1)      = mdBC{3};%% boundary condition; 
                         %%  BC=0 cts condition, BC=1 debonded (cavity-like);
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,area_j]   = ...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  Area_cavs(j,1)     = abs(area_j);
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
  JJ              = (1:Nunc)+(j-1)*Nunc;
  bn_coeffs(JJ)   = ip_d1chi*xyvecs'*[1;1i];
  if BC(j)==0
    bn_coeffs(JJ) = (1-m_irrs(j))*bn_coeffs(JJ);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  geom_stuff_inner(j,:) = ...
    {soL,Xvec{j},Yvec{j},LC_irrs(j),th_vec,d2xy_ds2};
    %[Xvec{j},Yvec{j}],pause
end
b1_coeffs   = imag(-bn_coeffs);%% -\xi_2: -\pa_s->\sin(\th)=n_1;
b2_coeffs   = real(bn_coeffs);%% \xi_1: -\pa_s->-\cos(\th)=n_2;
%%
xtra  = {b1_coeffs,b2_coeffs,Area_cavs,m_irrs,d_irrs};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALLOCATE MEMORY FOR MATRICES:
MU_dnW   = zeros(2*N_MP+Ntot,Np_otr);
MU_dsW   = MU_dnW;
MU_inner = zeros(2*N_MP+Ntot,Ntot);
FF       = MU_dnW(:,1);
%%
MV_dnW   = zeros(Ntot,Np_otr);
MV_dsW   = MV_dnW;
MV_inner = zeros(Ntot,Ntot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTER EQN:
%% 0=\pa_s\int_otr.G\pa_s'W.ds'
%%    -1/2.\pa_n'W -\int_otr.\pa_nG\pa_n'W.ds'
%%    +\pa_s\int_inr.G[\pa_s'W].ds'
%%    -\int_otr.\pa_nG[\pa_n'W].ds'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALC #1:
%% 1st row of MU_dsW, MU_dnW & MU_inner:
ip_left     = (MD_otr/LC_otr)*ipCS;
%% NB MD_otr/LC_otr applies \pa_s;
ip_right    = diag(LC_otr*wq_otr);
Mlog        = MD_otr*Mlog_otr*ipCS;
IP_stuff    = {ip_left,ip_right,Mlog};
GEOM_stuff  = {geom_stuff_outer,geom_stuff_outer};
JJ          = 1:2*N_MP;
JJ0         = 1:Np_otr;

if 0
   if BC==0
      load MPfib
      disp('*****testing fibres*****'), pause(.5);
   else
      disp('*****testing cavities*****'), pause(.5);
      load MPcav
      UU_coeffs_MP          = U2-1i*U1;
      Wcoeffs_MP_out        = mp_coeffs;
      Wcoeffs_MP_in         = 0*mp_coeffs;
      Wcoeffs_MP_in(1)      = 1i;
      Wcoeffs_MP_in(N_MP+1) = -1;
   end
   %%
   xyLth_otr = {xvec_otr,yvec_otr,[],theta_otr};
   xyLth_inr = geom_stuff_inner(1,2:5);
   %%
   [dnW_otr,dsW_otr,dnW_inr,dnW_inr0,dsW_inr,dsW_inr0] = ...
     get_needed_coeffs_MP(Wcoeffs_MP_out,Wcoeffs_MP_in,...
       xyLth_otr,xyLth_inr,A,B,N_MP,ip_d1chi);
   if BC==1
      UU_mpg = dsW_inr-dsW_inr0;
   else
      UU_mpg = dnW_inr-dnW_inr0;
   end
   VV_mpg    = dsW_inr-m_irrs(j)*dsW_inr0;
   YY        = {dnW_otr,dsW_otr,dnW_inr,dnW_inr0,dsW_inr,dsW_inr0};
   %%
   xvec_inr  = geom_stuff_inner{1,2};
   yvec_inr  = geom_stuff_inner{1,3};
   LC_inr    = geom_stuff_inner{1,4};
   theta_inr = geom_stuff_inner{1,5};
elseif 0
   dsW_otr   = exp(1i*theta_otr);
   dnW_otr   = dsW_otr/1i;
   %%
   theta_inr = geom_stuff_inner{1,5};
   LC_inr    = geom_stuff_inner{1,4};
   dsW_inr   = LC*ip_d1chi*( exp(1i*theta_inr) );
   dnW_inr   = dsW_inr/1i;
   %%
   dsW_inr0  = dsW_inr;
   dnW_inr0  = dnW_inr;
end

%%->+\pa_s\int_otr.G[\pa_s'W]ds'
MU_dsW(JJ,JJ0) = ...
   + BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);

%%->-1/2.\pa_n'W -\int_otr.\pa_nG\pa_n'W.ds'
IP_stuff       = {ipCS,ip_right,[]};
MU_dnW(JJ,JJ0) = -.5*ipCS+...
    -BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
%%
for j=1:Nirregs
  GEOM_stuff{2}   = geom_stuff_inner(j,:);
  JJ0             = (1:Nunc)+(j-1)*Nunc;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if BC(j)==0%% NORMAL BC:
    %%->-\int_inr.\pa_n'G[\pa_n'W]ds'
    IP_stuff         = {ipCS,ip_d1chi,[]};
    MatGn            = BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
    MU_inner(JJ,JJ0) = -MatGn;
    %%
%      IP_stuff={MD_otr/LC_otr*ipCS,ip_d1chi,[]};
%      MatGs=BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
  else%% FULLY DEBONDED:
    %%->+\pa_s\int_inr.G[\pa_s'W]ds'
    IP_stuff         = {ip_left,ip_d1chi,[]};
    %% NB ip_left has \pa_s accounted for;
    MatGs            = BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
    MU_inner(JJ,JJ0) = MatGs;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if 0
  tst_eqn_outer   = ...
    MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
     (MatGs*dsW_inr-MatGn*dnW_inr),pause
    [MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
      -MatGn*UU_coeffs_MP,...
     MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
      -MatGn*UU_coeffs_MP2];%,pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INNER EQN's:
%% (1-m_j)*e^(1i*th)=1/2(m_j+1)[W_n]
%%                   -(m_j-1)*\pa_s\int_otr.G\pa_s'W.ds'
%%                   +(m_j-1)*\pa_n\int_otr.G\pa_n'W.ds'
%%                   +(m_j-1)*\pa_n\int_inr.G[\pa_n'W].ds'
%% 1/L\int_inr.chi_m'(s/L).{******}.ds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALC #2:
%% rest of MU_dsW, MU_dnW & MU_inner
%%  (eqns for on the irregularities,
%%   from \int [G\pa_s'W]ds'):
ip_right_otr   = diag(LC_otr*wq_otr);
Mlog           = [];
JJ0            = 1:Np_otr;
BASIS          = 0;

for j=1:Nirregs
  LC  = LC_irrs(j);
  mj  = m_irrs(j);
  MD  = 1/LC*MD_inr;
  JJ  = (1:Nunc)+(j-1)*Nunc+(2*N_MP);
  %%
  if BASIS==0%% 'f' basis:
    ipl_1diff  = MD*(LC*ip_d1chi);
    ipl_0diff  = LC*ip_d1chi;
    Mlog       = MD_inr*Mlog_in;
    FF(JJ)     = MD*LC*bn_coeffs;
    Mdiag      = (mj+1)/2*eye(Nunc);
  else%% 'c' basis:
    ipl_1diff  = ip_d1chi;
    ipl_0diff  = -LC*ip_chi;
    Mlog       = Mlog_in;
    FF(JJ)     = bn_coeffs;
    Mdiag      = (mj+1)/2*TTf2c;
  end
  GEOM_stuff   = {geom_stuff_inner(j,:),...
               geom_stuff_outer};
  %%
  IP_stuff  = {ipl_0diff,ip_right_otr,[]};
  U_FACTOR  = 1+(BC(j)-1)*mj;
  V_FACTOR  = 1-BC(j)*mj;
  %%
  mGn       = BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
  %%->-\pa_n\int_otr.G\pa_s'W.ds'
  MV_dsW(JJ-2*N_MP,:)   = -V_FACTOR*mGn;
  %%->(m_j-1)*\pa_n\int_otr.G\pa_n'W.ds'
  MU_dnW(JJ,:)          = -U_FACTOR*mGn;
  %%
  IP_stuff              = {ipl_1diff,ip_right_otr,[]};
  mGs                   = BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
  %%->-\pa_s\int_otr.G\pa_s'W.ds'
  MV_dnW(JJ-2*N_MP,:)   = -V_FACTOR*mGs;
  %%->(m_j-1)*\pa_s\int_otr.G\pa_s'W.ds'
  MU_dsW(JJ,:)          = U_FACTOR*mGs;

  %%->(1/2)(m_j+1)[W_n]
  JJ0 = (1:Nunc)+(j-1)*Nunc;
  if BC(j)==0
    MU_inner(JJ,JJ0) = ...
       +Mdiag;
  else
    MV_inner(JJ-2*N_MP,JJ0)   = ...
       +Mdiag;
  end
  for r=1:Nirregs
    GEOM_stuff{2} = geom_stuff_inner(r,:);
    IP_stuff      = {ipl_0diff,ip_d1chi,[]};
    mGn_U         =   BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
%  TTf2c*mGn_U*dnW_inr,pause
    %%
    if j==r
      IP_stuff = {ipl_1diff,ip_d1chi,Mlog};
    else
      IP_stuff = {ipl_1diff,ip_d1chi,[]};
    end
    mGs_U   = BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
    %%
    JJ0     = (1:Nunc)+(r-1)*Nunc;
    %%->(m_j-1)*\pa_n\int_inr.G[\pa_n'W].ds'
    if BC(j)==0
      MU_inner(JJ,JJ0)        = ...
            MU_inner(JJ,JJ0)-(1-mj)*mGn_U;
      MV_inner(JJ-2*N_MP,JJ0) = -mGs_U;
    else
      MU_inner(JJ,JJ0)        = mGs_U;
      MV_inner(JJ-2*N_MP,JJ0) = ...
            MV_inner(JJ-2*N_MP,JJ0)-(1-mj)*mGn_U;%%%%%
    end
  end
  if BASIS==1
    MV_inner(JJ-2*N_MP,:)  = TTc2f*MV_inner(JJ-2*N_MP,:);
    MV_dsW(JJ-2*N_MP,:)    = TTc2f*MV_dsW(JJ-2*N_MP,:);
    MV_dnW(JJ-2*N_MP,:)    = TTc2f*MV_dnW(JJ-2*N_MP,:);
  end
  if 0%%test vs fibre results:
    [dsW_inr,MV_dnW*dnW_otr+MV_dsW*dsW_otr+...
        +MV_inner*(dnW_inr-dnW_inr0)];
    fn   = (1-mj)*ipl_0diff*exp(1i*theta_inr);
    [FF(JJ),fn,.5*(1+mj)*(dnW_inr-dnW_inr0)+...
        +.5*(1-mj)*(dnW_inr+dnW_inr0)];
    [fn,.5*(1+mj)*(dnW_inr-dnW_inr0)+...
        +(1-mj)*(mGs*dsW_otr-mGn*dnW_otr+...
           -mGn_U*(dnW_inr-dnW_inr0))];
    [MU_dsW(JJ,:)*dsW_otr,(1-mj)*mGs*dsW_otr]
    pause
  end
  if 0%% test vs cavity results:
    load MPGcav
%      MPGcav_mats=...
%         [FFcav,Mcav_Wn*dnW_otr,Mcav_Ws*dsW_otr,Mcav_U*dsW_inr];
    TTf2c*[.5*dnW_inr+mGn_U*dnW_inr,...
            MU_dnW(JJ,:)*dnW_otr/U_FACTOR+...
            MU_dsW(JJ,:)*dsW_otr/U_FACTOR+...
            mGs_U*dsW_inr];
    [.5*dnW_inr+mGn_U*dnW_inr,...
            MU_dnW(JJ,:)*dnW_otr/U_FACTOR+...
            MU_dsW(JJ,:)*dsW_otr/U_FACTOR+...
            mGs_U*dsW_inr];
    [dnW_inr,...
           MU_dnW(JJ,:)*dnW_otr/U_FACTOR+...
            MU_dsW(JJ,:)*dsW_otr/U_FACTOR+...
            mGs_U*(dsW_inr-dsW_inr0)];
    tst_Veqn   = [.5*(dsW_inr+dsW_inr0),...
          -mGs*dnW_otr-mGn*dsW_otr+...
          -mGn_U*(dsW_inr-dsW_inr0)];%,pause
    %%
    tst_Veqn   = [dsW_inr-mj*dsW_inr0,...
          .5*(1+mj)*(dsW_inr-dsW_inr0)+...
            +.5*(1-mj)*(dsW_inr+dsW_inr0),...
          .5*(1+mj)*(dsW_inr-dsW_inr0)+(1-mj)*(...
            -mGs*dnW_otr-mGn*dsW_otr+...
            -mGn_U*(dsW_inr-dsW_inr0)),...
          MV_dnW*dnW_otr+MV_dsW*dsW_otr+...
            +MV_inner*(dsW_inr-dsW_inr0),...
          b00_theory*(1-mj)*LC*ip_d1chi*cos(theta_inr)+...
             +MV_inner*(dsW_inr-dsW_inr0),...
          [MV_inner,(1-mj)*TTc2f*ip_d1chi*xvec_inr]*...
             [dsW_inr-dsW_inr0;b00_theory] ],pause
%        IP_stuff={-ip_d1chi,ip_d1chi,-Mlog_in};
%        mG_U=BONE_kernelmat_log(GEOM_stuff,IP_stuff,1);
%        mG_U*dsW_inr,pause

%        [dnW_inr,MU_dnW(JJ,:)*dnW_otr/U_FACTOR+...
%              MU_dsW(JJ,:)*dsW_otr/U_FACTOR+...
%              mGs_U*(dsW_inr-dsW_inr0)],pause
%  LC*MD*bn_coeffs
    tst_Ueqns  = [FF,...
        MU_dsW*dsW_otr+MU_dnW*dnW_otr+...
          MU_inner*(dsW_inr-dsW_inr0)],pause
     %%
     if 0%BC==1
       UU_mpg  = dsW_inr-dsW_inr0;
       VV_mpg  = dsW_inr-mj*dsW_inr0;
       MV_mpg  = [MV_inner,(1-mj)*TTc2f*ip_d1chi*xvec_inr];
       save MPGcav2 b00_theory UU_mpg VV_mpg MV_mpg Hmat
       tstVeqn = [VV_mpg,MV_mpg*[UU_mpg;b00_theory]];%,pause
       tst_b00_eqn=...
         [b00_theory,Hmat(1,2)+1i*(sum(Area_cavs)-Hmat(1,1)),...
            1i*sum(Area_cavs)+bn_coeffs.'*VV_mpg]
       Hmat
       [bn_coeffs.'*VV_mpg,Hmat(1,2)-1i*Hmat(1,1)]
       pause
     end
  end
end

if 0
  [dnW_inrT,dnW_inr0T,dsW_inrT,dsW_inr0T,...
    UU_mpgT,VV_mpgT,ipT,TTf2Tn_u,TTf2Tn_v,UU_P,Up_P]  = ...
      transform_exp_to_Tn(dnW_inr,dnW_inr0,dsW_inr,dsW_inr0,...
        LC,BASIS,mdBC,TTc2f,FF(JJ),irr_vars);
  if BC==1
    save MPGfib1 dnW_otr dsW_otr dnW_inrT dnW_inr0T dsW_inrT dsW_inr0T UU_mpgT VV_mpgT
    %%
%      disp('409: otr UU contrib:');
%      J_=(1:2*N_MP);
%      [MU_dsW(J_,:)*dsW_otr,...
%       MU_dnW(J_,:)*dnW_otr,... 
%       MU_inner(J_,:)*UU_mpg];
%       tst_otr=[ans,sum(ans,2)]
     disp('415: inr UU contrib:');
     TTf2Tn_u*[MU_dsW(JJ,:)*dsW_otr,...
     MU_dnW(JJ,:)*dnW_otr,... 
     MU_inner(JJ,:)*UU_mpg];
     tst_inr   = [TTf2Tn_u*FF(JJ),sum(ans,2)];
     %%
     TTf2Tn_v*[MV_dsW*dsW_otr,...
     MV_dnW*dnW_otr,...
     MV_inner*UU_mpg];
     Fac = 1+mj;
     M0  = Fac*TTf2Tn_v*...
        [MV_dsW/(1-mj)*dsW_otr, MV_dnW/(1-mj)*dnW_otr,...
         -mGn_U*UU_mpg];
     tstV   = [VV_mpgT,(1-mj)*UU_P/2+(1+mj)*Up_P/2,...
             (1-mj)*UU_P/2+sum(M0,2)];
     [Fac*Up_P/2,sum(M0,2)];
     %%
     M2     = TTf2Tn_v*MV_dnW;
     tst2   = M2(1:10,1:7);
     if 0
       fn   = MU_inner(JJ,:)*UU_mpg;
       test_plot(soL,fn,LC), hold on;
       %%
       fnT  = TTf2Tn_u*fn;
       test_plot_U(soL,fnT,LC), hold off;
     end
     %%
%       plot(soL,cos(theta_inr)), hold on;
%       test_plot(FF(JJ),soL,LC_irrs(j)), hold off;
  else
    save MPGfib0 dnW_otr dsW_otr dnW_inrT dnW_inr0T dsW_inrT dsW_inr0T UU_mpgT VV_mpgT
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DnW_inr,DnW_inr0,DsW_inr,DsW_inr0,...
  UU_mpg,VV_mpg,ip_d1chi,TT_f2Tn_u,TT_f2Tn_v,UU_P,Up_P]  = ...
    transform_exp_to_Tn(dnW_inr,dnW_inr0,dsW_inr,dsW_inr0,...
      LC,BASIS,mdBC,TTc2f,FF,irr_vars)

Nint     = 300;
Nterms   = .5*length(dnW_inr);
[soL,ip_stuff,xtra_ints,Mlog_in,wq_in] = ...
     BONE_ipmatrices_Tn(Nint,Nterms);
nvec        = (1:Nterms)';
ip_d1chi    = ip_stuff{1}{1}(nvec+1,:);
hnT         = ip_stuff{2}{1}(nvec+1);%% T_1 -> T_N
d1chi_vals  = ip_stuff{3}{1}(:,nvec+1);
Tn_vals     = d1chi_vals*diag(hnT);
ip_chi      = ip_stuff{1}{2}(nvec+1,:);
ipU         = diag(-nvec)*ip_chi;
ip1_Tn      = xtra_ints{1};
%%
hnU      = ip_stuff{2}{2}(nvec);
chi_vals = ip_stuff{3}{2}(:,nvec+1);
Un_vals  = chi_vals*diag(-nvec.*hnU);
%%
[ip_Pn,iplv_0diff_Pn_oL,iplv_1diff_Pn,...
  xtraP] = ...
    BONE_Pn_Stuff(soL,ip1_Tn,Nterms);
Pn_vals  = xtraP{1};
Pn_ends  = xtraP{2};
hnP      = xtraP{3};

if 0%% check if M1 below is invertible
    %% - it is - hooray!
  ip1 = xtra_ints{1};
  M1  = -chi_vals'*diag(ip1)*...
        ip_stuff{3}{1}(:,nvec)
  inv(M1)
  2/pi,pause
end
%%
kn       = (1:Nterms)'*pi;
CS       = [cos(soL*kn'),sin(soL*kn')]/LC;
IP0      = LC*diag(hnT)*ip_d1chi;
TT_f2Tn  = IP0*diag(sqrt(1-soL.^2))*CS;
TT_f2Pn  = LC*Pn_vals'*diag(ip1_Tn)*CS;
%%
TT_f2Tn_u   = -LC*ip_chi*CS;
TT_f2Tn_v   = LC*Pn_vals'*diag(ip1_Tn)*CS;
%%
if BASIS==1%% convert from cav basis to fib basis 1st;
  TT_f2Tn   = TT_f2Tn*TTc2f;
  TT_f2Tn_u = TT_f2Tn_u*TTc2f;
  TT_f2Tn_v = TT_f2Tn_v*TTc2f;
end
%%
DnW_inr  = TT_f2Tn*dnW_inr;
DnW_inr0 = TT_f2Tn*dnW_inr0;
DsW_inr  = TT_f2Tn*dsW_inr;
DsW_inr0 = TT_f2Tn*dsW_inr0;
%%
mj    = mdBC{1};
BC    = mdBC{3};
Hig   = dsW_inr-mj*dsW_inr0;
if BC==0
  UU_mpg = DnW_inr-DnW_inr0;
  UU     =dnW_inr-dnW_inr0;
  VV     =dsW_inr0;
  Up     =dsW_inr+dsW_inr0;
else
  UU_mpg = DsW_inr-DsW_inr0;
  UU     = dsW_inr-dsW_inr0;
  VV     = dsW_inr+mj*dsW_inr0;
  Up     = dsW_inr+dsW_inr0;
%    [Hig,2*mj/(1+mj)*UU+(1-mj)/(1+mj)*VV],pause
end
VV_mpg   = TT_f2Pn*VV;
UU_P     = TT_f2Pn*UU;
Up_P     = TT_f2Pn*Up;
%  [(1+mj)*Up_P/2,VV_mpg-(1-mj)*UU_P/2],pause

if 1
  xyvecs = BONE_get_rsdtheta_NRquick(...
           irr_vars,soL );
  yvec   = xyvecs(2,:)';
  xvec   = xyvecs(1,:)';
  %%

  b1_coeffs = -ip_d1chi*yvec;
  b2_coeffs = ip_d1chi*xvec;
  c1_coeffs = -ip_Pn*yvec;
  c2_coeffs = ip_Pn*xvec;
  bb        = 2*mj/(1+mj)*[b1_coeffs,b2_coeffs]';
  cc        = (1-mj)/(1+mj)*[c1_coeffs,c2_coeffs]';
  H0        = bb*UU_mpg + cc*VV_mpg;
  Hmat      = real(H0*[1i,1]);
  %%
%    {LC*CS'
  d1_coeffs = -(LC*CS')*diag(ip1_Tn)*yvec;
  d2_coeffs = (LC*CS')*diag(ip1_Tn)*xvec;
  dd        = [d1_coeffs,d2_coeffs];
  Hmat      = real(dd'*Hig*[1i,1]);

%    plot(soL,[-yvec,(LC*CS)*d1_coeffs,Tn_vals*b1_coeffs,Pn_vals*c1_coeffs]),pause
%   plot(soL,[xvec,(LC*CS)*d2_coeffs,0*Tn_vals*b2_coeffs,Pn_vals*c2_coeffs]),pause
end

if 0%% check forcing;
  ff  = exp(-1i*pi*soL)/1i;%%circle
%    plot(soL,real(ff/1i)), hold on;
  test_plot(soL,FF,LC), hold on; pause
  %%
%    fnU=ipU*ff;
%    ffU=Un_vals*fnU;pause
%    plot(soL,real(ffU/1i),'--r');pause
  %%
  fnT = -LC*ip_chi*ff;[fnT,TT_f2Tn_u*FF]
  ffT = Un_vals*diag(nvec/LC)*fnT;
  test_plot_U(soL,fnT,LC);
  hold off, pause;
end

if 0
  f1  = CS*(dsW_inr+mj*dsW_inr0);
  f2  = Pn_vals*(VV_mpg./hnP)/LC;
  plot(soL,real([f1,f2]/1i)); pause
end

Ye = {dnW_inr,dnW_inr0,dsW_inr,dsW_inr0};
Yt = {DnW_inr,DnW_inr0,DsW_inr,DsW_inr0};
for j=[]
  q1  = CS*Ye{j};
  p1  = q1.*sqrt(1-soL.^2);
  plot(soL,real(q1)), hold on;
  %%
  p2  = (Tn_vals*diag(1./(LC*hnT)))*Yt{j};
  q2  = p2./sqrt(1-soL.^2);
  [Yt{j},  LC*diag(hnT)*ip_d1chi*p1]
  plot(soL,real(q2),'--r'), hold off;
  pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnW_otr,dsW_otr,dnW_inr,dnW_inr0,dsW_inr,dsW_inr0]=...
  get_needed_coeffs_MP(Wcoeffs_MP_out,Wcoeffs_MP_in,...
  xyLth_otr,xyLth_inr,A,B,N_MP,ip_d1chi)

[dxW_otr,dyW_otr] = ...
     BONE_diff_multipoles_laplace_doubly_periodic(...
             xyLth_otr{1},xyLth_otr{2},N_MP,[A B],[0 0]);
cto      = diag(cos(xyLth_otr{4}));
sto      = diag(sin(xyLth_otr{4}));
dsW_otr  = (cto*dxW_otr+sto*dyW_otr);
dnW_otr  = (sto*dxW_otr-cto*dyW_otr);
dsW_otr  = [real(dsW_otr),imag(dsW_otr)]*Wcoeffs_MP_out;
dnW_otr  = [real(dnW_otr),imag(dnW_otr)]*Wcoeffs_MP_out;
%%
[dxW_inr,dyW_inr] = ...
     BONE_diff_multipoles_laplace_doubly_periodic(...
             xyLth_inr{1},xyLth_inr{2},N_MP,[A B],[0 0]);
cti      = diag(cos(xyLth_inr{4}));
sti      = diag(sin(xyLth_inr{4}));
dsW_inr  = (cti*dxW_inr+sti*dyW_inr);
dnW_inr  = (sti*dxW_inr-cti*dyW_inr);
dsW_inr  = [real(dsW_inr),imag(dsW_inr)]*Wcoeffs_MP_out;
dnW_inr  = [real(dnW_inr),imag(dnW_inr)]*Wcoeffs_MP_out;
%%
[dxW_inr0,dyW_inr0]  = ...
       BONE_diff_multipoles_laplace(...
             xyLth_inr{1},xyLth_inr{2},N_MP,[0 0]);
dsW_inr0 = (cti*dxW_inr0+sti*dyW_inr0);
dnW_inr0 = (sti*dxW_inr0-cti*dyW_inr0);
dsW_inr0 = [real(dsW_inr0),imag(dsW_inr0)]*Wcoeffs_MP_in;
dnW_inr0 = [real(dnW_inr0),imag(dnW_inr0)]*Wcoeffs_MP_in;
%%
LC       = xyLth_inr{3};
dnW_inr  = LC*ip_d1chi*dnW_inr;
dnW_inr0 = LC*ip_d1chi*dnW_inr0;
dsW_inr  = LC*ip_d1chi*dsW_inr;
dsW_inr0 = LC*ip_d1chi*dsW_inr0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_plot(soL,FF,LC)

Nterms   = length(FF)/2;
kn       = (1:Nterms)*pi;
CS       = [cos(soL*kn),sin(soL*kn)];
ff       = CS*(FF/LC);
plot(soL,real(ff/1i),'--g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_plot_U(soL,FF,LC)

Nterms   = length(FF)-1;
nvec     = 0:Nterms;
DD       = diag((nvec+1)/LC);
ff       = OP_interp_gegenbauer(soL,1,DD*FF);
plot(soL,real(ff/1i),'--m');
