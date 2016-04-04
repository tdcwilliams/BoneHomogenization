function [MU_dnW,MU_dsW,MU_inner,FF,...
           MV_dnW,MV_dsW,MV_inner,xtra]=...
  BONE_MPG_kernelmat_zone_fibres_Tn(...
    geom_stuff_outer,Irr_vars,NN,ipCS)

Nterms=NN(1);
Nint=150;
if length(NN)==2
  Nint=NN(2);
end
TST1=0; nn=3;
TST2=0;
USE_FS=1;

%  ipCS=IP_stuff{1};
N_MP=size(ipCS,1)/2;
kn=pi*(1:N_MP)';
DD=diag(kn);
MD_otr=[0*DD,DD;-DD,0*DD];%%differentiates F.series
DD=diag(-.5./kn);
Mlog_otr=[DD,0*DD;0*DD,DD];
%  ip_dCS=MD*ip_CS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GEOMETRY OF BUFFER ZONE:
%% geom_stuff_outer=...
%%   {soL,LC,xvec,yvec,LC_otr,theta_otr,A,B};
soL_otr=geom_stuff_outer{1};
xvec_otr=geom_stuff_outer{2};
yvec_otr=geom_stuff_outer{3};
LC_otr=geom_stuff_outer{4};
theta_otr=geom_stuff_outer{5};
Np_otr=length(soL_otr);
wq_otr=2/Np_otr+0*soL_otr;
A=geom_stuff_outer{7};
B=geom_stuff_outer{8};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSE GF:
if USE_FS
  GF=[];%%->log/2/pi
else
  GF={@GRN_laplace_doubly_periodic,A,B};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GEOMETRY OF FIBRES:
Nunc=Nterms;
[soL,ip_stuff,xtra_ints,Mlog_in,wq_in]=...
   BONE_ipmatrices_Tn(Nint,Nterms);
ip_d1chi=ip_stuff{1}{1}(2:end,:);
ip_chi=ip_stuff{1}{2}(2:end,:);
ip1_Tn=xtra_ints{1};
%%
[ip_Pn,iplv_0diff_Pn_oL,iplv_1diff_Pn,xtraP]=...
  BONE_Pn_Stuff(soL,ip1_Tn,Nterms);
Pn_vals=xtraP{1};
Pn_ends=xtraP{2};
%%
nvec=(1:Nterms)';
hnT=ip_stuff{2}{1}(nvec+1);%% T_1 -> T_N
Tn_vals=ip_stuff{3}{1}(:,nvec+1)*diag(hnT);
hnU=ip_stuff{2}{2}(nvec);%% U_0->U_{N-1}
%  {ip_stuff{3}{2}(:,nvec),diag(-nvec.*hnU)}
Un_vals=ip_stuff{3}{2}(:,nvec+1)*diag(-nvec.*hnU);
dTn_vals=Un_vals*diag(nvec);
  %% T_n'(t)=nU_{n-1}(t);
%%
Nirregs=size(Irr_vars,1);
Nunc=Nterms;
Ntot=Nirregs*Nunc;
bn_coeffs=zeros(Ntot,1);
cn_coeffs=zeros(Ntot,1);
FF=bn_coeffs;
%%
geom_stuff_inner=cell(Nirregs,6);

for j=1:Nirregs
  irr_vars=Irr_vars(j,1:3);
  centre(j,:)=irr_vars{3}{3};
  %%
  mdBC=Irr_vars{j,4};
  m_irrs(j,1)=mdBC{1};%% internal shear modulus (rel to host);
  d_irrs(j,1)=mdBC{2};%% internal density (rel to host);
  BC(j,1)=mdBC{3};%% boundary condition;
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,area_j]=...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  Area_irrs(j,1)=abs(area_j);
  if area_j>0
    %% NB normal should point INTO the cavity:
    xyvecs=fliplr(xyvecs);
    ds_dt=flipud(ds_dt);
    th_vec=flipud(th_vec)+pi;
    dtheta_ds=-flipud(dtheta_ds);
    d2s_dt2=-flipud(d2s_dt2);
    d2theta_ds2=flipud(d2theta_ds2);
    d2xy_ds2=flipud(d2xy_ds2);
    tvec=flipud(tvec);
  end
  Xvec{j}=xyvecs(1,:)';
  Yvec{j}=xyvecs(2,:)';
  LC_irrs(j)=LC;
  Xtra(j,:)={th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ=(1:Nunc)+(j-1)*Nunc;
  bn_coeffs(JJ)=ip_d1chi*xyvecs'*[1;1i];
  cn_coeffs(JJ)=ip_Pn*xyvecs'*[1;1i];
  mj=m_irrs(j);
  if BC(j)==0
    cn_coeffs(JJ)=(1-mj)*cn_coeffs(JJ);
    bn_coeffs(JJ)=0;
  else
    bn_coeffs(JJ) = 2*mj/(1+mj)*bn_coeffs(JJ);
    cn_coeffs(JJ) = (1-mj)/(1+mj)*cn_coeffs(JJ);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  geom_stuff_inner(j,:)=...
    {soL,Xvec{j},Yvec{j},LC_irrs(j),th_vec,d2xy_ds2};
    %[Xvec{j},Yvec{j}],pause
  %%
  soL_ends=[-1;1];
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2]=...
        BONE_get_rsdtheta_NRquick( irr_vars,soL_ends );
  for r=1:2
    geom_stuff_inner_ends{j,r}=...
      {soL_ends(r),xyvecs(1,r),xyvecs(2,r),LC_irrs(j),...
        th_vec(r),d2xy_ds2(r)};
  end
end
b1_coeffs=imag(-bn_coeffs);%% -\xi_2: -\pa_s->\sin(\th)=n_1;
b2_coeffs=real(bn_coeffs);%% \xi_1: -\pa_s->-\cos(\th)=n_2;
c1_coeffs=imag(-cn_coeffs);
c2_coeffs=real(cn_coeffs);
%%
xtra={{b1_coeffs,c1_coeffs},{b2_coeffs,c2_coeffs},...
        Area_irrs,m_irrs,d_irrs};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALLOCATE MEMORY FOR MATRICES:
MU_dnW=zeros(2*N_MP+Ntot,Np_otr);
MU_dsW=MU_dnW;
MU_inner=zeros(2*N_MP+Ntot,Ntot);
FF=MU_dnW(:,1);
%%
MV_dnW=zeros(Ntot,Np_otr);
MV_dsW=MV_dnW;
MV_inner=zeros(Ntot,Ntot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTER EQN:
%% 0=\pa_s\int_otr.G\pa_s'W.ds'
%%    -1/2.\pa_n'W -\int_otr.\pa_nG\pa_n'W.ds'
%%    +\pa_s\int_inr.G[\pa_s'W].ds'
%%    -\int_otr.\pa_nG[\pa_n'W].ds'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALC #1:
%% 1st row of MU_dsW, MU_dnW & MU_inner:
ip_left=(MD_otr/LC_otr)*ipCS;
%% NB MD_otr/LC_otr applies \pa_s;
ip_right=diag(LC_otr*wq_otr);
Mlog=MD_otr*Mlog_otr*ipCS;
IP_stuff={ip_left,ip_right,Mlog};
GEOM_stuff={geom_stuff_outer,geom_stuff_outer};
JJ=1:2*N_MP;
JJ0=1:Np_otr;

%%->+\pa_s\int_otr.G[\pa_s'W]ds'
MU_dsW(JJ,JJ0)=...
   + BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);

%%->-1/2.\pa_n'W -\int_otr.\pa_nG\pa_n'W.ds'
IP_stuff={ipCS,ip_right,[]};
MU_dnW(JJ,JJ0)=-.5*ipCS+...
    -BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
%%
for j=1:Nirregs
  GEOM_stuff{2}=geom_stuff_inner(j,:);
  JJ0=(1:Nunc)+(j-1)*Nunc;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if BC==0%% NORMAL BC:
    %%->-\int_inr.\pa_n'G[\pa_n'W]ds'
    IP_stuff={ipCS,ip_d1chi,[]};
    MatGn=BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
    MU_inner(JJ,JJ0)=-MatGn;
    %%
%      IP_stuff={MD_otr/LC_otr*ipCS,ip_d1chi,[]};
%      MatGs=BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
  else%% FULLY DEBONDED:
    %%->+\pa_s\int_inr.G[\pa_s'W]ds'
    IP_stuff={ip_left,ip_d1chi,[]};
    %% NB ip_left has \pa_s accounted for;
    IS_CLOSED=0;
    MatGs=BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,...
            IS_CLOSED,GF);
    MU_inner(JJ,JJ0)=MatGs;
    %%
    IP_stuff={ipCS,ip_d1chi,[]};
    MatGn=BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if 0%BC==1
  load MPGfib1
  disp('184: otr UU contrib:')
    [MU_dsW(JJ,:)*dsW_otr,...
     MU_dnW(JJ,:)*dnW_otr,... 
     MU_inner(JJ,:)*UU_mpgT];
  tst_otr=[ans,sum(ans,2)];

  tst_eqn_outer=...
    [MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
     (MatGs*dsW_inrT-MatGn*dnW_inrT),...
     MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
     MatGs*UU_mpgT]%pause
%      [MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
%        -MatGn*UU_coeffs_MP,...
%       MU_dnW(JJ,:)*dnW_otr+MU_dsW(JJ,:)*dsW_otr+...
%        -MatGn*UU_coeffs_MP2];%,pause
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
ip_right_otr=diag(LC_otr*wq_otr);
Mlog=[];
JJ0=1:Np_otr;
%%

for j=1:Nirregs
  LC=LC_irrs(j);
  mj=m_irrs(j);
  JJ=(1:Nunc)+(j-1)*Nunc+(2*N_MP);
  %%
  if 0
    iplv_0diff=LC*Tn_vals'*diag(ip1_Tn);
    iplv_1diff=-dTn_vals'*diag(ip1_Tn);
    OP_ends=[(-1).^nvec,1+0*nvec];
  else
    iplv_0diff=LC*iplv_0diff_Pn_oL;
    iplv_1diff=iplv_1diff_Pn;
    OP_ends_v=Pn_ends;
    OP_ends_u=[(-1).^nvec,1+0*nvec];
  end
  %%
  if BC(j)==1
    iplu_1diff=ip_d1chi;
    iplu_0diff=-LC*ip_chi;
    Mlog=Mlog_in(2:end,2:end);
    Mdiag0=ip_d1chi*Pn_vals;
    Mdiag=(1-mj)/2*Mdiag0';
    %%
    U_FACTOR=1;
    V_FACTOR=1+mj;
    FF(JJ)=iplu_1diff*U_FACTOR*[Xvec{j},Yvec{j}]*[1;1i];
  else
    iplu_0diff=LC*Tn_vals'*diag(ip1_Tn);
    iplu_1diff=-dTn_vals'*diag(ip1_Tn);
    Mlog=iplv_1diff*Tn_vals*Mlog_in(2:end,2:end);
    Mdiag=(mj+1)/2*eye(Nterms);
    %%
    U_FACTOR=1-mj;
    V_FACTOR=1;
    theta=Xtra{j,1};
    FF(JJ)=iplu_0diff*U_FACTOR*exp(1i*theta);
  end
  %%
  GEOM_stuff={geom_stuff_inner(j,:),geom_stuff_outer};
  IP_stuff={1,ip_right_otr,[]};
  %%
  mGn=BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
  MU_dnW(JJ,:)=-U_FACTOR*iplu_0diff*mGn;
  %%->-\pa_n\int_otr.G\pa_s'W.ds'
  MV_dsW(JJ-2*N_MP,:)=-V_FACTOR*iplv_0diff*mGn;
  %%->(m_j-1)*\pa_n\int_otr.G\pa_n'W.ds'

  %%
  IP_stuff={1,ip_right_otr,[]};
  IS_CLOSED=1;
  mGs=BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,...
        IS_CLOSED,GF);
  %%->(m_j-1)*\pa_s\int_otr.G\pa_s'W.ds'
  MU_dsW(JJ,:)=U_FACTOR*iplu_1diff*mGs;

  %%->-\pa_s\int_otr.G\pa_n'W.ds'
  MV_dnW(JJ-2*N_MP,:)=-V_FACTOR*iplv_1diff*mGs;

  %%
  IP_stuff={1,ip_right_otr,[]};
  GEOM_stuff_end=...
    {geom_stuff_inner_ends{j,2},geom_stuff_outer};
  Gp1=BONE_kernelmat_arbGF(GEOM_stuff_end,IP_stuff,...
        IS_CLOSED,GF);
  %%
  IP_stuff={1,ip_right_otr,[]};
  GEOM_stuff_end=...
    {geom_stuff_inner_ends{j,1},geom_stuff_outer};
  Gm1=BONE_kernelmat_arbGF(GEOM_stuff_end,IP_stuff,...
        IS_CLOSED,GF);

  for r=1:Nterms
    mGp1(r,:)=Gp1;
    mGm1(r,:)=Gm1;
  end
  MV_dnW(JJ-2*N_MP,:)=MV_dnW(JJ-2*N_MP,:)+...
        -V_FACTOR*( diag( OP_ends_v(:,2) )*mGp1 +...
                    -diag( OP_ends_v(:,1) )*mGm1 );
  %%->(1/2)(m_j+1)[W_n]
  JJ0=(1:Nunc)+(j-1)*Nunc;
  if BC==0
    MU_inner(JJ,JJ0)=Mdiag;
    MU_dsW(JJ,:)=MU_dsW(JJ,:)+...
        +U_FACTOR*( diag( OP_ends_u(:,2) )*mGp1 +...
                    -diag( OP_ends_u(:,1) )*mGm1 );
  else
    MV_inner(JJ-2*N_MP,JJ0)=Mdiag;
  end

  for r=1:Nirregs
    GEOM_stuff={geom_stuff_inner(j,:),...
                 geom_stuff_inner(r,:)};
%  TTf2c*mGn_U*dnW_inr,pause
    %%
    JJ0=(1:Nunc)+(r-1)*Nunc;
    %%->(m_j-1)*\pa_n\int_inr.G[\pa_n'W].ds'
    if BC(j)==0
      IP_stuff={iplu_0diff,ip_d1chi,[]};
      mGn_U=BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
      %%
      if j==r
        IP_stuff={iplv_1diff,ip_d1chi,Mlog};
      else
        IP_stuff={iplv_1diff,ip_d1chi,[]};
      end
      IS_CLOSED=0;
      mGs_U=BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,...
              IS_CLOSED,GF);
      %%
      MU_inner(JJ,JJ0)=...
            MU_inner(JJ,JJ0)-U_FACTOR*mGn_U;
      MV_inner(JJ-2*N_MP,JJ0)=-mGs_U;
      %%
      Mlog2=( diag(Mlog_in(2:end,2:end)).*OP_ends_u(:,2) )';
      IP_stuff={1,ip_d1chi,Mlog2};
      GEOM_stuff_end=...
        {geom_stuff_inner_ends{j,2},geom_stuff_inner(r,:)};
      Gp1=BONE_kernelmat_arbGF(GEOM_stuff_end,IP_stuff,...
        IS_CLOSED,GF);
      %%
      Mlog1=( diag(Mlog_in(2:end,2:end)).*OP_ends_u(:,1) )';
      IP_stuff={1,ip_d1chi,Mlog1};
      GEOM_stuff_end=...
        {geom_stuff_inner_ends{j,1},geom_stuff_inner(r,:)};
      Gm1=BONE_kernelmat_arbGF(GEOM_stuff_end,IP_stuff,...
        IS_CLOSED,GF);
      %%
      for it=1:Nterms
        mvGp1(it,:)=Gp1;
        mvGm1(it,:)=Gm1;
      end
      MV_inner(JJ-2*N_MP,JJ0)=MV_inner(JJ-2*N_MP,JJ0)+...
        -V_FACTOR*( diag( OP_ends_v(:,2) )*mvGp1 +...
                    -diag( OP_ends_v(:,1) )*mvGm1 );

    else
      IP_stuff={iplv_0diff,ip_d1chi,[]};
      mGn_U=BONE_kernelmat_dn_arbGF(GEOM_stuff,...
              IP_stuff,GF);
      %%
      if j==r
        IP_stuff={iplu_1diff,ip_d1chi,Mlog};
%          IP_stuff={ip_d1chi,ip_d1chi,Mlog};
      else
        IP_stuff={iplu_1diff,ip_d1chi,[]};
      end

      IS_CLOSED=0;
      mGs_U=BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,...
              IS_CLOSED,GF);
      %%
      MU_inner(JJ,JJ0)=mGs_U;
      MV_inner(JJ-2*N_MP,JJ0)=...
            MV_inner(JJ-2*N_MP,JJ0)-V_FACTOR*mGn_U;%%%%%
    end
  end
end

if 0%BC==1
  load MPGfib1
  %%
  J_=1:2*N_MP;%{MU_dsW(J_,:}
  [MU_dsW(J_,:)*dsW_otr,...
    MU_dnW(J_,:)*dnW_otr,...
     MU_inner(J_,:)*UU_mpgT];
  tst_otr=[ans,sum(ans,2)];

  tst_inr0=[MU_dsW(JJ,:)*dsW_otr,...
             MU_dnW(JJ,:)*dnW_otr,...
              MU_inner(JJ,:)*UU_mpgT];
  tst_inr=[FF(JJ),sum(tst_inr0,2)];
%    test_plot_U(soL,MU_inner(JJ,:)*UU_mpgT,LC)
  %%
  tstV0=[MV_dsW*dsW_otr,...
             MV_dnW*dnW_otr,...
              MV_inner*UU_mpgT];
  tstV1=(1-mj)/(1+mj)*tstV0;
  tst2=MV_dnW(1:10,1:7);
  tstV=[VV_mpgT,sum(tstV0,2)];
  tstV(:,1)-tstV(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y,Gn,fn0]=test_ig_Gnfn0(theta0,r0,z,nn)

%f=(x0+1i*y0)^n=r0^n*e^(1i*n*th0)
fn0=-nn*r0^(nn-1)*exp(1i*nn*theta0);
ds0=r0;
%  dx=real(z-r0*exp(theta0));
%  dy=imag(z-r0*exp(theta0));
R=abs(z-r0*exp(1i*theta0));
r=abs(z);
Gn=(1/2/pi./R.^2).*(r-r0/r*real(z*exp(-1i*theta0)));
  %%z=r*exp(1i*theta)=>z_n=+z_r=exp(1i*theta)
Y=Gn.*(fn0*ds0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_plot_U(soL,FF,LC)

Nterms=length(FF)-1;
nvec=0:Nterms;
DD=diag((nvec+1)/LC);
ff=OP_interp_gegenbauer(soL,1,DD*FF);
plot(soL,real(ff/1i),'k');