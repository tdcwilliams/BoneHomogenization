function C_eff=BONE_MPG_fibres_rect_cell_Tn(AB,Irr_vars_outer,...
  Irr_vars_inner,N_MP,Nterms);
%% CALL: C_eff=BONE_MPG_fibres_rect_cell(AB,Irr_vars_outer,...
%%              Irr_vars_inner,N_MP,Nterms);
%%
%% DESCRIPTION:
%% Works out the effective material properties for a
%%  periodic elastic medium (rectangular unit cell)
%%  when there are N smooth cavities per unit cell;
%% Solves the first order cell problem to get
%%  a homogenized zero-th order wave eqn;
%% Uses a Fourier series expansion to solve the problem;
%%
%% INPUTS:
%% AB=[A,B] gives the non-dim sizes of the unit cell
%%  (A is the height, B the width);
%% Irr_vars{j}={function handle, curve parameters,...
%%              [scalers, rotation angle, translation]}
%%   for the j-th cavity;
%% Nterms = order of Fourier approximation
%%  for the displacement around the crack;
%%
%% OUTPUTS:
%% Effective material properties
%%  C_eff=[C^*_44,C^*_45,C^*_55]

%  om=physvars(1);
%  mu=physvars(2);
%  rho=physvars(3);
%  %%
%  k0sq=rho*om^2/mu;
%  k0=sqrt(k0sq);
%  Lnd=1/k0;
%  ab=ab_dim/Lnd;

if nargin==2%%use some test inputs:
  Nterms=AB;
  N_MP=Nterms;
  BC=Irr_vars_outer;
  %%
  AB=.5*[1 1];
  if 1
    crk_fxn=@CURVEprof_circarc;
    crk_prams_out={[1,+1]};%% fraction of circle, anticlockwise;
    radius_out=.4
    srt_out={radius_out*[1 1],0,[0 0]};%area=pi*radius^2
    %%
    crk_prams_in={[1,-1]};%% fraction of circle, clockwise;
    radius_in=.225
    srt_in={radius_in*[1 1],0,[0 0]};%area=pi*radius^2
    mj=10;
    dj=1;
    fib_stuff={mj,dj,BC};
  end
  Irr_vars_outer={crk_fxn,crk_prams_out,srt_out};
  Irr_vars_inner{1}={crk_fxn,crk_prams_in,srt_in,fib_stuff};
%    Irr_vars={irr_vars};
%  Nterms=5;
end

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint=250;
[soL,ip_stuff]=...
   BONE_ipmatrices_cos(Nint,N_MP);
ipCS=ip_stuff{1}{1}(2:end,:);
%%
Nzones=length(Irr_vars_inner);
Ntot_outer=Nzones*(2*N_MP);
Ntot_inner=zeros(Nzones,1);
Nunc=2*N_MP;

USE_TN=1;

for j=1:Nzones
  %% put inner stuff aside for later:
  irr_vars_inner=Irr_vars_inner{j};
  Ntot_inner(j)=(2-USE_TN)*Nterms*size(irr_vars_inner,1);

  %% get geometries of boundaries of buffer zones
  %% NB should go anti-clockwise round these boundaries
  %%  ie normal should point outwards:
  irr_vars=Irr_vars_outer(j,:);
  centre(j,:)=irr_vars{3}{3};
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,Area_buffer]=...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  if Area_buffer<0
    xyvecs=fliplr(xyvecs);
    ds_dt=flipud(ds_dt);
    th_vec=flipud(th_vec)+pi;
    dtheta_ds=-flipud(dtheta_ds);
    d2s_dt2=-flipud(d2s_dt2);
    d2theta_ds2=flipud(d2theta_ds2);
    d2xy_ds2=flipud(d2xy_ds2);
    tvec=flipud(tvec);
  end
%    FAC_OUT(j,1)=-1+2*(Area_buffer>0);
%  tst=xyvecs(:,1:10)'
%  tst=ds_dt(1:10)
%  tst=exp(1i*th_vec(1:10))
%  tst=dtheta_ds(1:10)%!!-ve
%  tst=d2s_dt2(1:10)
%  tst=d2theta_ds2(1:10)
%  tst=d2xy_ds2(1:10,:)
%  tst=tvec(1:10)

  Xvec{j}=xyvecs(1,:)';
  Yvec{j}=xyvecs(2,:)';
  LC_irrs(j)=LC;
  Xtra(j,:)={th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  GEOM_stuff(j,:)=...
    {soL,xyvecs(1,:)',xyvecs(2,:)',LC,th_vec,...
       d2xy_ds2};
end

if 0%%plot the cavities:
  A=AB(1);
  B=AB(2);
  xy=[-B,-A;-B,A;B,A;B,-A;-B,-A];
  plot(xy(:,1),xy(:,2),'k'), hold on;
  %%
  for j=1:Nirregs
    plot(Xvec{j},Yvec{j},'r');
  end
  return;
end


%% CALCULATE KERNEL MATRIX:
Ntot=Ntot_outer+sum(Ntot_inner);
%  MK_outer=zeros(Ntot,Ntot_outer);
%  MK_inner=zeros(Ntot,sum(Ntot_inner));
MU=zeros(Ntot,Ntot);
MV=MU;
MV(1:Ntot_outer,:)=[];
%%
FF=zeros(Ntot,1);
Nj=Nunc+Ntot_inner;
J_inner=[];
%%
b1_coeffs=[];
b2_coeffs=[];
c1_coeffs=[];
c2_coeffs=[];
%%
Area_cavs=[];
m_irrs=[];
d_irrs=[];
%%
MP_stuff={[],AB};

for j=1:Nzones
  geom_stuff=[GEOM_stuff(j,:),{AB(1),AB(2)}];
  if USE_TN==0
    [MU_dnW,MU_dsW,MU_inner,Fj,...
     MV_dnW,MV_dsW,MV_inner,xtra]=...
       BONE_MPG_kernelmat_zone_fibres(...
         geom_stuff,Irr_vars_inner{j},Nterms,ipCS);
  else
    [MU_dnW,MU_dsW,MU_inner,Fj,...
     MV_dnW,MV_dsW,MV_inner,xtra]=...
       BONE_MPG_kernelmat_zone_fibres_Tn(...
         geom_stuff,Irr_vars_inner{j},Nterms,ipCS);
  end
  %%
  JJ=(1:Nj(j))+sum(Nj(1:j-1));
  JV=(1:Ntot_inner(j))+sum(Ntot_inner(1:j-1));
  JJ0_inner=(1:Ntot_inner(j))'+Nunc+sum(Nj(1:j-1));
  J_inner=[J_inner;JJ0_inner];
  %%
  FF(JJ)=Fj;
  MU(JJ,JJ0_inner)=MU_inner;
  MV(JV,JJ0_inner)=MV_inner;
  b1_coeffs=[b1_coeffs;xtra{1}{1}];
  b2_coeffs=[b2_coeffs;xtra{2}{1}];
  c1_coeffs=[c1_coeffs;xtra{1}{2}];
  c2_coeffs=[c2_coeffs;xtra{2}{2}];
  %%
  Area_cavs=[Area_cavs;xtra{3}];
  m_irrs=[m_irrs;xtra{4}];
  d_irrs=[d_irrs;xtra{5}];
  %%
  for r=1:Nzones
    JJ0_outer=(1:Nunc)+sum(Nj(1:r-1));
    MP_stuff{1}=centre(r,:);
    [dnW,dsW]=BONE_MP_kernelmatV2_dbl_periodic(...
              geom_stuff,MP_stuff,N_MP);
    MU(JJ,JJ0_outer)=MU_dnW*dnW+MU_dsW*dsW;
    MV(JV,JJ0_outer)=MV_dnW*dnW+MV_dsW*dsW;
%  tst=M_dnW*dnW+0*M_dsW*dsW;
%  tst=dnW;
%  rt=1;jt=(1:Nterms)+(rt-1)*Nterms;[tst(jt,1:3),tst(jt,end-2:end)]
%  [Nterms,size(M_dnW,1)]
  end
end

%% SOLVE SYSTEM:
MP_coeffs=MU\FF;
Area=4*prod(AB);
phi=Area_cavs/Area;
phi0=1-sum(phi);
%%
VV_coeffs=MV*MP_coeffs;
UU_coeffs=MP_coeffs(J_inner);
MP_coeffs(J_inner)=[];
%%
if 0
  load MPGfib1;
  [UU_mpgT,UU_coeffs]
  %%
  UU_coeffs=UU_mpgT;
  VV_coeffs=VV_mpgT;
end

%% CALCULATE \bfH MATRIX;
%%  Hmat=(1/Area)*\oint(\bfn\bfW^T)\rmd s
%%  NB \bfn is the normal pointing INTO the cavity;
if 0
  V2_coeffs=real(VV_coeffs);
  V1_coeffs=-imag(VV_coeffs);
  %%
  U1_coeffs=-imag(UU_coeffs);
  U2_coeffs=real(UU_coeffs);
  %%
  Hmat(1,1)=b1_coeffs.'*U1_coeffs/Area+...
            +c1_coeffs.'*V1_coeffs/Area;
  Hmat(1,2)=b1_coeffs.'*U2_coeffs/Area+...
            +c1_coeffs.'*V2_coeffs/Area;
  Hmat(2,1)=b2_coeffs.'*U1_coeffs/Area+...
            +c2_coeffs.'*V1_coeffs/Area;
  Hmat(2,2)=b2_coeffs.'*U2_coeffs/Area+...
            +c2_coeffs.'*V2_coeffs/Area;
else
  bb=[b1_coeffs,b2_coeffs]';
  cc=[c1_coeffs,c2_coeffs]';
  H0=bb*UU_coeffs + cc*VV_coeffs;
  Hmat=real(H0*[1i,1])/Area;
end
%  Hmat
if 0
  U1mpg=U1_coeffs;
  U2mpg=U2_coeffs;
  b00_theory=Hmat(1,2)+1i*(1-phi0-Hmat(1,1))
  save UUmpg U1mpg U2mpg b00_theory Hmat
end


%% CALCULATE EMPs:
rho_eff=phi0+phi'*d_irrs;
if 0
  C55eff=1+Hmat(1,1)+phi'*(m_irrs-1);
  C44eff=1+Hmat(2,2)+phi'*(m_irrs-1);
  C45eff=(Hmat(1,2)+Hmat(2,1))/2;
  C_eff=[C55eff,C45eff,C44eff,rho_eff];
else
  m_av=phi0+phi'*m_irrs;
  Mmat=m_av*eye(2)+.5*(Hmat+Hmat');
  %% Mmat=[m_x,m_xy;m_xy,m_y];
  C_eff={Mmat,rho_eff};
end

if 0%%energy conservation test:
  Etest=imag(bn_coeffs'*UU_coeffs(1:Ntot))
elseif 0
  C_eff
  C_eff2=BONE_MP_cavs_rect_cell(AB,Irr_vars_inner{1},N_MP)
end
