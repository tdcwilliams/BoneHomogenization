function C_eff=BONE_cavities_rect_cell(AB,Irr_vars,Nterms);
%% CALL: C_eff=BONE_cavities_rect_cell(AB,Irr_vars,Nterms);
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
  col=Irr_vars;
  %%
  AB=.5*[1 1];
  if 1
    crk_fxn=@CURVEprof_circarc;
    crk_prams={1};%% fraction of circle
    radius=.225;
    srt={radius*[1 1],0,[0 0]};%area=pi*radius^2
  end
  Irr_vars={crk_fxn,crk_prams,srt};
%    Irr_vars={irr_vars};
%  Nterms=5;
end

GFxn=@GRN_laplace_doubly_periodic;
GF_args={AB(1),AB(2),[],[]};

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint=300;
Nunc=2*Nterms;
[soL,ip_stuff,Mlog,w_qd]=...
   BONE_ipmatrices_cos(Nint,Nterms);
IP_stuff{1,length(ip_stuff)+1}=Mlog;
IP_stuff(1:end-1)=ip_stuff;

%% GET (x,y) AND \theta EVALUATED AT THE QUAD POINTS:
nvec=(1:Nterms)';
nvec2=(1:2*Nterms)';
ip_d1chi=ip_stuff{1}{1}(2:end,:);
ip_chi=ip_stuff{1}{2}(2:end,:);
%%
Nirregs=size(Irr_vars,1);
Ntot=Nirregs*Nunc;
bn_coeffs=zeros(Ntot,1);

for j=1:Nirregs
  irr_vars=Irr_vars(j,:);
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,Area_cavs(j,1)]=...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  if Area_cavs(j)>0
    %% NB should go clockwise round boundaries
    %%  ie normal should point inwards:
    xyvecs=fliplr(xyvecs);
    ds_dt=flipud(ds_dt);
    th_vec=flipud(th_vec)+pi;
    dtheta_ds=-flipud(dtheta_ds);
    d2s_dt2=-flipud(d2s_dt2);
    d2theta_ds2=flipud(d2theta_ds2);
    d2xy_ds2=flipud(d2xy_ds2);
    tvec=flipud(tvec);
  else
    Area_cavs(j)=-Area_cavs(j);
  end
  Xvec{j}=xyvecs(1,:)';
  Yvec{j}=xyvecs(2,:)';
  LC_irrs(j)=LC;
  Xtra(j,:)={th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ=nvec2+(j-1)*Nunc;
  bn_coeffs(JJ)=-ip_d1chi*xyvecs'*[1;1i];
end

b1coeffs=imag(bn_coeffs);
b2coeffs=-real(bn_coeffs);
Area=4*prod(AB);
Area0=Area-sum(Area_cavs);
phi0=Area0/Area;%%volume fraction of host material;
b00_const=1i*(1-phi0);

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
MK=zeros(Ntot,Ntot);
for j=1:Nirregs
  for r=1:Nirregs
    KEEPSING=(j~=r);
%      GF_args{3}={Xtra{j,1},1};
    GF_args{4}=KEEPSING;
    %% GF_args={ A,B,{th_vec,no of deriv's},KEEPSING };

    geom_stuff={soL,Xvec{j},Yvec{j},...
      LC_irrs(j),Xtra(j,:)};
    geom0_stuff={soL,Xvec{r},Yvec{r},...
      LC_irrs(r),Xtra(r,:)};
    GEOM_stuff={geom_stuff,geom0_stuff};
    %%
    JJ=nvec2+(j-1)*Nunc;
    JJ0=nvec2+(r-1)*Nunc;
    MK(JJ,JJ0)=BONE_kernelmat_cav_cos(GEOM_stuff,...
                 IP_stuff,GFxn,GF_args);
  end
end%,MK
MK(:,Ntot+1)=-b1coeffs;
MK(Ntot+1,:)=[-b1coeffs'/Area,1];%% CORRECT!!

%% CALCULATE FORCING TERMS:
%  Nunc=Nterms+1;
%  Ntot=Nirregs*Nunc;
FF=[-bn_coeffs/2; b00_const];
%  F1_coeffs=FF;
%  F2_coeffs=FF;

%  nvec=(1:Nterms+1)';
for j=1:Nirregs
  LC=LC_irrs(j);
  JJ=nvec2+(j-1)*Nunc;
  geom_stuff={soL,Xvec{j},Yvec{j},LC,Xtra(j,:)};
  %%
  for r=1:Nirregs
    KEEPSING=(j~=r);
    GF_args{4}=KEEPSING;
    %%
    LC0=LC_irrs(r);
    geom0_stuff={soL,Xvec{r},Yvec{r},LC0,Xtra(r,:)};
    GEOM_stuff={geom_stuff,geom0_stuff};
    %%
    FF(JJ)=FF(JJ)+...
      -BONE_forcing_cav_cos(GEOM_stuff,...
         {ip_stuff,w_qd},GFxn,GF_args);
  end
%    th_vec=Xtra_Tn{j,1};
%    f1_coeffs=LC*ip_chi*sin(th_vec);
%    f2_coeffs=-LC*ip_chi*cos(th_vec);
%    F1_coeffs(JJ)=f1_coeffs;
%    F2_coeffs(JJ)=f2_coeffs;
end

if 1
  load UUmpg;
  Umpg=U2mpg-1i*U1mpg;
  tst_b00eqn=...
    [MK(Ntot+1,:)*[Umpg;b00_theory],b00_const];
  tst_eqns=[FF,MK*[Umpg;b00_theory]];
  %%
  tst_MK=MK(1:Ntot,(1:10));
  tstFF=FF(1:Ntot)
  pause
end

%% SOLVE INTEGRAL EQUATION:
UU_coeffs=MK\FF;
U1_coeffs=-imag(UU_coeffs(1:Ntot));
U2_coeffs=real(UU_coeffs(1:Ntot));

if 0%% look at solution
  d1chi_vals=ip_stuff{3}{1}(:,2:end);
  U1=d1chi_vals*U1_coeffs;
  U2=d1chi_vals*U2_coeffs;
[max(U1),max(U2)],ans(1)/ans(2)
  subplot(1,2,1), plot(soL,U1,col), hold on;
  subplot(1,2,2), plot(soL,U2,col), hold on;
elseif 0%% look at forcing
  f1_coeffs=-imag(FF(1:Ntot));%[MK*U1_coeffs,f1_coeffs]
  f2_coeffs=real(FF(1:Ntot));%[MK*U2_coeffs,f2_coeffs]
  d1chi_vals=ip_stuff{3}{1}(:,2:end);
  f1=d1chi_vals*f1_coeffs;
  f2=d1chi_vals*f2_coeffs;
[max(f1),max(f2)],ans(1)/ans(2)
  subplot(1,2,1), plot(soL,f1,col), hold on;
  subplot(1,2,2), plot(soL,f2,col), hold on;
end

%% CALCULATE \bfH MATRIX;
Hmat(1,1)=b1coeffs.'*U1_coeffs/Area;%[b1coeffs,U1_coeffs]
Hmat(1,2)=b1coeffs.'*U2_coeffs/Area;%[b1coeffs,U2_coeffs]
Hmat(2,1)=b2coeffs.'*U1_coeffs/Area;%[b2coeffs,U1_coeffs]
Hmat(2,2)=b2coeffs.'*U2_coeffs/Area;%[b2coeffs,U2_coeffs]
%  Hmat(1,1)=F1_coeffs.'*U1_coeffs/Area;
%  Hmat(1,2)=F1_coeffs.'*U2_coeffs/Area;
%  Hmat(2,1)=F2_coeffs.'*U1_coeffs/Area;
%  Hmat(2,1)=F2_coeffs.'*U2_coeffs/Area;

%% CALCULATE EMPs:
rho_eff=phi0;
C44eff=phi0+Hmat(1,1);
C55eff=phi0+Hmat(2,2);
C45eff=(Hmat(1,2)+Hmat(2,1))/2;
C_eff=[C44eff,C45eff,C55eff,rho_eff];

if 1%%energy conservation test:
  Etest=imag(bn_coeffs'*UU_coeffs(1:Ntot))
end