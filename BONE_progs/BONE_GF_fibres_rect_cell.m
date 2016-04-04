function C_eff=BONE_GF_fibres_rect_cell(AB,Irr_vars,Nterms);
%% CALL: C_eff=BONE_GF_fibres_rect_cell(AB,Irr_vars,Nterms);
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
%%              [scalers, rotation angle, translation],...
%%              {m_j,d_j,BC}  }
%%   for the j-th fibre/cavity (BC=0 for fibre, BC=1 for cavity);
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
  Nterms = AB;
  BC     = Irr_vars;
  %%
  AB=.5*[1 1];
  if 1
    crk_fxn    = @CURVEprof_circarc;
    crk_prams  = {1};%% fraction of circle
    radius     = .225;
    srt        = {radius*[1 1],0,[0 0]};%area=pi*radius^2
    mj         = 10;
    dj         = 1;
    fib_stuff  = {mj,dj,BC};
  end
  Irr_vars  = {crk_fxn,crk_prams,srt,fib_stuff};
%    Irr_vars={irr_vars};
%  Nterms=5;
end

%  GFxn=@GRN_laplace_doubly_periodic;
%  GF_args={AB(1),AB(2),[],[]};
GF={@GRN_laplace_doubly_periodic,AB(1),AB(2)};

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint=300;
Nunc=2*Nterms;
[soL,ip_stuff,Mlog,w_qd]=...
   BONE_ipmatrices_cos(Nint,Nterms);
IP_stuff{1,length(ip_stuff)+1}=Mlog;
IP_stuff(1:end-1)=ip_stuff;
%%
kn=pi*(1:Nterms)';
DD=diag(1./kn);
TTf2c=[0*DD,-DD;DD,0*DD];%%fib-cav map
DD=diag(kn);
MD_inr=[0*DD,DD;-DD,0*DD];%% differentiates F.series
TTc2f=MD_inr;%%cav->fib map

%% GET (x,y) AND \theta EVALUATED AT THE QUAD POINTS:
nvec=(1:Nterms)';
nvec2=(1:2*Nterms)';
ip_d1chi=ip_stuff{1}{1}(2:end,:);
ip_chi=ip_stuff{1}{2}(2:end,:);
%%
Nirregs=size(Irr_vars,1);
Ntot=Nirregs*Nunc;
bn_coeffs=zeros(Ntot,1);
vU_b00=bn_coeffs;
vV_b00=bn_coeffs;

for j=1:Nirregs
  irr_vars=Irr_vars(j,1:3);
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,Area_fibs(j,1)]=...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  %%
  mdBC         = Irr_vars{j,4};
  m_irrs(j,1)  = mdBC{1};%% internal shear modulus (rel to host);
  d_irrs(j,1)  = mdBC{2};%% internal density (rel to host);
  BC(j,1)      = mdBC{3};%% boundary condition (BC=0 for fibre, BC=1 for cavity);
  %%
  if Area_fibs(j)>0
    %% NB should go clockwise round boundaries
    %%  ie normal should point inwards:
    xyvecs        = fliplr(xyvecs);
    ds_dt         = flipud(ds_dt);
    th_vec        = flipud(th_vec)+pi;
    dtheta_ds     = -flipud(dtheta_ds);
    d2s_dt2       = -flipud(d2s_dt2);
    d2theta_ds2   = flipud(d2theta_ds2);
    d2xy_ds2      = flipud(d2xy_ds2);
    tvec          = flipud(tvec);
  else
    Area_fibs(j)  = -Area_fibs(j);
  end
  Xvec{j}      = xyvecs(1,:)';
  Yvec{j}      = xyvecs(2,:)';
  LC_irrs(j)   = LC;
  Xtra(j,:)    = {th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ  = nvec2+(j-1)*Nunc;
  bn  = ip_d1chi*xyvecs'*[1;1i];
  if BC(j)==0
    bn_coeffs(JJ) = (1-m_irrs(j))*bn;
    vU_b00(JJ)    = (1-m_irrs(j))*imag(bn);
    vV_b00(JJ)    = TTc2f*real(bn);
  else
    bn_coeffs(JJ) = bn;
    vU_b00(JJ)    = imag(bn);
    vV_b00(JJ)    = TTc2f*(1-m_irrs(j))*real(bn);
  end
  %%
  Geom_Stuff(j,:) = ...
    {soL,Xvec{j},Yvec{j},LC_irrs(j),th_vec,d2xy_ds2};
end

b1coeffs    = imag(-bn_coeffs);
b2coeffs    = real(bn_coeffs);
Area        = 4*prod(AB);
phi         = Area_fibs/Area;
phi0        = 1-sum(phi);%%volume fraction of host material;
b00_const   = 1i*(1-phi0-m_irrs.'*phi);

if 0%%plot the cavities:
  A   = AB(1);
  B   = AB(2);
  xy  = [-B,-A;-B,A;B,A;B,-A;-B,-A];
  plot(xy(:,1),xy(:,2),'k'), hold on;
  %%
  for j=1:Nirregs
    plot(Xvec{j},Yvec{j},'r');
  end
  return;
end


%% CALCULATE KERNEL MATRIX:
MU = zeros(Ntot,Ntot);
MV = MU;
FF = MU(:,1);
%%


ipl_0diff   = -LC*ip_chi;
ipl_1diff   = ip_d1chi;
Mlog        = Mlog(2:end,2:end);
%%
for j=1:Nirregs
  JJ     = nvec2+(j-1)*Nunc;
  FF(JJ) = bn_coeffs(JJ);
  %%
  mj     = m_irrs(j);
  Mdiag  = (mj+1)/2*TTf2c;
  if BC(j)==0
    MU(JJ,JJ)  = Mdiag;
  else
    MV(JJ,JJ)  = Mdiag;
  end
  %%
  for r=1:Nirregs
    KEEPSING   = (j~=r);
%      GF_args{3}={Xtra{j,1},1};
    GF_args{4} = KEEPSING;
    %% GF_args={ A,B,{th_vec,no of deriv's},KEEPSING };

    geom_stuff    = {soL,Xvec{j},Yvec{j},...
      LC_irrs(j),Xtra(j,:)};
    geom0_stuff   = {soL,Xvec{r},Yvec{r},...
      LC_irrs(r),Xtra(r,:)};
    GEOM_stuff    = {geom_stuff,geom0_stuff};
    GEOM_stuff    = {Geom_Stuff(j,:),Geom_Stuff(r,:)};
    %%
    IP_stuff   = {ipl_0diff,ip_d1chi,[]};
    mGn_U      = BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
    %%
    if j==r
      IP_stuff = {ipl_1diff,ip_d1chi,Mlog};
    else
      IP_stuff = {ipl_1diff,ip_d1chi,[]};
    end
    mGs_U   = BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,1,GF);
    %%
    JJ0  = nvec2+(r-1)*Nunc;
    if BC(j)==0
      MU(JJ,JJ0)  = ...
            MU(JJ,JJ0)-(1-mj)*mGn_U;
      MV(JJ,JJ0)  = TTc2f*( -mGs_U );
    else
      MU(JJ,JJ0)  = mGs_U;
      MV(JJ,JJ0)  = TTc2f*(...
            MV(JJ,JJ0)-(1-mj)*mGn_U );%%%%%
    end
  end
end%,MK
MU             = [MU,vU_b00];
MV             = [MV,vV_b00];
mb00           = -b1coeffs'/Area*MV;
mb00(end)      = mb00(end)+1;
MU(Ntot+1,:)   = mb00;%% CORRECT!!
FF             = [FF;b00_const];


%% SOLVE INTEGRAL EQUATION:
UU_coeffs   = MU\FF;
VV_coeffs   = MV*UU_coeffs;
V1_coeffs   = -imag(VV_coeffs);
V2_coeffs   = real(VV_coeffs);

if 0
  load MPGcav2
  unc       = [UU_mpg;b00_theory];
  tstVeqn   = ...
    [VV_mpg,MV_mpg*unc,MV*unc]
  [FF,MU*unc]
  V1_coeffs = -imag(VV_mpg);
  V2_coeffs = real(VV_mpg);

  [b00_theory,...
    1i*(1-phi0-m_irrs.'*phi)+b1coeffs'/Area*VV_mpg],
  pause

elseif 0
  load UUmpg;
  Umpg         = U2mpg-1i*U1mpg;
  tst_b00eqn   = ...
    [MK(Ntot+1,:)*[Umpg;b00_theory],b00_const];
  tst_eqns     = [FF,MK*[Umpg;b00_theory]]
  %%
  tst_MK = MK(1:Ntot,(1:10));
  tstFF  = FF(1:Ntot);
  pause
end


if 0%% look at solution
  d1chi_vals   = ip_stuff{3}{1}(:,2:end);
  U1           = d1chi_vals*U1_coeffs;
  U2           = d1chi_vals*U2_coeffs;
[max(U1),max(U2)],ans(1)/ans(2)
  subplot(1,2,1), plot(soL,U1,col), hold on;
  subplot(1,2,2), plot(soL,U2,col), hold on;
elseif 0%% look at forcing
  f1_coeffs    = -imag(FF(1:Ntot));%[MK*U1_coeffs,f1_coeffs]
  f2_coeffs    = real(FF(1:Ntot));%[MK*U2_coeffs,f2_coeffs]
  d1chi_vals   = ip_stuff{3}{1}(:,2:end);
  f1           = d1chi_vals*f1_coeffs;
  f2           = d1chi_vals*f2_coeffs;
[max(f1),max(f2)],ans(1)/ans(2)
  subplot(1,2,1), plot(soL,f1,col), hold on;
  subplot(1,2,2), plot(soL,f2,col), hold on;
end

%% CALCULATE \bfH MATRIX;
Hmat(1,1)   = b1coeffs.'*V1_coeffs/Area;%[b1coeffs,U1_coeffs]
Hmat(1,2)   = b1coeffs.'*V2_coeffs/Area;%[b1coeffs,U2_coeffs]
Hmat(2,1)   = b2coeffs.'*V1_coeffs/Area;%[b2coeffs,U1_coeffs]
Hmat(2,2)   = b2coeffs.'*V2_coeffs/Area;%[b2coeffs,U2_coeffs]

%  [bn_coeffs,b2coeffs-1i*b1coeffs]
%  Hmat+phi'*m_irrs*eye(2)
%  Hmat(1,1)=F1_coeffs.'*U1_coeffs/Area;
%  Hmat(1,2)=F1_coeffs.'*U2_coeffs/Area;
%  Hmat(2,1)=F2_coeffs.'*U1_coeffs/Area;
%  Hmat(2,1)=F2_coeffs.'*U2_coeffs/Area;

%% CALCULATE EMPs:
if 0
  rho_eff   = phi0+phi'*d_irrs;
  C55eff    = 1+Hmat(1,1)+phi'*(m_irrs-1);
  C44eff    = 1+Hmat(2,2)+phi'*(m_irrs-1);
  C45eff    = (Hmat(1,2)+Hmat(2,1))/2;
  C_eff     = [C55eff,C45eff,C44eff,rho_eff];
else
  rho_eff   = phi0+phi'*d_irrs;
  m_av      = phi0+phi'*m_irrs;
  Mmat      = m_av*eye(2)+.5*(Hmat+Hmat');
  C_eff     = {Mmat,rho_eff};
end

if 0%%energy conservation test:
  Etest  = imag(bn_coeffs'*UU_coeffs(1:Ntot))
end
