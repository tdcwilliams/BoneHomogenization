function C_eff=BONE_cracks_rect_cell(AB,Irr_vars,Npolys);
%% CALL: C_eff=BONE_cracks_rect_cell(AB,Irr_vars,Npolys);
%%
%% DESCRIPTION:
%% Works out the effective material properties for a periodic elastic
%%  medium (rectangular unit cell)
%%  when there are N cracks per unit cell;
%% Solves the first order cell problem to get a homogenized zero-th
%%  order wave eqn;
%% Allows for singularities at the crack end-points by using
%%  weighted Chebyshev polynomials in a Galerkin scheme;
%%
%% INPUTS:
%% AB=[A,B] gives the non-dim sizes of the unit cell
%%  (A is the height, B the width);
%% Irr_vars{j}={function handle, curve parameters,...
%%              [scalers, rotation angle, translation]}
%%   for the j-th crack;
%% Npolys = order of chebyshev approximation for the jump
%%  in displacement across the crack;
%%
%% OUTPUTS:
%% Effective material properties
%%  C_eff=[C^*_44,C^*_45,C^*_55]

%  om=physvars(1);
%  mu=physvars(2);
%  rho=physvars(3);
%  k0sq=rho*om^2/mu;
%  k0=sqrt(k0sq);
%  Lnd=1/k0;
%  ab=ab_dim/Lnd;

if nargin==0%%use some test inputs:
  Npolys=10;
  AB=[.5 .5];
  if 0
    crk_fxn=@CURVEprof_circarc;
    crk_prams={.5};%% fraction of circle
    srt={.1*[1 1],0,[0 0]};
    Irr_vars={crk_fxn,crk_prams,srt};
  else
    a=.25;
    Irr_vars=...
      CURVEget_strtline(-a,a);
    A=AB(1);
    B=AB(2);
    Ceff_test=[1,0,1/(1-2*A/pi/B*log(cos(pi*a/2/A))),1]
  end
end

GFxn=@GRN_laplace_doubly_periodic;
GF_args={AB(1),AB(2),[],[]};

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
Nint=300;
Nterms=Npolys;
Nunc=Nterms;
[soL_Tn,ip_stuff,xtra_ints,Mlog,w_Tn]=...
   BONE_ipmatrices_Tn(Nint,Nterms);
IP_stuff{1,length(ip_stuff)+1}=Mlog;
IP_stuff(1:end-1)=ip_stuff;
%%
[soL_Pn,w_Pn]=OP_numint_legendre(Nint);

%% GET (x,y) AND \theta EVALUATED AT THE QUAD POINTS:
nvec=(1:Nterms)';
ip_d1chi=ip_stuff{1}{1}(2:end,:);
ip_chi=ip_stuff{1}{2}(2:end,:);
%%
Nirregs=size(Irr_vars,1);
Ntot=Nirregs*Nunc;
bn_coeffs=zeros(Ntot,1);
%%
for j=1:length(Nirregs)
  irr_vars=Irr_vars(j,:);
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC]=...
      BONE_get_rsdtheta_NRquick( irr_vars,soL_Tn );
  Xvec_Tn{j}=xyvecs(1,:)';
  Yvec_Tn{j}=xyvecs(2,:)';
  LC_irrs(j)=LC;
  Xtra_Tn(j,:)={th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ=nvec+(j-1)*Nunc;
  bn_coeffs(JJ)=-ip_d1chi*xyvecs'*[1; 1i];
    %% Chebyshev coefficients of [x,-y], to compute H;
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC]=...
      BONE_get_rsdtheta_NRquick( irr_vars,soL_Pn );
  Xvec_Pn{j}=xyvecs(1,:)';
  Yvec_Pn{j}=xyvecs(2,:)';
  Xtra_Pn(j,:)={th_vec,dtheta_ds,ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
end
b1_coeffs=imag(bn_coeffs);%% -\xi_2; -\pa_s->n_1=\sin(\th)
b2_coeffs=-real(bn_coeffs);%% \xi_1; -\pa_s->n_2=-\cos(\th)
FF=[-bn_coeffs;0];

%% CALCULATE KERNEL MATRIX:
MK=zeros(Ntot+1,Ntot+1);
MK(Ntot+1,Ntot+1)=1;
nvec=(1:Nterms)';%%don't need the n=0 term:
for j=1:Nirregs
  for r=1:Nirregs
    KEEPSING=(j~=r);
    GF_args{4}=KEEPSING;
    %%
    geom_stuff={soL_Tn,Xvec_Tn{j},Yvec_Tn{j},LC_irrs(j),Xtra_Tn(j,:)};
    geom0_stuff={soL_Tn,Xvec_Tn{r},Yvec_Tn{r},LC_irrs(r),Xtra_Tn(r,:)};
    GEOM_stuff={geom_stuff,geom0_stuff};
    %%
    JJ=nvec+(j-1)*Nunc;
    JJ0=nvec+(r-1)*Nunc;
    MK(JJ,JJ0)=BONE_kernelmat_crk_Tn(GEOM_stuff,...
                 ip_stuff,GFxn,GF_args);
  end
end

Area=4*prod(AB);
MK(1:end-1,Ntot+1)=-b1_coeffs;
MK(Ntot+1,1:end-1)=-b1_coeffs/Area;

%% SOLVE INTEGRAL EQUATION:
UU_coeffs=MK\FF;
b00=UU_coeffs(end);
UU_coeffs(end)=[];
%%
U1_coeffs=-imag(UU_coeffs);
U2_coeffs=real(UU_coeffs);

%% CALCULATE \bfH MATRIX;
Hmat(1,1)=b1_coeffs.'*U1_coeffs/Area;
Hmat(1,2)=b1_coeffs.'*U2_coeffs/Area;
Hmat(2,1)=b2_coeffs.'*U1_coeffs/Area;
Hmat(2,2)=b2_coeffs.'*U2_coeffs/Area;
%  Hmat(1,1)=F1_coeffs.'*U1_coeffs/Area;
%  Hmat(1,2)=F1_coeffs.'*U2_coeffs/Area;
%  Hmat(2,1)=F2_coeffs.'*U1_coeffs/Area;
%  Hmat(2,1)=F2_coeffs.'*U2_coeffs/Area;

%% CALCULATE EMPs:
rho_eff=1;
C44_eff=1+Hmat(1,1);
C55_eff=1+Hmat(2,2);
C45_eff=(Hmat(1,2)+Hmat(2,1))/2;
C_eff=[C44_eff,C45_eff,C55_eff,rho_eff];