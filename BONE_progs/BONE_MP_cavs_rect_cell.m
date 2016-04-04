function C_eff=BONE_MP_cavs_rect_cell(AB,Irr_vars,Nterms);
%% CALL: C_eff=BONE_MP_cavs_rect_cell(AB,Irr_vars,Nterms);
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
%% Irr_vars(j,:)={function handle, curve parameters,...
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
    crk_prams={[1 1]};%% fraction of circle
    radius=.225
    srt={radius*[1 1],0,[0 0]};%area=pi*radius^2
  end
  Irr_vars={crk_fxn,crk_prams,srt};
%    Irr_vars={irr_vars};
%  Nterms=5;
end

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint=150;
Nunc=2*Nterms;
[soL,ip_stuff]=...
   BONE_ipmatrices_cos(Nint,Nterms);

%% GET (x,y) AND \theta EVALUATED AT THE QUAD POINTS:
nvec=(1:Nterms)';
nvec2=(1:2*Nterms)';
ip_d1chi=ip_stuff{1}{1}(2:end,:);
%%
Nirregs=size(Irr_vars,1);
Ntot=Nirregs*Nunc;
bn_coeffs=zeros(Ntot,1);
FAC=bn_coeffs;

for j=1:Nirregs
  irr_vars=Irr_vars(j,:);
  centre(j,:)=irr_vars{3}{3};
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,Area_j]=...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  if Area_j>0
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
  end
  Xvec{j}=xyvecs(1,:)';
  Yvec{j}=xyvecs(2,:)';
  Area_cavs(j,1)=abs(Area_j);
%    fac=-sign(Area_j);
  LC_irrs(j)=LC;%plot(th_vec),pause
  Xtra(j,:)={th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ=nvec2+(j-1)*Nunc;
  bn_coeffs(JJ)=ip_d1chi*exp(1i*th_vec);
%    FAC(JJ)=1;%fac;
end

b1coeffs=imag(bn_coeffs);%%\sin(\th)=n1
b2coeffs=-real(bn_coeffs);%%-\cos(\th)=n2
%% NB normal points INTO the cavity
Area=4*prod(AB);
Area0=Area-sum(Area_cavs);
phi0=Area0/Area;%%volume fraction of host material;
b00_const=-1i*(1-phi0);

if 0%%plot the cavities:
  A=AB(1);
  B=AB(2);
  xy=[-B,-A;-B,A;B,A;B,-A;-B,-A];
  plot(xy(:,1),xy(:,2),'k'), hold on;
  %%
  for j=1:Nirregs
    plot(Xvec{j},Yvec{j},'k'), hold on;
    for r=1:15:Nint
      ll=.1;
      u1=cos(th_vec(r));
      u2=sin(th_vec(r));
      x1=Xvec{j}(r);
      y1=Yvec{j}(r);
      plot([x1,x1+ll*u1],[y1,y1+ll*u2],'g');
      plot([x1,x1+ll*u2],[y1,y1-ll*u1],'r');
    end
    daspect([1 1 1]);
  end
  return;
end


%% CALCULATE KERNEL MATRIX:
MK_dnW=zeros(Ntot,Ntot);
MK_W=MK_dnW;
MP_stuff={[],AB};

for j=1:Nirregs
  for r=1:Nirregs
    geom_stuff={soL,Xvec{j},Yvec{j},...
      LC_irrs(j),Xtra(j,:)};
    %%
    JJ=nvec2+(j-1)*Nunc;
    JJ0=nvec2+(r-1)*Nunc;
    %%
    MP_stuff{1}=centre(r,:);
    [M_dnW,M_W,xtra]=BONE_MP_kernelmat_dbl_periodic(...
              geom_stuff,ip_stuff,MP_stuff,Nterms);
    MK_dnW(JJ,JJ0)=M_dnW;
    MK_W(JJ,JJ0)=LC_irrs(j)*M_W;
  end
end

%% SOLVE INTEGRAL EQUATION:
FF=bn_coeffs;
W_coeffs=MK_dnW\FF;
UU_coeffs=MK_W*W_coeffs;
U1_coeffs=-imag(UU_coeffs(1:Ntot));
U2_coeffs=real(UU_coeffs(1:Ntot));

if 0
  W=xtra{1}*W_coeffs;
  wn=ip_d1chi*W;
  CS=ip_stuff{3}{1}(:,2:end);
  W_ap=CS*wn;
  %%
  subplot(3,2,1), plot(soL,real(W));% W_2
  hold on, plot(soL,real(W_ap),'--r'), hold off;
  subplot(3,2,2), plot(soL,-imag(W));%% W_1
  hold on, plot(soL,-imag(W_ap),'--r'), hold off;
  %%
  W_s=xtra{3}*W_coeffs;
  lc=LC_irrs(1);
  kn=(1:Nterms)'*pi/lc;
  DD=diag(kn);
  MD=[0*DD,DD;-DD,0*DD];
  dwn=MD*wn;
  dW_ap=CS*dwn;
  subplot(3,2,3), plot(soL,real(W_s));
  hold on, plot(soL,real(dW_ap),'--r'), hold off;
  subplot(3,2,4), plot(soL,-imag(W_s));
  hold on, plot(soL,-imag(dW_ap),'--r'), hold off;
  %%
  W_n=xtra{2}*W_coeffs;
  subplot(3,2,5), plot(soL,real(W_n));
  hold on, plot(soL,cos(th_vec),'--r'), hold off;
  subplot(3,2,6), plot(soL,-imag(W_n));
  hold on, plot(soL,-sin(th_vec),'--r'), hold off;
end

%% CALCULATE \bfH MATRIX;
%%  Hmat=(1/Area)*\oint(\bfn\bfW^T)\rmd s
%%  NB \bfn is the normal pointing INTO the cavity;
Hmat(1,1)=b1coeffs.'*U1_coeffs/Area;%[b1coeffs,U1_coeffs]
Hmat(1,2)=b1coeffs.'*U2_coeffs/Area;%[bn_coeffs(:,1),U2_coeffs]
Hmat(2,1)=b2coeffs.'*U1_coeffs/Area;%[b2coeffs,U1_coeffs]
Hmat(2,2)=b2coeffs.'*U2_coeffs/Area;%[b2coeffs,U2_coeffs]
%  Hmat=Area0/Area*Hmat;

%% CALCULATE EMPs:
if 0
  rho_eff=phi0;
  C44eff=phi0+Hmat(1,1);
  C55eff=phi0+Hmat(2,2);
  C45eff=(Hmat(1,2)+Hmat(2,1))/2;
  C_eff=[C44eff,C45eff,C55eff,rho_eff];
else
  rho_eff=phi0;
  m_av=phi0;
  Mmat=m_av*eye(2)+.5*(Hmat+Hmat');
  C_eff={Mmat,rho_eff};
end
if 0%%energy conservation test:
  Etest=imag(bn_coeffs'*UU_coeffs(1:Ntot))
end

if 0%% outputs to help fix DPGF method:
  b00_theory=Hmat(1,2)+1i*( 1-phi0-Hmat(1,1) )
  return;
  %%
  A=AB(1);
  B=AB(2);
  Y=B*soL;
  X=A+0*Y;
  centre=[0 0];
  M=Nterms;
  Wx0=BONE_diff_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre);
  Wx=[real(Wx0),imag(Wx0)]*W_coeffs;
  WpA0=BONE_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre);
  WpA=[real(WpA0),imag(WpA0)]*W_coeffs;
  b00=ip_stuff{1}{1}(1,:)*Wx
  a00=ip_stuff{1}{1}(1,:)*WpA
  %%
  plot(Y,imag(Wx)), hold on;
  plot(Y,imag(b00)+0*Y);
  %%
  X=-A+0*Y;
  Wx0=BONE_diff_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre);
  Wx=[real(Wx0),imag(Wx0)]*W_coeffs;
  WmA0=BONE_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre);
  WmA=[real(WmA0),imag(WmA0)]*W_coeffs;
%    a00=ip_stuff{1}{1}(1,:)*WpA
  b00=ip_stuff{1}{1}(1,:)*Wx
  plot(Y,imag(Wx),'--r');
  plot(Y,imag(b00)+0*Y,'--r'), hold off;
  %%
  intr=2*LC*ip_stuff{1}{1}(1,:);
  Wcirc=xtra{1}*W_coeffs;
  Acirc=abs(Area_j);%[Acirc,LC^2/pi]
  tst1=[b00*A,1/(4*B)*intr*...
          ( sin(th_vec).*Wcirc-(Xvec{1}+A).*exp(1i*th_vec) )]
  tst2=[b00*A,1/(4*B)*intr*...
          ( sin(th_vec).*Wcirc-1i*Xvec{1}.*sin(th_vec) )]
  tst3=[b00,1i*(1-phi0)+1/(4*A*B)*intr*( sin(th_vec).*Wcirc )]
  tst4=[b00,1i*(1-phi0)+Hmat(1,2)-1i*Hmat(1,1)]
  tst5=[intr*( sin(th_vec).*Wcirc ),...
     (4*A*B)*(Hmat(1,2)-1i*Hmat(1,1))]
  tst6=[intr*(Xvec{1}.*sin(th_vec)),-Acirc]
end

if 0
  mp_coeffs=W_coeffs;
  B1=b1coeffs;
  B2=b2coeffs;
  U1=U1_coeffs;
  U2=U2_coeffs;
  b00_theory=Hmat(1,2)+1i*( 1-phi0-Hmat(1,1) );
  Hmat
  save MPcav mp_coeffs B1 B2 U1 U2 b00_theory Hmat
end