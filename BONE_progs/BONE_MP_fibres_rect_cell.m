function C_eff=BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
%% CALL: C_eff=BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
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
  Nterms = AB;
  col    = Irr_vars;
  %%
  AB  = .5*[1 1];
  if 1
    crk_fxn    = @CURVEprof_circarc;
    crk_prams  = {1};%% fraction of circle
    radius     = .225
    srt        = {radius*[1 1],0,[0 0]};%area=pi*radius^2
    midi       = [10 1]
  end
  Irr_vars  = {crk_fxn,crk_prams,srt,midi};
%    Irr_vars={irr_vars};
%  Nterms=5;
end

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint           = 300;
Nunc           = 2*Nterms;
[soL,ip_stuff] = ...
   BONE_ipmatrices_cos(Nint,Nterms);

%% GET (x,y) AND \theta EVALUATED AT THE QUAD POINTS:
nvec     = (1:Nterms)';
nvec2    = (1:2*Nterms)';
ip_d1chi = ip_stuff{1}{1}(2:end,:);
w_qd     = 2*ip_stuff{1}{1}(1,:);
%%
Nirregs     = size(Irr_vars,1);
Ntot        = Nirregs*Nunc;
bn_coeffs   = zeros(Ntot,1);

for j=1:Nirregs
  irr_vars     = Irr_vars(j,1:3);
  centre(j,:)  = irr_vars{3}{3};
  %%
  midi         = Irr_vars{j,4};
  m_irrs(j,1)  = midi(1);
  d_irrs(j,1)  = midi(2);
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,Area_cavs(j,1)] = ...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  if Area_cavs(j)>0
    %% NB normal should point INTO the cavity:
    xyvecs        = fliplr(xyvecs);
    ds_dt         = flipud(ds_dt);
    th_vec        = flipud(th_vec)+pi;
    dtheta_ds     = -flipud(dtheta_ds);
    d2s_dt2       = -flipud(d2s_dt2);
    d2theta_ds2   = flipud(d2theta_ds2);
    d2xy_ds2      = flipud(d2xy_ds2);
    tvec          =  flipud(tvec);
  end
  Xvec{j}      = xyvecs(1,:)';
  Yvec{j}      = xyvecs(2,:)';%[Xvec{j},Yvec{j}],pause
  LC_irrs(j)   = LC;
  Xtra(j,:)    = {th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  %%
  JJ              = nvec2+(j-1)*Nunc;
  bn_coeffs(JJ)   = ip_d1chi*(1-midi(1))*exp(1i*th_vec);
end

if 0
  DD           = diag(nvec*pi/LC);
  MD           = [0*DD,DD;-DD,0*DD];
  MDinv        = inv(MD);
  bnV_coeffs   = (1-midi(1))*ip_d1chi*xyvecs'*[1;1i];
  b1Vcoeffs    = -imag(bnV_coeffs);
  b2Vcoeffs    = real(bnV_coeffs);
  %%
%    [bn_coeffs,MD*bnV_coeffs]
%    [imag(bn_coeffs),-MD*b1Vcoeffs]
%     [-real(bn_coeffs),-MD*b2Vcoeffs]
end

b1coeffs    = imag(bn_coeffs);
b2coeffs    = -real(bn_coeffs);
Area        = 4*prod(AB);
Area0       = Area-sum(Area_cavs);
phi0        = Area0/Area;%%volume fraction of host material;
phi         = Area_cavs/Area;
b00_const   = -1i*(1-phi0);

if 0%%plot the cavities:
  A   = AB(1);
  B   = AB(2);
  xy  = [-B,-A;-B,A;B,A;B,-A;-B,-A];
  plot(xy(:,1),xy(:,2),'k'), hold on;
  %%
  for j=1:Nirregs
    plot(Xvec{j},Yvec{j},'r');
    TH   = Xtra{j,1};
    for r=1:5:length(TH)
      x  = Xvec{j}(r)+.05*[cos(TH(r));0;sin(TH(r))];
      y  = Yvec{j}(r)+.05*[sin(TH(r));0;-cos(TH(r))];
      plot(x,y,'g');
    end
  end
  daspect([1 1 1]);
  pause;
  return;
end

%% CALCULATE KERNEL MATRIX:
MK_dnWout   = zeros(Ntot,Ntot);
MK_Wout     = MK_dnWout;
MK_Win      = MK_dnWout;
MK_dnWin    = MK_dnWout;
MP_stuff    = {[],AB};

for j=1:Nirregs
  for r=1:Nirregs
    geom_stuff = {soL,Xvec{j},Yvec{j},...
      LC_irrs(j),Xtra(j,:)};
    %%
    JJ   = nvec2+(j-1)*Nunc;
    JJ0  = nvec2+(r-1)*Nunc;
    %%
    MP_stuff{1}                  = centre(r,:);
    [M_dnWout,M_Wout,xtra_out]   = ...
       BONE_MP_kernelmat_dbl_periodic(...
           geom_stuff,ip_stuff,MP_stuff,Nterms);
    MK_dnWout(JJ,JJ0)            = M_dnWout;
    MK_Wout(JJ,JJ0)              = LC_irrs(j)*M_Wout;
  end
  [M_dnWin,M_Win,xtra_in]  = BONE_MP_kernelmat(...
            geom_stuff,ip_stuff,centre(j,:),Nterms);
  MK_dnWin(JJ,JJ) = m_irrs(j)*M_dnWin;
  MK_Win(JJ,JJ)   = LC_irrs(j)*M_Win;
end

%% SOLVE INTEGRAL EQUATION:
FF       = [bn_coeffs;0*bn_coeffs];
MK       = [MK_dnWout,-MK_dnWin;MK_Wout,-MK_Win];
W_coeffs = MK\FF;
%%
JJout       = 1:Ntot;
JJin        = JJout+Ntot;
UU_coeffs   = MK_Win*W_coeffs(JJin);
%  UU_coeffs   = MK_Wout*W_coeffs(JJ);

if 0
  bnV_coeffs   = (1-midi(1))*ip_d1chi*xyvecs'*[1;1i];
  b1Vcoeffs    = -imag(bnV_coeffs);
  b2Vcoeffs    = real(bnV_coeffs);
  %%
  dsW          = xtra_in{3}*W_coeffs(JJin);
  VV_coeffs_MP = LC*ip_d1chi*dsW;
  %%
  dnW_in       = xtra_in{2}*W_coeffs(JJin);
  dnW_out      = xtra_out{2}*W_coeffs(JJout);
  UU_coeffs_MP = LC*ip_d1chi*(dnW_out-dnW_in);
  %%
  Wcoeffs_MP_in   = W_coeffs(JJin);
  Wcoeffs_MP_out  = W_coeffs(JJout);
  %%
  save MPfib Wcoeffs_MP_in Wcoeffs_MP_out VV_coeffs_MP UU_coeffs_MP b1Vcoeffs b2Vcoeffs
end

if 0
  W1     = xtra_out{1}*W_coeffs(JJout);
  wn1    = ip_d1chi*W1;[wn1,UU_coeffs];%,LC,pause
  CS     = ip_stuff{3}{1}(:,2:end);
  W1_ap  = CS*wn1;
  W2     = xtra_in{1}*W_coeffs(JJin);%[W1,W2],pause
  %%
%    subplot(3,2,1), plot(soL,real(W1));% W_2
%    hold on, plot(soL,real(W2),'g');
%    plot(soL,real(W1_ap),'--r'), hold off;
%    subplot(3,2,2), plot(soL,-imag(W1));%% W_1
%    hold on, plot(soL,-imag(W2),'g');
%    plot(soL,-imag(W1_ap),'--r'), hold off;
  %%
  W_s1   = xtra_out{3}*W_coeffs(JJout);
  W_s2   = xtra_in{3}*W_coeffs(JJin);
  tstWs  = LC_irrs(1)*ip_d1chi*[W_s1,W_s2];%%
  lc     = LC_irrs(1)
  kn     = (1:Nterms)'*pi/lc;
  DD     = diag(kn);
  MD     = [0*DD,DD;-DD,0*DD];
  dwn    = MD*wn1;%%!!INNER dsW coeffs
  dW_ap  = CS*dwn;
%    subplot(3,2,3), 
plot(soL,real(W_s1)), hold on;
%    hold on, plot(soL,real(dW_ap),'--r'), hold off;
%    subplot(3,2,4), plot(soL,-imag(W_s));
%    hold on, plot(soL,-imag(dW_ap),'--r'), hold off;
  %%
  W_n1   = xtra_out{2}*W_coeffs(JJout);
  W_n2   = xtra_in{2}*W_coeffs(JJin);
  Y      = W_n1-midi(1)*W_n2;[Y,(1-midi(1))*exp(1i*th_vec)]

  [dxW_p,dyW_p]   = ...
     BONE_diff_multipoles_laplace_doubly_periodic(...
             Xvec{1},Yvec{1},Nterms,AB,[0 0]);
   ct    = cos(th_vec);%th_vec/pi,pause
   st    = sin(th_vec);
   dsW_p = (diag(ct)*dxW_p+diag(st)*dyW_p);%dsW_p(1:10,1:4),pause
   dsW_p = ...
     LC_irrs(1)*ip_d1chi*[real(dsW_p),imag(dsW_p)]*W_coeffs(JJout);
   %[tstWs,dsW_p]


%    subplot(3,2,5), plot(soL,real(Y));
%    hold on, plot(soL,(1-midi(1))*cos(th_vec),'--r'), hold off;
%    subplot(3,2,6), plot(soL,-imag(Y));
%    hold on, plot(soL,-(1-midi(1))*sin(th_vec),'--r'), hold off;
  %%
%    ig=sin(th_vec).*(-imag(W2));
%    mx=1+(midi(1)-1)/Area*(lc*w_qd)*ig
end

%% CALCULATE \bfH MATRIX;
if 0
  U1_coeffs = -imag(UU_coeffs);
  U2_coeffs = real(UU_coeffs);
  %%
  Hmat(1,1) = b1coeffs.'*U1_coeffs/Area;
  Hmat(1,2) = b1coeffs.'*U2_coeffs/Area;
  Hmat(2,1) = b2coeffs.'*U1_coeffs/Area;
  Hmat(2,2) = b2coeffs.'*U2_coeffs/Area;
else
  H0     = [b1coeffs,b2coeffs]'/Area*UU_coeffs;
  Hmat   = real([1i*H0,H0]);
end

if 0
%    [VV_coeffs_MP,MD*UU_coeffs]
  H11 = b1Vcoeffs.'*imag(-VV_coeffs_MP)/Area
  H12 = b1Vcoeffs.'*real(VV_coeffs_MP)/Area
  H21 = b2Vcoeffs.'*imag(-VV_coeffs_MP)/Area
  H22 = b2Vcoeffs.'*real(VV_coeffs_MP)/Area
end

%% CALCULATE EMPs:
if 0
  rho_eff   = phi0+phi'*d_irrs;
  C44eff    = 1+Hmat(1,1)+phi'*(m_irrs-1);
  C55eff    = 1+Hmat(2,2)+phi'*(m_irrs-1);
  C45eff    = (Hmat(1,2)+Hmat(2,1))/2;
  C_eff     = [C44eff,C45eff,C55eff,rho_eff];
else
  rho_eff   = phi0+phi'*d_irrs;
  m_av      = phi0+phi'*m_irrs;
  Mmat      = m_av*eye(2)+.5*(Hmat+Hmat');
  C_eff     = {Mmat,rho_eff};
end
if 0%%energy conservation test:
  Etest  = imag(bn_coeffs'*UU_coeffs(1:Ntot))
end

if 0%% generate some outputs for GF prog:
  Wcoeffs_MP   = UU_coeffs;
  %%
  yy     = B*soL;
  xx     = A+0*yy;
  centre = [0 0];
  W_x    = BONE_diff_multipoles_laplace_doubly_periodic(...
             xx,yy,Nterms,AB,centre)*W_coeffs(JJ_out);
  b00_MP = (1/2/B)*w_qd*W_x;
  %%
  Ucoeffs_MP   = MK_dnWout*W_coeffs(JJout)+...
             -MK_dnWin*W_coeffs(JJin)/m_irrs(1);
  save MPfibres b00_MP U_coeffs_MP W_coeffs_MP
end
