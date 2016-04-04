function C_eff=BONE_MPG_cracks_hexcell(BDrho,Irr_vars_outer,...
  Irr_vars_inner,N_MP,Nterms);
%% CALL: C_eff=BONE_MPG_cracks_rect_cell(BDrho,Irr_vars_outer,...
%%               Irr_vars_inner,N_MP,Nterms)
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
%% BDrho=[A,B] gives the non-dim sizes of the unit cell
%%  (A is the height, B the width);
%% Irr_vars_outer{j}={function handle, curve parameters,...
%%              [scalers, rotation angle, translation]}
%%   for the j-th subregion;
%% Irr_vars_inner{j}{r}={function handle, curve parameters,...
%%              [scalers, rotation angle, translation]}
%%   for the r-th crack in the j-th subregion;
%% Nterms = order of polynomial approximation
%%  for the displacement jump along the crack;
%% N_MP = no of multipoles to use around each subregion;
%%
%% OUTPUTS:
%% Effective material properties
%% Mmat=[C44eff,C45eff;C45eff,C55eff];
%% C_eff={Mmat,rho_eff};
%% rho_eff=1 is the effective density

if nargin==1%%use some test inputs:
   Nterms = BDrho;
elseif nargin==0
   Nterms = 10;
end
if nargin<=1
%  col   = Irr_vars;
   col   = 'k';
   %%
   BDrho = .5*[sqrt(3) 2 1];
   if 1%%inner circle
     crk_fxn   = @CURVEprof_circarc;
     crk_prams =  {.95};%% fraction of circle
     radius    = .225
     srt       = {radius*[1 1],0,[0 0]};%area=pi*radius^2
   end
   Irr_vars_inner{1} = {crk_fxn,crk_prams,srt};

   %% outer surrounding circle;
   radius            = .35;
   srt               = {radius*[1 1],0,[0 0]};
   Irr_vars_outer    = {crk_fxn,crk_prams,srt};
   %%
   N_MP  = Nterms;
 %    Irr_vars={irr_vars};
 %  Nterms=5;
end

%%unit cell stuff;
B           = BDrho(1);
D           = BDrho(2);
E           = BDrho(3)*D;
A           = D+E;
uc_verts_x  = [D E -E -D -E E]';
uc_verts_y  = [0 B B 0 -B -B]';
xy          = [uc_verts_x  uc_verts_y];
Area        = GEN_area(xy);

%if nargin==2%%use some test inputs:
%  Nterms = AB;
%  N_MP   = Nterms;
%  col    = Irr_vars_outer;
%  %%
%  AB  = .5*[1 1];
%  if 1
%    a                   = .25;
%    ang                 = .5*pi;
%    rot                 = exp(1i*ang);
%    Irr_vars_inner{1}   = ...
%      CURVEget_strtline(-rot*a,rot*a);
%    A          = AB(1);
%    B          = AB(2);
%    Ceff_test  = [1,0,1/(1-2*A/pi/B*log(cos(pi*a/2/A))),1]
%    %%
%    crk_fxn          = @CURVEprof_circarc;
%    crk_prams_out    = {[1,+1]};%% fraction of circle, anticlockwise;
%    radius_out       = .4
%    srt_out          = {radius_out*[1 1],0,[0 0]};%area=pi*radius^2
%    Irr_vars_outer   = {crk_fxn,crk_prams_out,srt_out};
%  else
%    crk_fxn          = @CURVEprof_circarc;
%    crk_prams_out    = {[1,+1]};%% fraction of circle, anticlockwise;
%    radius_out       = .4
%    srt_out          = {radius_out*[1 1],0,[0 0]};%area=pi*radius^2
%    Irr_vars_outer   = {crk_fxn,crk_prams_out,srt_out};
%    %%
%    crk_prams_in        = {[1,-1]};%% fraction of circle, clockwise;
%    radius_in           = .225
%    srt_in              = {radius_in*[1 1],0,[0 0]};%area=pi*radius^2
%    Irr_vars_inner{1}   = {crk_fxn,crk_prams_in,srt_in};
%  end
%%    Irr_vars={irr_vars};
%%  Nterms=5;
%end

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint           = 250;
[soL,ip_stuff] = ...
   BONE_ipmatrices_cos(Nint,N_MP);
ipCS           = ip_stuff{1}{1}(2:end,:);
%%
Nzones      = length(Irr_vars_inner);
Ntot_outer  = Nzones*(2*N_MP);
Ntot_inner  = zeros(Nzones,1);
Nunc        = 2*N_MP;

for j=1:Nzones
  irr_vars_inner  = Irr_vars_inner{j};
  Ntot_inner(j)   = Nterms*size(irr_vars_inner,1);
%    Ntot_inner(j)   = 2*Nterms*size(irr_vars_inner,1);
  %%
  irr_vars     = Irr_vars_outer(j,:);
  centre(j,:)  = irr_vars{3}{3};
  %%
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
      tvec,Area_cavs(j,1)] = ...
        BONE_get_rsdtheta_NRquick( irr_vars,soL );
  Xvec{j}      = xyvecs(1,:)';
  Yvec{j}      = xyvecs(2,:)';
  LC_irrs(j)   = LC;
  Xtra(j,:)={th_vec,dtheta_ds,...
    ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2};
  GEOM_stuff(j,:) = ...
    {soL,xyvecs(1,:)',xyvecs(2,:)',LC,...
      th_vec,d2xy_ds2,A,B};
end

if 0%%plot the cavities:
   xy = [xy;xy(1,:)];
   plot(xy(:,1),xy(:,2),'k'), hold on;
   %%
   for j=1:length(Xvec)
      plot(Xvec{j},Yvec{j},'r');
   end
   hold on;
end


%% CALCULATE KERNEL MATRIX:
Ntot  = Ntot_outer+sum(Ntot_inner);
%  MK_outer=zeros(Ntot,Ntot_outer);
%  MK_inner=zeros(Ntot,sum(Ntot_inner));
MK       = zeros(Ntot,Ntot);
FF       = zeros(Ntot,1);
Nj       = Nunc+Ntot_inner;
J_inner  = [];
%%
b1_coeffs   = [];
b2_coeffs   = [];
Area_cavs   = 0;
%%
MP_stuff    = {[],BDrho};

for j=1:Nzones
  geom_stuff   = GEOM_stuff(j,:);
  [M_dnW,M_dsW,Fj,MK_inner,xtra] = ...
     BONE_MPG_kernelmat_zone_cracks(...
       geom_stuff,Irr_vars_inner{j},Nterms,ipCS);
  JJ        = (1:Nj(j))+sum(Nj(1:j-1));
  JJ0_inner = (1:Ntot_inner(j))'+Nunc+sum(Nj(1:j-1));
  J_inner   = [J_inner;JJ0_inner];
  %%
  FF(JJ)             = Fj;
  MK(JJ,JJ0_inner)   = MK_inner;
  b1_coeffs          = [b1_coeffs;xtra{1}];
  b2_coeffs          = [b2_coeffs;xtra{2}];
  Area_cavs          = Area_cavs+xtra{3};
  for r=1:Nzones
    JJ0_outer        = (1:Nunc)+sum(Nj(1:r-1));
    MP_stuff{1}      = centre(r,:);
    [dnW,dsW]        = BONE_MP_kernelmatV2_hexperiodic(...
              geom_stuff,MP_stuff,N_MP);%%**
    MK(JJ,JJ0_outer) = M_dnW*dnW+M_dsW*dsW;
  end
end

%% SOLVE SYSTEM:
MP_coeffs   = MK\FF;
%Area        = 4*prod(AB);%%**
phi0        = 1-Area_cavs/Area;
%%
UU_coeffs            = MP_coeffs(J_inner);
U1_coeffs            = -imag(UU_coeffs);
U2_coeffs            = real(UU_coeffs);
MP_coeffs(J_inner)   = [];

%% CALCULATE \bfH MATRIX;
%%  Hmat=(1/Area)*\oint(\bfn\bfW^T)\rmd s
%%  NB \bfn is the normal pointing INTO the cavity;
Hmat(1,1)   = b1_coeffs.'*U1_coeffs/Area;%[bn_coeffs(:,1),U1_coeffs]
Hmat(1,2)   = b1_coeffs.'*U2_coeffs/Area;%[bn_coeffs(:,1),U2_coeffs]
Hmat(2,1)   = b2_coeffs.'*U1_coeffs/Area;%[bn_coeffs(:,2),U1_coeffs]
Hmat(2,2)   = b2_coeffs.'*U2_coeffs/Area;%[bn_coeffs(:,2),U2_coeffs]
%  Hmat     = Area0/Area*Hmat;

if 0
  U1mpg        = U1_coeffs;
  U2mpg        = U2_coeffs;
  b00_theory   = Hmat(1,2)-1i*Hmat(1,1)
  save UUmpg U1mpg U2mpg b00_theory Hmat
end

%% CALCULATE EMPs:
rho_eff  = phi0;
C44eff   = phi0+Hmat(1,1);
C55eff   = phi0+Hmat(2,2);
C45eff   = (Hmat(1,2)+Hmat(2,1))/2;

Mmat     = [C44eff,C45eff;C45eff,C55eff];
C_eff    = {Mmat,rho_eff};

if 0%%energy conservation test:
  Etest  = imag(bn_coeffs'*UU_coeffs(1:Ntot))
end
