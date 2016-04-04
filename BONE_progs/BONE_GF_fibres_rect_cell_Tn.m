function C_eff=BONE_GF_fibres_rect_cell_Tn(AB,Irr_vars,Nterms);
%% CALL: C_eff=BONE_GF_fibres_rect_cell_Tn(AB,Irr_vars,Nterms);
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
  BC     = Irr_vars;
  %%
  AB  = .5*[1 1];
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
%    Irr_vars  = {irr_vars};
%  Nterms   = 5;
end

%  GFxn     = @GRN_laplace_doubly_periodic;
%  GF_args  = {AB(1),AB(2),[],[]};
GF = {@GRN_laplace_doubly_periodic,AB(1),AB(2)};

%% GET QUADRATURE POINTS AND INNER PRODUCT MATRICES:
%%  ip_stuff={ {ip_d1chi,ip_chi},hn=[2,1,1...,1]',...
%%                {d1chi_vals,chi_vals} };
Nint  = 300;
Nunc  = Nterms;

[soL,ip_stuff,xtra_ints,Mlog_in,wq_in] = ...
   BONE_ipmatrices_Tn(Nint,Nterms);
ip_d1chi = ip_stuff{1}{1}(2:end,:);
ip_chi   = ip_stuff{1}{2}(2:end,:);
ip1_Tn   = xtra_ints{1};
%%
[ip_Pn,iplv_0diff_Pn_oL,iplv_1diff_Pn,xtraP] = ...
   BONE_Pn_Stuff(soL,ip1_Tn,Nterms);
Pn_vals  = xtraP{1};
Pn_ends  = xtraP{2};
%%
nvec     = (1:Nterms)';
hnT      = ip_stuff{2}{1}(nvec+1);%% T_1 -> T_N
Tn_vals  = ip_stuff{3}{1}(:,nvec+1)*diag(hnT);
hnU      = ip_stuff{2}{2}(nvec);%% U_0->U_{N-1}
%  {ip_stuff{3}{2}(:,nvec),diag(-nvec.*hnU)}
Un_vals  = ip_stuff{3}{2}(:,nvec+1)*diag(-nvec.*hnU);
dTn_vals = Un_vals*diag(nvec);

%  [soL,ip_stuff,Mlog,w_qd]=...
%     BONE_ipmatrices_cos(Nint,Nterms);
%  IP_stuff{1,length(ip_stuff)+1}=Mlog;
%  IP_stuff(1:end-1)=ip_stuff;
%  %%
%  kn=pi*(1:Nterms)';
%  DD=diag(1./kn);
%  TTf2c=[0*DD,-DD;DD,0*DD];%%fib-cav map
%  DD=diag(kn);
%  MD_inr=[0*DD,DD;-DD,0*DD];%% differentiates F.series
%  TTc2f=MD_inr;%%cav->fib map

%% GET (x,y) AND \theta EVALUATED AT THE QUAD POINTS:
%  nvec=(1:Nterms)';
%  nvec2=(1:2*Nterms)';
%  ip_d1chi=ip_stuff{1}{1}(2:end,:);
%  ip_chi=ip_stuff{1}{2}(2:end,:);
%%
Nirregs     = size(Irr_vars,1);
Ntot        = Nirregs*Nunc;
bn_coeffs   = zeros(Ntot,1);
cn_coeffs   = zeros(Ntot,1);

for j=1:Nirregs
   irr_vars  = Irr_vars(j,1:3);
   [xyvecs,ds_dt,th_vec, dtheta_ds,...
     d2s_dt2,d2theta_ds2,d2xy_ds2,LC,...
       tvec,Area_fibs(j,1)] = ...
         BONE_get_rsdtheta_NRquick( irr_vars,soL );
   %%
   mdBC         = Irr_vars{j,4};
   m_irrs(j,1)  = mdBC{1};%% internal shear modulus (rel to host);
   d_irrs(j,1)  = mdBC{2};%% internal density (rel to host);
   BC(j,1)      = mdBC{3};%% boundary condition;
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
   JJ  = (1:Nunc)+(j-1)*Nunc;
   bn  = ip_d1chi*xyvecs'*[1;1i];
   cn  = ip_Pn*xyvecs'*[1;1i];

   if BC(j)==0%% FIBRE CONDITION
      bn_coeffs(JJ) = 0;
      cn_coeffs(JJ) = (1-m_irrs(j))*cn;
%       vU_b00(JJ)=(1-m_irrs(j))*imag(bn);
%       vV_b00(JJ)=TTc2f*real(bn);

   else%% CAVITY/CRACK CONDITION
      bn_coeffs(JJ) = 2*m_irrs(j)/(1+m_irrs(j))*bn;
      cn_coeffs(JJ) = (1-m_irrs(j))/(1+m_irrs(j))*cn;
%       vU_b00(JJ)=imag(bn);
%       vV_b00(JJ)=TTc2f*(1-m_irrs(j))*real(bn);
   end
   %%
   Geom_Stuff(j,:) = ...
      {soL,Xvec{j},Yvec{j},LC_irrs(j),th_vec,d2xy_ds2};
   %%
   soL_ends  = [-1;1];
   [xyvecs,ds_dt,th_vec, dtheta_ds,...
     d2s_dt2,d2theta_ds2,d2xy_ds2]   = ...
         BONE_get_rsdtheta_NRquick( irr_vars,soL_ends );
   for r=1:2
      geom_stuff_inner_ends{j,r}   = ...
         {soL_ends(r),xyvecs(1,r),xyvecs(2,r),LC_irrs(j),...
            th_vec(r),d2xy_ds2(r)};
   end
end

bb = [imag(-bn_coeffs),real(bn_coeffs)]';
b1 = bb(1,:);
cc = [imag(-cn_coeffs),real(cn_coeffs)]';
c1 = cc(1,:);
%%
Area        = 4*prod(AB);
phi         = Area_fibs/Area;
phi0        = 1-sum(phi);%%volume fraction of host material;
b00_const   = 1i*(1-phi0-m_irrs.'*phi);

if 0%%plot the irregularities;
   A   = AB(1);
   B   = AB(2);

   %%plot the unit cell:
   xy  = [-B,-A;-B,A;B,A;B,-A;-B,-A];
   plot(xy(:,1),xy(:,2),'k'), hold on;

   %%plot the fibres:
   col = {'k','r'};
   for j=1:Nirregs
      plot(Xvec{j},Yvec{j},col{BC(j)+1});
   end
   return;
end


%% CALCULATE KERNEL MATRIX:
MU       = zeros(Ntot,Ntot);
MV       = MU;
FF       = MU(:,1);
vU_b00   = FF;
vV_b00   = FF;
%%
iplv_0diff  = LC*iplv_0diff_Pn_oL;%%to expand in terms of P_n/(h_n*L)
iplv_1diff  = iplv_1diff_Pn;%%to expand in terms of P_n/(h_n*L)
OP_ends_v   = Pn_ends;
OP_ends_u   = [(-1).^nvec,1+0*nvec];
%%
for j=1:Nirregs
   JJ           = (1:Nterms)+(j-1)*Nunc;
   LC           = LC_irrs(j);
   mj           = m_irrs(j);
   geom_stuff   = Geom_Stuff(j,:);
   if BC(j)==0%%FIBRE CONDITION:
      if 1%% multiply by T_n to make (1/2)*(1+m_j)[W_n] term in IE
          %% produce a diagonal matrix;
         iplu_0diff = LC*Tn_vals'*diag(ip1_Tn);
            %% multiply by T_m(s/L) on the left;
         iplu_1diff = -dTn_vals'*diag(ip1_Tn);
            %% multiply by T'_m(s/L) on the left;
         MU(JJ,JJ)  = (mj+1)/2*eye(Nterms);
            %% diagonal comes from (1/2)*(1+m_j)[W_n] term in IE for [W_n];

      else%% multiply by P_n since forcing term in IE is smooth,
          %% so integrals and (1/2)*(1+m_j)[W_n] should combine to be smooth;

         iplu_0diff = LC*Tn_vals'*diag(ip1_Tn);
           %%multiply by T_m(s/L) on the left;
         iplu_1diff = -dTn_vals'*diag(ip1_Tn);
           %%multiply by T'_m(s/L) on the left;
         MU(JJ,JJ)  = (mj+1)/2*eye(Nterms);
           %%diagonal comes from (1/2)*(1+m_j)[W_n] term in IE for [W_n];
      end
      Mlog  = iplv_1diff*Tn_vals*Mlog_in(2:end,2:end);
         %% for int_{-L}^LW.[f=d/ds(\xi1+i*xi2)]ds
         %%     =-int_{-L}^L[W_s(\xi1+i*xi2)]ds
         %% to expand W_s=a_nP_n(s/L), multiply it by iplv_1diff;

      U_FACTOR   = 1-mj;
      V_FACTOR   = 1;
      theta      = Xtra{j,1};
      %%
      FF(JJ)     = iplu_0diff*U_FACTOR*exp(1i*theta);
      vU_b00(JJ) = imag(FF(JJ));
      vV_b00(JJ) = iplv_0diff*V_FACTOR*cos(theta);
   else%%CRACK/CAVITY CONDITION:
      iplu_1diff = ip_d1chi;%%multiply by weighted T_m's;
      iplu_0diff = -LC*ip_chi;%%multiply by weighted U_m's;
      Mlog       = Mlog_in(2:end,2:end);%%diagonal matrix - inner prod with weighted T_m's twice;
      Mdiag0     = ip_d1chi*Pn_vals;
      MV(JJ,JJ)  = (1-mj)/2*Mdiag0';%%??
      %%
      U_FACTOR   = 1;
      V_FACTOR   = 1+mj;%%??
      theta      = Xtra{j,1};
      %%
      FF(JJ)     = iplu_1diff*U_FACTOR*[Xvec{j},Yvec{j}]*[1;1i];
      vU_b00(JJ) = iplu_0diff*U_FACTOR*sin(theta);
      vV_b00(JJ) = iplv_0diff*V_FACTOR*cos(theta);
   end
   %%
   for r=1:Nirregs
      geom0_stuff   = Geom_Stuff(r,:);
      GEOM_stuff    = {geom_stuff,geom0_stuff};
%   T Tf2c*mGn_U*dnW_inr,pause
      %%
      JJ0  = (1:Nunc)+(r-1)*Nunc;
      %%->(m_j-1)*\pa_n\int_inr.G[\pa_n'W].ds'

      %%???????????????????????????????????????????????????????
      if BC(j)==0%%why no dependence on BC(r)??
      %%???????????????????????????????????????????????????????

         IP_stuff = {iplu_0diff,ip_d1chi,[]};
         mGn_U    = BONE_kernelmat_dn_arbGF(GEOM_stuff,IP_stuff,GF);
         %%
         if j==r
            IP_stuff  = {iplv_1diff,ip_d1chi,Mlog};
         else
            IP_stuff  = {iplv_1diff,ip_d1chi,[]};
         end
         IS_CLOSED   = 0;
         mGs_U       = BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,...
                          IS_CLOSED,GF);
             %% want to expand int_{-L}^LG_s[W_n]ds = sum_n[a_nP_n];
             %%  => expand int_{-L}^LG[W_n]ds = sum_n[a'_nT_n];
             %% ip1v_1diff*[a'_n]=[a_n]
         %%
         MU(JJ,JJ0)  = ...
               MU(JJ,JJ0)-U_FACTOR*mGn_U;
         MV(JJ,JJ0)  = -mGs_U;
         %%
         Mlog2          = ( diag(Mlog_in(2:end,2:end)).*OP_ends_u(:,2) )';
         IP_stuff       = {1,ip_d1chi,Mlog2};
         GEOM_stuff_end = ...
            {geom_stuff_inner_ends{j,2},geom0_stuff};
         Gp1   = BONE_kernelmat_arbGF(GEOM_stuff_end,IP_stuff,...
                    IS_CLOSED,GF);
         %%
         Mlog1          = ( diag(Mlog_in(2:end,2:end)).*OP_ends_u(:,1) )';
         IP_stuff       = {1,ip_d1chi,Mlog1};
         GEOM_stuff_end = ...
            {geom_stuff_inner_ends{j,1},geom0_stuff};
         Gm1   = BONE_kernelmat_arbGF(GEOM_stuff_end,IP_stuff,...
                    IS_CLOSED,GF);
         %%
         for it=1:Nterms
            mvGp1(it,:)  = Gp1;
            mvGm1(it,:)  = Gm1;
         end
         MV(JJ,JJ0)  = MV(JJ,JJ0)+...
            -V_FACTOR*( diag( OP_ends_v(:,2) )*mvGp1 +...
                        -diag( OP_ends_v(:,1) )*mvGm1 );

      else%% if BC==1
         IP_stuff = {iplv_0diff,ip_d1chi,[]};
         mGn_U    = BONE_kernelmat_dn_arbGF(GEOM_stuff,...
                       IP_stuff,GF);
         %%
         if j==r
            IP_stuff  = {iplu_1diff,ip_d1chi,Mlog};
%             IP_stuff  = {ip_d1chi,ip_d1chi,Mlog};
         else
            IP_stuff  = {iplu_1diff,ip_d1chi,[]};
         end

         IS_CLOSED   = 0;
         mGs_U       = BONE_kernelmat_arbGF(GEOM_stuff,IP_stuff,...
                          IS_CLOSED,GF);
         %%
         MU(JJ,JJ0)  = mGs_U;
         MV(JJ,JJ0)  = ...
               MV(JJ,JJ0)-V_FACTOR*mGn_U;%%%%%
      end
   end
end

mb00  = -1/Area*(b1+c1*MV);
mb00  = [mb00,1-1/Area*c1*vV_b00];
%%
MU = [ [MU,vU_b00];mb00 ];
MV = [MV,vV_b00];
FF = [FF;b00_const];


%% SOLVE INTEGRAL EQUATION:
UU_coeffs      = MU\FF;
VV_coeffs      = MV*UU_coeffs;
b00            = UU_coeffs(end);
UU_coeffs(end) = [];
%%
H0    = bb*UU_coeffs+cc*VV_coeffs;
Hmat  = real([1i*H0,H0])/Area;
%%
rho_eff  = phi0+phi'*d_irrs;
m_av     = phi0+phi'*m_irrs;
Mmat     = m_av*eye(2)+.5*(Hmat+Hmat');
C_eff    = {Mmat,rho_eff};


if 0%%energy conservation test:
  Etest  = imag(bn_coeffs'*UU_coeffs(1:Ntot))
end
