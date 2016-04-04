%% fig_crack_lim.m
%% elliptical cavities (major axis=a)
%% -> straight crack (length 2a) as eccentricity is increased;
%%
%% also can compare Achenbach & Li (1986)
%% & and Nemat-Nasser (**) results for straight cracks;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
cav_solver  = @BONE_GF_cavs_hex_cell;
crk_solver  = @BONE_GF_cracks_hex_cell;

%% equilateral hexagon with sides of length 1;
E              = .5;%% half top side (and a radius);
B              = E*cos(pi/3);%% half height;
D              = 2*E;%% half width = radius;
BDrho          = [B D E/D];
unitcell_prams = BDrho;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np       = 50;
aoA_vec  = (1:np)'/np*.5;

%% Achenbach & Li (1986);
Ach_Li   = 1./(1-2*A/pi/B*log(cos(pi/2*aoA_vec)));

Npolys   = 10;
ecc_vec  = [.05,.1,.25];
Nvec     = [20 25 35];
%  m_y=0*Ach_Li;
for j=1:np
  a         = aoA_vec(j)*A;
  disp([num2str(j),' configurations done (out of ',num2str(np),')']);
  %%
  Irr_vars  = CURVEget_strtline(-a,a);
  C_eff     = feval(crk_solver,unitcell_prams,Irr_vars,Npolys);
  m_y(j,1)  = C_eff{1}(4);%[m_y(j),Ach_Li(j)]
  m_x(j,1)  = 1;
  %%
  for r=1:length(ecc_vec)
    crk_fxn    = @CURVEprof_circarc;
    crk_prams  = {1};%% fraction of circle
    radius     = a;
    ecc        = ecc_vec(r);
    srt        = {radius*[1 ecc],0,[0 0]};
    Irr_vars   = {crk_fxn,crk_prams,srt};
    %%
    Nterms     = Nvec(r);
    C_eff      = feval(cav_solver,unitcell_prams,Irr_vars,Nterms);
    m_y(j,r+1) = C_eff{1}(4);
    m_x(j,r+1) = C_eff{1}(1);
  end
end

subplot(1,2,1);
GEN_set_plotorder;
plot(aoA_vec,m_y,'k'), hold on;
n0 = 1;
JJ = n0:2*n0:np;
plot(aoA_vec(JJ),Ach_Li(JJ),'.');
%%
subplot(1,2,2);
GEN_set_plotorder;
plot(aoA_vec,m_x,'k')
%%
if 1
  subplot(1,2,1);
  xlim([0 .5]);
  GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
  %%
  subplot(1,2,2);
  xlim([0 .5]);
  GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
  %%
  saveas(gcf,'crack_lim.fig')
  saveas(gcf,'crack_lim.jpg')
end
