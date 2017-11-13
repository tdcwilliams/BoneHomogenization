%% table1_WA06_cavs.m
%% - effect of changing radii of cavities:
clear;

outdir   = 'outPA06'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
basedir  = 'outPA06/MPcav_datfiles/';
if ~exist(baseir,'dir')
   mkdir(baseir);
end

Nvec  = [1 3 5 11 19];% 49 51];
rvec  = [.05 .1 .15 .2 .25 .3 .35 .4 .45];% .495];
%%
AB          = .5*[1 1];
FXN         = @BONE_MP_cavs_rect_cell;
FXN2        = @BONE_MPG_cavs_rect_cell;
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of circle

for j=1:length(Nvec)%%use some test inputs:
  Nterms = Nvec(j);
  N_MP   = Nterms;
  %%
  for r=1:length(rvec)
    radius     = rvec(r);
    srt        = {radius*[1 1],0,[0 0]};%area=pi*radius^2
    Irr_vars   = {crk_fxn,crk_prams,srt};
    C_eff      = feval(FXN,AB,Irr_vars,Nterms);
    MSTAR(j,r) = C_eff{1}(1);
    %%
    srt_out       = {1.05*radius*[1 1],0,[0 0]};
    Irr_vars_out  = {crk_fxn,crk_prams,srt_out};
    C_eff2        = feval(FXN2,AB,Irr_vars_out,{Irr_vars},N_MP,Nterms);
    MSTAR2(j,r)   = C_eff2{1}(1);
  end
end

disp('Table 1 (cavs):');
disp(MSTAR);
disp(MSTAR2);
%%
save([basedir,'tab1-cavs-mp.dat'],'MSTAR','-ascii','-tabs');
