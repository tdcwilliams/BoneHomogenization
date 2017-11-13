%% table1_WA06.m
%% - effect of changing radii of fibres:
clear;

outdir   = 'outPA06'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
basedir  = 'outPA06/MPfib_datfiles/';
if ~exist(baseir,'dir')
   mkdir(baseir);
end

Nvec  = [1 3 5 11 19];% 49 51];
rvec  = [.05 .1 .15 .2 .25 .3 .35 .4 .45 .495];
%%
AB          = .5*[1 1];
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of circle

for j=1:length(Nvec)%%use some test inputs:
  Nterms = Nvec(j)
  %%
  for r=1:length(rvec)
    radius     = rvec(r);
    srt        = {radius*[1 1],0,[0 0]};%area=pi*radius^2
    midi       = [10 1];
    Irr_vars   = {crk_fxn,crk_prams,srt,midi};
    C_eff      = BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
    MSTAR(j,r) = C_eff{1}(1);
  end
end

disp('Table 1:');
disp(MSTAR);
%%
save([basedir,'tab1.dat'],'MSTAR','-ascii','-tabs');
