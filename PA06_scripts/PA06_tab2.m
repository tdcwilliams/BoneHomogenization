%% table2_WA06.m
%% - effect of changing contrast of fibres:

outdir   = 'outPA06'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
basedir  = 'outPA06/MPfib_datfiles/';
if ~exist(baseir,'dir')
   mkdir(baseir);
end

mvec     = [6 20 120 400];% 49 51];
rvec     = [.356825,.418414,.472035,.4886025];
Nterms   = 15;
%%
AB          = .5*[1 1];
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of circle

for j=1:length(mvec)%%use some test inputs:
  midi   = [mvec(j) 1];
  %%
  for r=1:length(rvec)
    radius     = rvec(r);
    srt        = {radius*[1 1],0,[0 0]};%area=pi*radius^2
    Irr_vars   = {crk_fxn,crk_prams,srt,midi};
    C_eff      = BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
    MSTAR(j,r) = C_eff{1}(1);
  end
end

disp('Table 2:');
disp(MSTAR);
%%
save([basedir,'tab2.dat'],'MSTAR','-ascii','-tabs');
