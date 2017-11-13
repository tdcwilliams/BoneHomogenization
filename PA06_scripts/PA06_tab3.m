%% table3_WA06.m
%% - effect of having 2 different fibres:
clear;

outdir   = 'outPA06'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
basedir  = 'outPA06/MPfib_datfiles/';
if ~exist(baseir,'dir')
   mkdir(baseir);
end

Nvec  = [1 3 5 7 9 13];
rvec  = [.1 .2 .3 .4 .45 .495];
%%
AB          = .5*[2 1];
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of circle

if 1%% value in PA06:
  midi1  = [1e2 1];
  midi2  =[1e-2 1];
else%% check against orthotropic case:
  disp('CHECKING AGAINST ORTHOTROPIC CASE:');
  midi1  = [10 1];
  midi2  = [10 1];
end

for j=1:length(Nvec)%%use some test inputs:
  Nterms = Nvec(j)
  %%
  for r=1:length(rvec)
    radius        = rvec(r);
    srt1          = {radius*[1 1],0,[-0.5 0]};
    srt2          = {radius*[1 1],0,[0.5 0]};
    Irr_vars(1,:) = {crk_fxn,crk_prams,srt1,midi1};
    Irr_vars(2,:) = {crk_fxn,crk_prams,srt2,midi2};
    C_eff         = BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
    MxSTAR(j,r)   = C_eff{1}(1,1);
    MySTAR(j,r)   = C_eff{1}(2,2);
  end
end

disp('Table 3, m_x:');
disp(MxSTAR);
disp('Table 3, m_y:');
disp(MySTAR);
%%
save([basedir,'tab3a.dat'],'MxSTAR','-ascii','-tabs');
save([basedir,'tab3b.dat'],'MySTAR','-ascii','-tabs');
