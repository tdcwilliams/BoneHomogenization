%% PA06_tab5or6.m
%% - effect of changing eccentricity of fibres:
clear;
warning off;
disp('square unit cell, single elliptical fibre (b is constant):')
disp('~|a');
disp('___');
disp('N|~');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir   = 'outPA06'
if ~exist(outdir,'dir')
   mkdir(outdir);
end

METHOD   = 3
if METHOD==1%%MP
   basedir  = 'outPA06/MPfib_datfiles/';
elseif METHOD==2%%GF
   basedir  = 'outPA06/GFfib_datfiles/';
elseif METHOD==3%%MPG
   basedir  = 'outPA06/MPGfib_datfiles/';
end
if ~exist(baseir,'dir')
   mkdir(baseir);
end

NO = 6;%% TABLE NO:
%%
Nvec  = [1 3 5 11 15 19];% 49 51];
avec  = [.1 .2 .3 .4];
%%
AB          = .5*[1 1];
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of ellipse
midi        = [10 1];

if NO==5
  b   = .1;
else
  b   = .25;
end

for j=1:length(Nvec)
  Nterms = Nvec(j);
  %%
  for r=1:length(avec)
    a = avec(r);
%      srt={[a b],0,[0 0]};
    srt           = {[a b],0,[0 0]};
    if METHOD==1
       Irr_vars   = {crk_fxn,crk_prams,srt,midi};
       C_eff      = BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
    elseif METHOD==2
       BC         = 0;
       Irr_vars   = {crk_fxn,crk_prams,srt,{midi(1),midi(2),BC}};
       C_eff      = BONE_GF_fibres_rect_cell(AB,Irr_vars,Nterms);
    elseif METHOD==3
       BC                  = 0;
       Irr_vars_inner{1}   = {crk_fxn,crk_prams,srt,{midi(1),midi(2),BC}};
       srt_out             = srt;
       %%
       r0               = max(a,b);
       r1               = min(AB);
       rad_out          = r0+.6*(r1-r0);
       srt_out{1}       = rad_out*[1 1];
       Irr_vars_outer   = {crk_fxn,crk_prams,srt_out};
       %%
       N_MP    = 10;
       C_eff   = BONE_MPG_fibres_rect_cell(AB,Irr_vars_outer,...
                     Irr_vars_inner,N_MP,Nterms);
    end
    MxSTAR(j,r)   = C_eff{1}(1);
    MySTAR(j,r)   = C_eff{1}(4);
%      disp(sprintf('m_xy  = %d',C_eff(2)));
  end
end

if NO==5
   disp('Table 5, m_x:');
   disp(MxSTAR);
   disp('Table 5, m_y:');
   disp(MySTAR);
   %%
   save([basedir,'tab5a.dat'],'MxSTAR','-ascii','-tabs');
   save([basedir,'tab5b.dat'],'MySTAR','-ascii','-tabs');
else
   disp('Table 6, m_x:');
   disp(MxSTAR);
   disp('Table 6, m_y:');
   disp(MySTAR);
   %%
   save([basedir,'tab6a.dat'],'MxSTAR','-ascii','-tabs');
   save([basedir,'tab6b.dat'],'MySTAR','-ascii','-tabs');
end
