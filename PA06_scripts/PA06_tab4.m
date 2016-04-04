%% table4_WA06.m
%% - effect of moving position of 1 out of 4 fibres:
clear;

METHOD   = 3
if METHOD==1%%MP
   basedir  = 'outPA06/MPfib_datfiles/';
elseif METHOD==2%%GF
   basedir  = 'outPA06/GFfib_datfiles/';
elseif METHOD==3%%MPG
   basedir  = 'outPA06/MPGfib_datfiles/';
end

Nvec  = [1 3 5 7 11 13 15];
rvec  = [.1 .2 .3 .35 .4];
%%
AB          = .5*[2 2];
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of circle
%%
midi  = [10 1];
if 1%% value in WA06:
  pert   = [.1 -.2];
else%% check against orthotropic case [DONE -> OK]:
  disp('CHECKING AGAINST ORTHOTROPIC CASE:');
  pert   = [0 0];
end
centres  = [.5*[1 -1];.5*[1 1];.5*[-1 1];...
         .5*[-1 -1]+pert ];

for j=1:length(Nvec)%%use some test inputs:
  Nterms = Nvec(j);
  %%
  for r=1:length(rvec)
    radius  = rvec(r);
    if METHOD==1
       for q=1:4
         srt            = {radius*[1 1],0,centres(q,:)};
         Irr_vars(q,:)  = {crk_fxn,crk_prams,srt,midi};
       end
       C_eff      = BONE_MP_fibres_rect_cell(AB,Irr_vars,Nterms);
    elseif METHOD==2
       BC         = 0;
       for q=1:4
         srt            = {radius*[1 1],0,centres(q,:)};
         Irr_vars(q,:)  = {crk_fxn,crk_prams,srt,{midi(1),midi(2),BC}};
       end
       C_eff      = BONE_GF_fibres_rect_cell(AB,Irr_vars,Nterms);
    elseif METHOD==3
       BC                  = 0;
       for q=1:4
          srt               = {radius*[1 1],0,centres(q,:)};
          Irr_vars_inner{q} = {crk_fxn,crk_prams,srt,{midi(1),midi(2),BC}};
          %%
          srt_out             = srt;
          %%
          r0                  = radius;
          r1                  = min(AB);
          rad_out             = r0+.4*(r1-r0);
          srt_out{1}          = rad_out*[1 1];
          Irr_vars_outer(q,:) = {crk_fxn,crk_prams,srt_out};
       end
       %%
       N_MP    = 15;
       C_eff   = BONE_MPG_fibres_rect_cell(AB,Irr_vars_outer,...
                     Irr_vars_inner,N_MP,Nterms);
    end
    MxSTAR(j,r)   = C_eff{1}(1);
    MxySTAR(j,r)  = C_eff{1}(2);
    MySTAR(j,r)   = C_eff{1}(4);
  end
end

disp('Table 4, m_x:');
disp(MxSTAR);
disp('Table 4, m_y:');
disp(MySTAR);
disp('Table 4, m_xy:');
disp(MxySTAR);
%%
save([basedir,'tab4a.dat'],'MxSTAR','-ascii','-tabs');
save([basedir,'tab4b.dat'],'MySTAR','-ascii','-tabs');
save([basedir,'tab4c.dat'],'MxySTAR','-ascii','-tabs');
