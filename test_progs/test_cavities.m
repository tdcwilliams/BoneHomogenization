%% test_cavities.m

AB=.5*[1 1];
Nterms=10;

crk_fxn_out=@CURVEprof_circarc;
crk_prams_out={1};%% fraction of circle
radius_out=.4
rot_out=0;
srt_out={radius_out*[1 1],rot_out,[0 0]};
Irr_vars_out={crk_fxn_out,crk_prams_out,srt_out};

if 1
  crk_fxn=@CURVEprof_crooked_egg;
  crk_prams={[1 1 1]};
  scaler=.25*[1 1];
  rot=70;
  srt={scaler,rot,[0 0]};%area=pi*radius^2
  midi=[10 1]
else
  crk_fxn=@CURVEprof_circarc;
  crk_prams={1};%% fraction of circle
  radius=.225
  rot=30
  srt={radius*[1.25 1],rot,[0 0]};%area=pi*radius^2
  midi=[10 1]
end

Irr_vars={crk_fxn,crk_prams,srt,midi};
C_eff=BONE_MP_cavs_rect_cell(AB,Irr_vars,Nterms);
Mmat_MP=C_eff{1}
[mm,angs]=BONE_proc_M(Mmat_MP);
{mm(1),mm(2),angs(1)*180/pi},pause
%%
BC=1;
fib_stuff={midi(1),midi(2),BC};
srt{2}=rot+90;
Irr_vars2={crk_fxn,crk_prams,srt,fib_stuff};
C_eff=BONE_GF_fibres_rect_cell(AB,Irr_vars2,Nterms);
Mmat_GF=C_eff{1}
%%
C_eff=BONE_GF_fibres_rect_cell_Tn(AB,Irr_vars2,2*Nterms);
Mmat_GF2=C_eff{1}
%%
N_MP=Nterms;
C_eff=BONE_MPG_fibres_rect_cell(AB,Irr_vars_out,...
  {Irr_vars2},N_MP,Nterms);
Mmat_MPG=C_eff{1}
%%
C_eff=BONE_MPG_fibres_rect_cell_Tn(AB,Irr_vars_out,...
  {Irr_vars2},N_MP,Nterms);
Mmat_MPG2=C_eff{1}

[mm,angs]=BONE_proc_M(Mmat_GF);
{mm(1),mm(2),angs(1)*180/pi}