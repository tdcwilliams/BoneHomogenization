%% test_BONE_MPG_cavs_rect_cell.m
clear;
Nterms=30
N_MP=Nterms;
DO_ELLIPSE_out=1;
DO_ELLIPSE_in=1;
%%
AB=.5*[1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shape of buffer zone:
if DO_ELLIPSE_out
  crk_fxn=@CURVEprof_circarc;
%    crk_prams_out={[1,+1]};%% fraction of circle, anticlockwise;
%    crk_prams_out={1};
  crk_prams_out={[1,-1]};
  radius_out=.4;
  scaler_out=radius_out*[1 1]
  disp(sprintf('outer ellipse eccentricity %d',...
    scaler_out(1)/scaler_out(2)));
  srt_out={scaler_out,0,[0 0]};
  Irr_vars_outer={crk_fxn,crk_prams_out,srt_out}
else
  crk_fxn=@CURVEprof_crooked_egg;
  crk_prams_out={[.25, 1, -1]};
  srt_out={[1 .8],0,[.075 -.1]};
  Irr_vars_outer={crk_fxn,crk_prams_out,srt_out}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shape of cavity:
if DO_ELLIPSE_in
  crk_fxn_in=@CURVEprof_circarc
%    crk_prams_in={[1,-1]};%% fraction of circle, clockwise;
  crk_prams_in={1};
  radius_in=.225;
  scaler_in=radius_in*[1.5 1.5]
  srt_in={scaler_in,0,[0 0]};%area=pi*radius^2
  Irr_vars_inner{1}={crk_fxn_in,crk_prams_in,srt_in};
else
  crk_fxn_in=@CURVEprof_crooked_egg;
  crk_prams_in={[.25, 1, 1]};
  srt_in={[1 .8]/2,0,[.075 -.1]/2};
  Irr_vars_inner{1}={crk_fxn_in,crk_prams_in,srt_in};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_eff_hybrid=BONE_MPG_cavs_rect_cell(AB,Irr_vars_outer,...
        Irr_vars_inner,N_MP,Nterms)
C_eff_MP=BONE_MP_cavs_rect_cell(AB,Irr_vars_inner{1},N_MP)
%  C_eff_G=BONE_cavities_rect_cell(AB,Irr_vars_inner{1},N_MP)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_ELLIPSE_in==1
  crk_prams_in={[1,-1]};
  Irr_vars_inner{1}={crk_fxn_in,crk_prams_in,srt_in};
else
  crk_prams_in={[.25 1 -1]};
  Irr_vars_inner{1}={crk_fxn_in,crk_prams_in,srt_in};
end
if DO_ELLIPSE
  crk_prams_out={1};
else
  crk_prams_out={[.25 1 1]};
end
Irr_vars_outer={crk_fxn,crk_prams_out,srt_out};
C_eff=BONE_MPG_cavs_rect_cell(AB,Irr_vars_outer,...
        Irr_vars_inner,N_MP,Nterms)