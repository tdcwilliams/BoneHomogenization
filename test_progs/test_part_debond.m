%% test_part_debond.m
clear;
Nterms0  = 5;
Nterms   = 7;
NtermsTn = 20;
scaler   = [1.2 1];
m_fib    = 5;
%%
crk_fxn_out    = @CURVEprof_circarc;
crk_prams_out  = {1};%% fraction of circle
radius_out     = .4

if 0
  crk_fxn   = @CURVEprof_crooked_egg;
  crk_prams = {[1 1 1]};
  scaler    = .25*[1 1];
  rot       = 70;
  srt       = {scaler,rot,[0 0]};%area=pi*radius^2
  midi      = [10 1]
else
  %% 1x2, 1x3 or 1x4 unit cell;
  Nfib      = 2;
  crk_fxn   = @CURVEprof_circarc;
  AB        = .5*[1 1];
  alpc      = .5;%% fraction of circle that is debonded
  crk_prams = {{alpc},{1-alpc}};%% fraction of circle
  radius    = {.15,.15};
  centre    = {[0,0],[0 0]};
  rot       = {0,180};
  bc        = {1,0};%%bc=1 debonded, bc=0 cts
  midi_gf   = {{m_fib 1 bc{1}},{m_fib 1 bc{2}}};
  for r=1:Nfib
    srt              = {radius{r}*scaler,rot{r},centre{r}};%area=pi*radius^2
    Irr_vars_gf(r,:) = {crk_fxn,crk_prams{r},srt,midi_gf{r}};
    %%
    Irr_vars_mpg_inner{r}  = {crk_fxn,crk_prams{r},srt,midi_gf{r}};
    %%
    srt  = {radius_out*[1 1],rot{r},centre{r}};%area=pi*radius^2
    Irr_vars_mpg_outer(r,:)   = {crk_fxn_out,crk_prams_out,srt};
  end
  Irr_vars_mpg_inner = {Irr_vars_gf};
  srt                = {radius_out*[1 1],0,centre{1}};%area=pi*radius^2
  Irr_vars_mpg_outer = {crk_fxn_out,crk_prams_out,srt};
  

  %% if m_fib=1 then can compare with crack results
  crk_fxn   = @CURVEprof_circarc;
  crk_prams = {alpc};%% fraction of circle
  radius    = {.15};
  centre    = {[0,0]};
  rot       = {0};
  bc        = 0;
  Nfib0     = 1
  for r=1:Nfib0
    srt                 = {radius{r}*scaler,rot{r},centre{r}};%area=pi*radius^2
    Irr_vars_gf0(r,:)   = {crk_fxn,crk_prams,srt};
    %%
    Irr_vars_mpg_inner0{r}    = {crk_fxn,crk_prams,srt};
    srt                       = {radius_out*[1 1],rot{r},centre{r}};%area=pi*radius^2
    Irr_vars_mpg_outer0(r,:)  = {crk_fxn_out,crk_prams_out,srt};
  end
end

if 1
  C_eff        = BONE_GF_fibres_rect_cell_Tn(AB,Irr_vars_gf,NtermsTn);
  Mmat_GF_Tn   = C_eff{1} 

% C_eff        = BONE_MPG_fibres_rect_cell_Tn(...
%   AB,Irr_vars_mpg_outer,Irr_vars_mpg_inner,7,NtermsTn);
% Mmat_MPG_Tn  = C_eff{1} 

  [mm,angs] = BONE_proc_M(Mmat_GF_Tn);
  {mm(1),mm(2),angs(1)*180/pi},return;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  C_eff        = BONE_GF_cracks_rect_cell(AB,Irr_vars_gf0,NtermsTn);
  Mmat_GF_crk  = C_eff{1}

  C_eff        = BONE_MPG_cracks_rect_cell(...
    AB,Irr_vars_mpg_outer0,Irr_vars_mpg_inner0,7,NtermsTn);
  Mmat_MPG_crk = C_eff{1} 

  [mm,angs] = BONE_proc_M(Mmat_MPG_crk);
  {mm(1),mm(2),angs(1)*180/pi}
end
return;
%%
BC          = [0;0];
Irr_vars2   = Irr_vars;
for r=1:Nfib
  midi            = Irr_vars{r,4};
  fib_stuff       = {midi(1),midi(2),BC(r)};
  Irr_vars2{r,4}  = fib_stuff;
end

if 0
  C_eff     = BONE_GF_fibres_rect_cell(AB,Irr_vars2,Nterms);
  Mmat_GF   = C_eff{1}
end
%%
if 0
  C_eff     = BONE_GF_fibres_rect_cell_Tn(AB,Irr_vars2,2*Nterms);
  Mmat_GF2  = C_eff{1}
end
%%
crk_fxn     = @CURVEprof_circarc;
crk_prams   = {1};%% fraction of circle
radius      = {.225,.225};
centre      = {[-.25,0],[.25,0]};
rot         = {0,0};
for r=1:Nfib
  srt                   = {radius{r}*[1.25 1],rot{r},centre{r}};%area=pi*radius^2
  Irr_vars_outer(r,:)   = {crk_fxn,crk_prams,srt};
  Irr_vars_inner{r}     = Irr_vars2(r,:);
end

return;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_MP     = Nterms;
C_eff    = BONE_MPG_fibres_rect_cell(AB,Irr_vars_outer,...
  Irr_vars_inner,N_MP,Nterms);
Mmat_MPG = C_eff{1}
%%
C_eff       = BONE_MPG_fibres_rect_cell_Tn(AB,Irr_vars_out,...
  {Irr_vars2},N_MP,Nterms);
Mmat_MPG2   = C_eff{1}

%[mm,angs]  = BONE_proc_M(Mmat_GF);
%{mm(1),mm(2),angs(1)*180/pi}
