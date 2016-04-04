%%WP06_fig10.m

np=50;
r1vec=(1:np)'/np*.5;
mx=0*r1vec;
%%
AB=.5*[1 1];
Nterms=7;

for j=1:np%nargin==0%%use some test inputs:
  if 1
    crk_fxn=@CURVEprof_circarc;
    crk_prams={1};%% fraction of circle
    radius=r1vec(j);
    srt={radius*[1 1],0,[0 0]};%area=pi*radius^2
  end
  Irr_vars={crk_fxn,crk_prams,srt};
  C_eff=BONE_MP_cavs_rect_cell(AB,Irr_vars,Nterms);
  mx(j)=C_eff(1);
%    Irr_vars={irr_vars};
%  Nterms=5;
end

plot(r1vec,mx,'r')