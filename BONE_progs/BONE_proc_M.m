function [mm,angs]=BONE_proc_M(Mmat)

a=Mmat(1);
d=Mmat(4);
b=Mmat(2);
%%
sig=(a+d)/2;
del=(a-d)/2;
%%
m_max=sig+sqrt(b^2+del^2);
m_min=sig-sqrt(b^2+del^2);
mm=[m_min,m_max]';
%%
if b~=0
  ang_min=atan((a-m_min)/b);
elseif a>=d
  ang_min=pi/2;
else
  ang_min=0;
end
ang_max=ang_min+pi/2;
angs=[ang_min,ang_max]';