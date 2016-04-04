function [rvecs,ds_dt,th_vec, dtheta_ds,d2s_dt2,d2theta_ds2,...
   d2r_ds2,LC,tvec,Area,xtra_out]=...
     BONE_get_rsdtheta_NRquick(crk_vars,svec_on_LC);

fxn_handle  = crk_vars{1};%%which routine to use for crack shape
crk_prams   = crk_vars{2};%%parameters needed for given shape
srt         = crk_vars{3};%%scaling,rotation (in degrees),translations
Nint        = length(svec_on_LC);

%% calc length of curve:
[LC,dsLC,TW]   = ...
  feval(@arc_length,1,fxn_handle,...
          crk_prams,srt);
LC = LC/2;%% halve length so s\in[-LC,LC];


if nargout==10%% calc area bounded by curve
  Area   = get_area(fxn_handle,crk_prams,srt);
end

%% now find which values of t give desired value of s,
%% and thus find the required quadrature points.
tvec     = 0*svec_on_LC;
interval = [-1 1];

%% DON'T NEED TO SEARCH FOR t IF s/L=1,-1:
jp1         = find(svec_on_LC==1);
tvec(jp1)   = 1;
jm1         = find(svec_on_LC==-1);
tvec(jm1)   = -1;
jnot        = find( (svec_on_LC~=1)&(svec_on_LC~=-1) );
Npts        = length(jnot);
%%
if Npts>0
  s_target=( 1+svec_on_LC(jnot) )*LC;
  guess=svec_on_LC(jnot);

  if 0%% test arc_length_quick:
    LC_=arc_length_quick(1,fxn_handle,crk_prams,srt,TW);
   [LC,LC_/2]
  elseif 0
    s_=arc_length_quick(guess,fxn_handle,crk_prams,srt,TW);
    for j=1:length(guess)
      s_(j,2)=arc_length(guess(j),fxn_handle,crk_prams,srt);
    end
%    s_
  end

  [tvec(jnot),flag]=GEN_findroot_NRsafe2(...
    @arc_length_quick,interval,guess,s_target,fxn_handle,crk_prams,srt,TW);
end

%for j=1:length(tvec)
%  tst(j,:)=[svec_on_LC(j)*LC,...
%    feval(@arc_length,tvec(j),fxn_handle,crk_prams,srt,0)-LC];
%end
%tst,LC,return

%% calc r,dr,d2r for basic function
if 0
  [dr,rvecs,d2r,d3r,phi_vec]=feval(fxn_handle,tvec,crk_prams);
else
  [dr,rvecs,d2r,d3r]=feval(fxn_handle,tvec,crk_prams);
end

%% THEN ALLOW FOR SCALING, ROTATIONS AND TRANSLATIONS
scaler   = srt{1};
M_sc     = diag(scaler);
th_rot   = pi*srt{2}/180;%% angle in radians:
M_rot    = [cos(th_rot),-sin(th_rot);sin(th_rot),cos(th_rot)];
v_trans  = srt{3};
rvecs    = diag(v_trans)*ones(size(rvecs))+M_rot*M_sc*rvecs;
%%
dr    = M_rot*M_sc*dr;
dx    = dr(1,:)';
dy    = dr(2,:)';
ds_dt = sqrt(dx.^2+dy.^2);

%% CALC \theta & MAKE IT CONTINUOUS:
if 0
  th_vec = th_rot+...
    angle( scaler(1)*cos(phi_vec)+i*scaler(2)*sin(phi_vec) );
else
  th_vec = angle(dx+1i*dy);
end

for j=2:length(th_vec)
  if th_vec(j)-th_vec(j-1)>1.75*pi
    th_vec(j:end) = th_vec(j:end)-2*pi;
  end
  if th_vec(j)-th_vec(j-1)<-1.75*pi
    th_vec(j:end) = th_vec(j:end)+2*pi;
  end
end


%[ds_dt,th_vec]
%tangvecs=[ dx./ds_dt dy./ds_dt ]';
%%
d2r         = M_rot*M_sc*d2r;
d2x         = d2r(1,:)';
d2y         = d2r(2,:)';
dtheta_ds   = ( d2y.*dx - d2x.*dy )./ds_dt.^3;
% - LB (27.08.09) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! - %
% dtheta_ds=( d2y.*dx - d2x.*dy )./ds_dt.^2; dtheta_ds=dtheta_ds/LC;
% ----------------- %
d2s_dt2  = (d2x.*dx+d2y.*dy)./ds_dt;
%%
d3r      = M_sc*d3r;
d3x      = d3r(1,:).';
d3y      = d3r(2,:).';
d2theta_dt2 = (d3y.*dx-d3x.*dy)./ds_dt.^2-2*d2s_dt2.*dtheta_ds;
d2theta_ds2 = (d2theta_dt2-d2s_dt2.*dtheta_ds)./ds_dt.^2;
%%
d2r_ds2  = [(d2x.*ds_dt-d2s_dt2.*dx)./ds_dt.^3,...
                (d2y.*ds_dt-d2s_dt2.*dy)./ds_dt.^3];

%% integrate wrt tau=s/LC
%%  => ds/d(tau)=LC & d^2s/dtau^2=d(LC)/d(tau)=0:
%  dS_dt=0*svec_on_LC;
%Svec=LC+dS_dt;%%%%???????????

if nargout==11
  xx  = rvecs(1,:).';
  yy  = rvecs(2,:).';
  [XX0,XX]  = meshgrid(xx,xx);
  [YY0,YY]  = meshgrid(yy,yy);
  [SS0,SS]  = meshgrid(svec_on_LC,svec_on_LC);
  %%
  ZZ  = (XX-XX0)+1i*(YY-YY0);
  Theta  = angle( ZZ./(SS-SS0) );
  %%
  for j=1:Nint
    Theta(j,j) = th_vec(j);
    for r=j-1:-1:1
      if Theta(j,r)-Theta(j,r+1)>1.75*pi
        Theta(j,1:r) = Theta(j,1:r)-2*pi;
      end
      if Theta(j,r)-Theta(j,r+1)<-1.75*pi
        Theta(j,1:r) = Theta(j,1:r)+2*pi;
      end
    end
    for r=j+1:Nint
      if Theta(j,r)-Theta(j,r-1)>1.75*pi
        Theta(j,r:end)  = Theta(j,r:end)-2*pi;
      end
      if Theta(j,r)-Theta(j,r-1)<-1.75*pi
        Theta(j,r:end)  = Theta(j,r:end)+2*pi;
      end
    end
%  Theta(j,1:20)/pi
%  plot(svec_on_LC,Theta(j,:).'/pi), ylim([.49 2.51]), pause;
  end
  xtra_out  = {TW,Theta};
end


function [s,ds,TW]=...
           arc_length(t,fxn_handle,crk_prams,srt)
%% f=s-s_target => df=ds
%% t a scalar

scaler   = srt{1};%% don't need to worry about rotations
              %% & translations when calc'ing arc length:
M_sc  = diag(scaler);
dr    = M_sc*feval(fxn_handle,t,crk_prams);
ds    = sqrt( dr(1,:).^2+dr(2,:).^2 );
%%
Ngl   = 0;
tol   = 1e-14;
err   = 2*tol;
s0    = 0;

while err>tol
  Ngl       = Ngl+50;
  [tgl,wgl] = OP_numint_legendre(Ngl,[-1 t]);
  dr        = M_sc*feval(fxn_handle,tgl,crk_prams);
  ds_       = sqrt( dr(1,:).^2+dr(2,:).^2 );
  s         = ds_*wgl;
  %%
  err = abs(1-s0/s);
  s0  = s;
end

%% get quadrature points & weights relative to [-1,1]:
TW = {(1+2*tgl-t)/(1+t),2/(1+t)*wgl};

function [s,ds]=...
  arc_length_quick(tvec,fxn_handle,crk_prams,srt,TW)
%% s_target can be a vector of same size as tvec;

if length(TW)==1
  Npts         = TW{1};
  [tgl0,wgl0]  = OP_numint_legendre(Npts);
else
  tgl0   = TW{1};
  wgl0   = TW{2};
end

scaler   = srt{1};%% don't need to worry about rotations
              %% & translations when calc'ing arc length:
M_sc  = diag(scaler);
dr    = M_sc*feval(fxn_handle,tvec,crk_prams);
ds    = sqrt( dr(1,:).^2+dr(2,:).^2 ).';
s     = 0*ds;
%%
for j=1:length(tvec)
  t      = tvec(j);
  wgl    = (1+t)/2*wgl0;
  tau    = ( (1+t)*tgl0+t-1 )/2;
  dr     = M_sc*feval(fxn_handle,tau,crk_prams);
  s(j)   = sqrt( dr(1,:).^2+dr(2,:).^2 )*wgl;
end

function Area=get_area(fxn_handle,crk_prams,srt)
%% f=s-s_target => df=ds
%% t a scalar

scaler   = srt{1};%% don't need to worry about rotations
                 %% & translations when calc'ing area:
%  M_sc=diag(scaler);
%  [dr,rvecs]=M_sc*feval(fxn_handle,t,crk_prams);
%  ds=sqrt( dr(1,:).^2+dr(2,:).^2 );
%%
Ngl   = 0;
tol   = 1e-8;
err   = 2*tol;
Area0 = 0;
while err>tol
  Ngl          = Ngl+50;
  [tgl,wgl]    = OP_numint_legendre(Ngl,[-1 1]);
  [dr,rvecs]   = feval(fxn_handle,tgl,crk_prams);
  %%
  ig     = prod(scaler)/2*(  rvecs(1,:).*dr(2,:)+...
              - rvecs(2,:).*dr(1,:) );
  Area   = ig*wgl;
  %%=.5*\int_{-1}^1[x(t).y'(t)-y(t).x'(t)]dt

  err    = abs(1-Area0/Area);
  Area0  = Area;
end

return;

Np    = 1e4;
tol   = 1e-8;
err   = 2*tol;
Area0 = 0;
while err>tol
  Np        = Np+50;
  tt        = -1+(0:Np-1)'*2/Np;
  [dr,xy0]  = feval(fxn_handle,tt,crk_prams);
  xy        = xy0'*M_sc;
  %%
  Area   = GEN_area(xy);
%  plot(xy(:,1),xy(:,2)),pause
  err    = abs(1-Area0/Area);%{Np,Area,pi/100,err},pause
  Area0  = Area;
end
