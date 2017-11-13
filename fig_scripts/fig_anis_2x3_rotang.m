%% fig_comp_achli_nemnas.m
%% - change 'height' of unit (rectangular) cell;
%% - compare to Delameter et al (1975);
%% - check effect of orthotropic host medium;

clear;
outdir   = 'out'
if ~exist(outdir,'dir')
   mkdir(outdir);
end

filename0   = 'out/anis_vs_rotang';
filename1   = 'out/anis_vs_phase';
filename    = 'out/anis_2x3_rotang';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A     = .5;
Bhat  = .5;
ahat  = A/2;
%%
nt          = 100;
theta_vec   = linspace(0,pi/2,nt)';
%%
alp_vec  = [.5 .75 1 1.25 1.5 1.75 2];
na       = length(alp_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist([filename0,'.mat'])
   disp('straight crack calc');
   %%
   Npolys   = 50;
   m1       = zeros(nt,na);
   m2       = m1;
   pdir     = m1;
   %%
   for j=1:nt
   %for j=nt
      for r=1:na
         th_hat   = theta_vec(j);
         alp      = alp_vec(r);
         Dalp     = diag([1 1/alp]);
         %%
         B  = Bhat*alp;
         th = atan(alp*tan(th_hat));
         a  = ahat*sqrt(cos(th_hat)^2+alp^2*sin(th_hat)^2);
         z  = a*exp(1i*th);
         AB = [A,B];
         %%
         Irr_vars    = CURVEget_strtline(-z,z);
         C_eff       = feval(crk_solver,AB,Irr_vars,Npolys);
         Mhat        = Dalp*C_eff{1}*Dalp;
         %%
         [U,Lam]  = eig(Mhat);
         m1(j,r)  = min(diag(Lam));
         m2(j,r)  = max(diag(Lam));
         %%
         jj          = find(diag(Lam)==m1(j,r));
         pdir(j,r)   = angle(U(1,jj)+1i*U(2,jj));
         %%
         if 0%alp==1
            [th_hat,th]
            [A,Bhat,ahat]
            [A,B,a]
            C_eff{1}
            Mhat
            [m1(j,r),m2(j,r),pdir(j,r)]
            pause
         end
         if j>1
            if pdir(j,r)-pdir(j-1,r)>.5*pi
               pdir(j,r)   = pdir(j,r)-pi;
            elseif pdir(j,r)-pdir(j-1,r)<-.5*pi
               pdir(j,r)   = pdir(j,r)+pi;
            end
         end
      end
      disp([num2str(j),' configurations done (out of ',num2str(nt),')']);
   end
   save([filename0,'.mat'],'alp_vec','theta_vec','A','Bhat','ahat','m1','m2','pdir');
else
   load([filename0,'.mat']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylab  = {'m1','m2','chi/pi'};
%yy    = {fliplr(m1),fliplr(m2),fliplr(pdir)/pi};
yy    = {fliplr(m1),fliplr(m2),fliplr(pdir)/pi-0.5};
j2    = 5;
lw    = 1.5;
%%
for j=1:3
   subplot(2,3,j);
   GEN_set_plotorder;
   plot(theta_vec/pi,yy{j}(:,1:4),'k');
   hold on;
   plot(theta_vec/pi,yy{j}(:,j2),'r','linewidth',lw);
   GEN_proc_fig('th/pi',ylab{j});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
crk_solver  = @BONE_GF_cracks_rect_cell;
crk_fxn     = @CURVEprof_cos;
%%
A     = .5;
Bhat  = .5;
a     = A/2;
boB   = .25;
%%
nt          = 100;
phase_vec   = linspace(0,1,nt)';
%%
alp_vec  = [.5 .75 1 1.25 1.5 1.75 2];
na       = length(alp_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist([filename1,'.mat'])
   disp('cos crack calc');
   %%
   Npolys   = 50;
   m1       = zeros(nt,na);
   m2       = m1;
   pdir     = m1;
   %%
   for j=1:nt
   %for j=nt
      for r=1:na
         phase = phase_vec(j);
         alp   = alp_vec(r);
         Dalp  = diag([1 1/alp]);
         %%
         B  = Bhat*alp;
         b  = boB*B;
%        th = atan(alp*tan(th_hat));
%        a  = ahat*sqrt(cos(th_hat)^2+alp^2*sin(th_hat)^2);
%        z  = a*exp(1i*th);
         AB = [A,B];
         %%
%        Irr_vars    = CURVEget_strtline(-z,z);
         crk_vars = {1,phase};
         srt      = {[a,b],0,[0 0]};
         Irr_vars = {crk_fxn,crk_vars,srt};
         C_eff    = feval(crk_solver,AB,Irr_vars,Npolys);
         Mhat     = Dalp*C_eff{1}*Dalp;
         %%
         [U,Lam]  = eig(Mhat);
         m1(j,r)  = min(diag(Lam));
         m2(j,r)  = max(diag(Lam));
         %%
         jj          = find(diag(Lam)==m1(j,r));
         pdir(j,r)   = angle(U(1,jj)+1i*U(2,jj));
         %%
         if 0%alp==1
            [th_hat,th]
            [A,Bhat,ahat]
            [A,B,a]
            C_eff{1}
            Mhat
            [m1(j,r),m2(j,r),pdir(j,r)]
            pause
         end
         if j>1
            if pdir(j,r)-pdir(j-1,r)>.5*pi
               pdir(j,r)   = pdir(j,r)-pi;
            elseif pdir(j,r)-pdir(j-1,r)<-.5*pi
               pdir(j,r)   = pdir(j,r)+pi;
            end
         end
      end
      disp([num2str(j),' configurations done (out of ',num2str(nt),')']);
   end
   save([filename1,'.mat'],'alp_vec','phase_vec','A','Bhat','a','boB','m1','m2','pdir');
else
   load([filename1,'.mat']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylab  = {'m1','m2','chi/pi'};
%yy    = {fliplr(m1),fliplr(m2),fliplr(pdir)/pi};
yy    = {fliplr(m1),fliplr(m2),fliplr(pdir)/pi-0.5};
j2    = 5;
%%
for j=1:3
   subplot(2,3,j+3);
   GEN_set_plotorder;
   plot(phase_vec,yy{j}(:,1:4),'k');
   hold on;
   plot(phase_vec,yy{j}(:,j2),'r','linewidth',lw);
   GEN_proc_fig('phi',ylab{j});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lets  = {'(a)','(b)','(c)','(d)','(e)','(f)'};
Y0 = [.2 .7 .475 .2 .85 .5];
Y1 = [1 1.02 .6 1 .95 .6];
if 1
   Y1([3 6])   = Y1([3 6])-.5;
   Y0([3 6])   = Y0([3 6])-.5;
end
X1 = [.5 .5 .5 1 1 1];
ax = .1;
ay = [.87 .12 .12 .87 .87 .87];

for j=1:6
   subplot(2,3,j), hold on;
   x  = ax*X1(j);
   y  = Y0(j)+ay(j)*(Y1(j)-Y0(j));
   ylim([Y0(j),Y1(j)]);

   if mod(j,3)==0
      set(gca,'YTick',0:.03:.09);
   end

   txt   = text(x,y,lets{j});
   set(txt,'FontName','Times','FontSize',16);
   hold off;
end


saveas(gcf,[filename,'.eps'],'epsc');
%!gv out/anis_2x3_rotang.eps &
