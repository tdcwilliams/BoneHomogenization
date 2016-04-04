%% fig_comp_achli_nemnas.m
%% - change 'height' of unit (rectangular) cell;
%% - compare to Delameter et al (1975);
%% - check effect of orthotropic host medium;

clear;
filename = 'out/comp-achli-nemnas';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A  = .5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np       = 25;
%aoB_vec  = (1:np)'/np*10;
BoA_vec  = (1:np)'/np*2;
%aoA_vec  = .2:.2:.8;
aoA_vec  = .2:.3:.8;
%%
x0       = .05;
aoA_vec2 = x0+(1:np)'/np*(1-2*x0);
%aoB_vec2 = [.5, 1 2 4];
%BoA_vec2 = [.1 .25 .5 1];
BoA_vec2 = [.1 .5 1];
nc       = length(aoA_vec);

Ach_Li         = zeros(np,nc,2);
DD             = zeros(np,nc,2);
selfcon        = zeros(np,nc,2);
differential   = zeros(np,nc,2);
NemNas         = zeros(np,nc,2);

if ~exist([filename,'.mat'])
   Npolys   = 50;
   %%
   for j=1:np
      for r=1:nc
         a  = aoA_vec(r)*A;
         B  = BoA_vec(j)*A;
         AB = [A,B];
         %%
         Irr_vars    = CURVEget_strtline(-a,a);
         C_eff       = feval(crk_solver,AB,Irr_vars,Npolys);
         m_y(j,r,1)  = C_eff{1}(4);%[m_y(j),Ach_Li(j)]
         %%
         %Ach_Li(j,r,1)   = 1./(1-2*A/pi/B*log(cos(pi/2*a/A)));
         Ach_Li(j,r,1)        = BONEapp_AchLi(a,A,B);
         DD(j,r,1)            = BONEapp_DD(a,A,B);
         selfcon(j,r,1)       = BONEapp_selfcon(a,A,B);
         differential(j,r,1)  = BONEapp_differential(a,A,B);
         NemNas(j,r,1)        = BONEapp_NemNas(a,A,B);
         %%
         a  = aoA_vec2(j)*A;
         B  = A*BoA_vec2(r);
         %B  = a/aoB_vec2(1);
         AB = [A,B];
         %%
         Irr_vars    = CURVEget_strtline(-a,a);
         C_eff       = feval(crk_solver,AB,Irr_vars,Npolys);
         m_y(j,nc+1-r,2)  = C_eff{1}(4);%[m_y(j),Ach_Li(j)]
         %%
         %Ach_Li(j,r,2)   = 1./(1-2*A/pi/B*log(cos(pi/2*a/A)));
         Ach_Li(j,nc+1-r,2)         = BONEapp_AchLi(a,A,B);
         selfcon(j,nc+1-r,2)        = BONEapp_selfcon(a,A,B);
         differential(j,nc+1-r,2)   = BONEapp_differential(a,A,B);
         NemNas(j,nc+1-r,2)         = BONEapp_NemNas(a,A,B);
      end
      disp([num2str(j),' configurations done (out of ',num2str(np),')']);
   end
   save([filename,'.mat'],'BoA_vec','aoA_vec','BoA_vec2','aoA_vec2','A',...
         'm_y','Ach_Li','NemNas','selfcon','differential','DD');
else
   load([filename,'.mat']);
end

xlab  = {'B/A','a/A'};
xx    = {BoA_vec,aoA_vec2};
%XX    = {[2*x0 1-2*x0],aoA_vec2([1 end])};
XX    = {[2*x0 1],aoA_vec2([1 end])};
ax    = .70*[1 1 1];
ay    = [.15 .85 .15];
Y     = [.1 1];
lets  = {'(a)','(b)','(c)'};
%%
%lcols = {'k','b','r','g'};
lcols = {'k','r','b','m'};
%%
if 0
   for j=1:2
      subplot(1,2,j);
      GEN_set_plotorder;
      %%
      for r = 1:nc
         %Yplot = [m_y(:,r,j),NemNas(:,r,j),Ach_Li(:,r,j),differential(:,r,j),selfcon(:,r,j)];
         Yplot = [m_y(:,r,j),NemNas(:,r,j),Ach_Li(:,r,j),DD(:,r,j)];
         plot(xx{j},Yplot,lcols{r});
         hold on;
         plot(xx{j},differential(:,r,j),['^',lcols{r}]);
         plot(xx{j},selfcon(:,r,j),['x',lcols{r}]);
      end
   %  plot(xx{j},Ach_Li(:,:,j),'r');
   %  hold on;
   %  plot(xx{j},NemNas(:,:,j),'g');
   %  plot(xx{j},m_y(:,:,j),'k');
      %%
      X     = XX{j};
      x     = X(1)+ax(j)*(X(2)-X(1));
      y     = Y(1)+ay(j)*(Y(2)-Y(1));
      txt   = text(x,y,lets{j});
      set(txt,'FontName','times','FontSize',18);
      %%
      GEN_proc_fig(xlab{j},'my');
      xlim(X);
      ylim(Y);
   end
   saveas(gcf,[filename,'.eps'],'epsc');
else
   ar = 2.20;
   for j=1:2
      subplot(1,3,j);
      GEN_set_plotorder;
      %%
      for r = 1:nc
         %Yplot = [m_y(:,r,j),NemNas(:,r,j),Ach_Li(:,r,j),differential(:,r,j),selfcon(:,r,j)];
         %Yplot = [m_y(:,r,j),NemNas(:,r,j),Ach_Li(:,r,j),DD(:,r,j)];
         Yplot = [m_y(:,r,j),NemNas(:,r,j),Ach_Li(:,r,j)];
         plot(xx{j},Yplot,lcols{r});
         hold on;
         plot(xx{j},differential(:,r,j),['^',lcols{r}]);
         plot(xx{j},selfcon(:,r,j),['x',lcols{r}]);
      end
   %  plot(xx{j},Ach_Li(:,:,j),'r');
   %  hold on;
   %  plot(xx{j},NemNas(:,:,j),'g');
   %  plot(xx{j},m_y(:,:,j),'k');
      %%
      X     = XX{j};
      x     = X(1)+ax(j)*(X(2)-X(1));
      y     = Y(1)+ay(j)*(Y(2)-Y(1));
      txt   = text(x,y,lets{j});
      set(txt,'FontName','times','FontSize',18);
      %%
      if j==2
         set(gca,'XTick',[.2 .5 .8]);
      end
      GEN_proc_fig(xlab{j},'my');
      xlim(X);
      ylim(Y);
      GEN_daspect(ar);
   end
   %%
   %BB = [1 2];
   BB = 1;
   colin = {'k','b'};
   Ktop  = 1.25;
   for j=1:length(BB)
      B  = BB(j);
      cc = colin{j};
      %%
      A        = 1;
      a        = 0.5;
      beta     = 2*A/pi/B*log(cos(pi*a/2/A));
      my_eff   = 1/(1-beta);
      %beta     = A/pi/B*log(cos(pi*a/2/A));
      %my_eff   = 1/(1-2*beta);
      %%
      Irr_vars = CURVEget_strtline(-a,a);
      C_eff    = feval(crk_solver,[A B],Irr_vars,15);
      my_eff2  = C_eff{1}(4);
      %%
      np2      = 500;
      K        = (0:np2)'/np2*(Ktop*pi);%%=2kB
      cQ       = cos(K)+(beta/2)*K.*sin(K);
      %cQ       = cos(K)+beta*K.*sin(K);

      %Qeff_al  = acos(cQ);
      Qeff_al  = GEN_acos(cQ);%%more accurate than acos.m near cQ=1;
      %%
      %Kzero = RTS_imag_root_wtr_matlab(beta,1,pi*.4,pi*.6)
      %Qeff_al(find(imag(Qeff_al)~=0))  = NaN;
      [(Qeff_al(2)-Qeff_al(1))/(K(2)-K(1)),1/sqrt(my_eff),1/sqrt(my_eff2),sqrt(1-2*beta)]
      %%
      subplot(1,3,3);
      plot(Qeff_al/pi,K/pi,cc);
      hold on;
      %plot(0.5,Kzero/pi,'or');
      jj = find(K<.95*pi);
      %plot(sqrt(1/my_eff)*K(jj)/pi,K(jj)/pi,['--',cc]);
      plot(sqrt(1/my_eff2)*K(jj)/pi,K(jj)/pi,['--r']);
      %plot(sqrt(1/my_eff)*K(jj)/pi,K(jj)/pi,['--r']);
%     %plot(K/pi,[cQ cos(Qeff_al) 0*K]);
%     plot(K/pi,[cos(Qeff_al) 0*K]);
%     hold on, plot(Kzero/pi,0,'or');
%     hold off;
   end
   X  = [0 1];
   Y  = [0 Ktop];
   axis([X Y]);
   %%
   x  = ax(3)*X(2);
   y  = ay(3)*Y(2);
   txt   = text(x,y,lets{3});
   set(txt,'FontName','times','FontSize',18);
   GEN_proc_fig('2qB/p','2kB/p');
   GEN_daspect(ar);
   GEN_setsize_eps(.15,20,[]);
   %%
   figname  = [filename,'-vkeff.eps']
   saveas(gcf,figname,'epsc');
   %!gv out/comp-achli-nemnas-vkeff.eps &
end
