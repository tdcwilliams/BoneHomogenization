%%fig_alp_crit_rot_strt.m
fname = 'out/alp_crit_cos_phase.mat';
if exist(fname)
   load(fname);
else
   aoA      = .5;
   boB      = .5;
   Nr       = 50;
   phi_vec  = linspace(0,1,Nr)';
   alpc_vec = 0*phi_vec;
   fmin_vec = 0*phi_vec;

   for j=1
      psi                  = phi_vec(j);
      [alpc,fvalue,Mhat]   = BONE_alpcrit_cos_phase(aoA,boB,psi);
      alpc_vec(j)          = alpc;
      fmin_vec(j)          = fvalue;
   end

   for j=2:Nr
      psi                  = phi_vec(j);
      [alpc,fvalue,Mhat]   = BONE_alpcrit_cos_phase(aoA,boB,psi,alpc);
      alpc_vec(j)          = alpc;
      fmin_vec(j)          = fvalue;
      GEN_progrep([j Nr]);
   end
   save(fname,'alpc_vec','fmin_vec','phi_vec');
end

[ax,hh{1},hh{2}]  = plotyy(phi_vec,alpc_vec,phi_vec,fmin_vec);
if BW==0
   lc = {[0 0 0],[1 0 0]};
else
   lc = {[0 0 0],.35*[1 1 1]};
end
ls = {'-','--'};
lw = {1,1.5};

if 0
   yl = {'ac','Del'};
else
   yl = {'','Del'};
end

for j=1:2
   yyl   = get(ax(j),'Ylabel');
   set(yyl,'String',yl{j},'fontname','times','fontsize',16);
   if j==1
      xxl2  = get(ax(j),'Xlabel');
      set(xxl2,'String','psi','fontname','times','fontsize',16);
      xlp2  = get(xxl2,'position');
      set(ax(j),'xtick',[0 .5 1]);
   else
      set(ax(j),'xtick',[]);
      yyl2  = yyl;
      ylp2  = get(yyl,'position');
   end
   %set(get(ax(j),'fontname'),'String','times');
   set(ax(j),'fontname','times','fontsize',16,'ycolor',lc{j},...
      'xlim',[0 1]);
   %set(han,'string');
   %set(,'fontname','times','fontsize',16);
   %set(ax(j),'FontName','times','fontsize',16,...
      %'xlabel','theta, deg','ylabel',yl{j});
   set(hh{j},'color',lc{j},'linestyle',ls{j},'linewidth',lw{j});
end
