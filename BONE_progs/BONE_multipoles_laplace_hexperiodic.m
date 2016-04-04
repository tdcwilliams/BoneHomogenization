function F=BONE_multipoles_laplace_hexperiodic(...
            X,Y,M,BDrho,centre)
%% CALL: F=BONE_multipoles_laplace_hexperiodic(...
%%          X,Y,M,BDrho,centre)
%%
%% X,Y are x&y coordinates of points multipoles are to be
%%  eval'd at;
%% M is no of multipoles required;
%% centre=[c1,c2] contains coordinates of point
%%  where multipoles are to be centred about;
%%
%% BDrho=[B,D,rho=E/D] are three parameters needed
%%  to describe the unit cell;
%%
%% hex unit cell - top side and bottom side are
%% located at y=\pm B,
%%  & are both 2E long;
%% unit cell is symmetric about y=0,
%%  & is 2D wide;
%% -periodicity in X direction is 2A, A=D+E;
%% -periodicity in Y direction is 2B;
%% -there is periodicity also at an angle \phi=atan(B/A) direction
%%   & repeating length is C=sqrt(A^2+B^2);

do_test  = 0;
if nargin==0%%check look of unit cell;
   do_test  = 1;
   if 1%%regular hexagon, sides of unit length;
      B     = sqrt(3)/2;
      D     = 1;
      rho   = .5/D;
      BDrho = [B D rho];
      %%
      if 0
         X  = 1;
         Y  = 1;
      else
         xx = linspace(-1,2.5,150);
         yy = linspace(-1,3,150);
         [YY,XX]  = meshgrid(yy,xx);
         X  = XX(:);
         Y  = YY(:);
         %{X,Y}, return;
      end
      M        = 5;
      centre   = [0 0];
   end
end
B  = BDrho(1);
D  = BDrho(2);
E  = D*BDrho(3);
%%
A  = D+E;
AB = [A,B];
F  = BONE_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre);
%%
%phi   = atan(B/A);
%C     = sqrt(B^2+A^2);
%%
if 1
   centre2(1)  = centre(1)+A;
   centre2(2)  = centre(2)+B;
   F           = F + BONE_multipoles_laplace_doubly_periodic(...
                   X,Y,M,AB,centre2);
else
   X2 = X-A;
   Y2 = Y-B;
   F  = F + BONE_multipoles_laplace_doubly_periodic(...
          X2,Y2,M,AB,centre);
end

if do_test
   for j=1:2
      subplot(1,2,j);
      FF    = 0*XX;
      FF(:) = real((-1i)^j*F(:,M));
      pcolor(XX,YY,FF);
      caxis(10*[-1 1]);
      daspect([1 1 1]);
      hold on;
      'paused',pause
      %%
      zp = exp(2i*pi/6*(0:5)');
      xp = real(zp);
      yp = imag(zp);
      plot(xp,yp,'or');
      hold on;
   %  plot([-xp(2) xp(2)],B*[1 1],'k');
   %  plot([-xp(2) xp(2)],-B*[1 1],'--k');
      lw = 1.5;
      xt = centre(1)+E*[-1 1];
      yt = centre(2)+B*[1 1];
      xb = centre(1)+E*[-1 1];
      yb = centre(2)-B*[1 1];
      xc = centre(1)+D*[-1 1];
      yc = centre(2)+0*[-1 1];
      x1 = [xt,xc(2),fliplr(xb),xc(1),xt(1)];
      y1 = [yt,yc(2),fliplr(yb),yc(1),yt(1)];
      %plot(E*[-1 1],B*[1 1],'k','linewidth',lw);
      %plot(E*[-1 1],-B*[1 1],'--k','linewidth',lw);
      %plot(D*[-1 1],0*[1 1],'r','linewidth',lw);
      plot(x1,y1,'k','linewidth',lw);
      %disp('centre cell'),pause;
      %%
      xt = centre(1)+E*[-1 1];
      yt = centre(2)+B*[1 1]+2*B;
      xb = centre(1)+E*[-1 1];
      yb = centre(2)-B*[1 1]+2*B;
      xc = centre(1)+D*[-1 1];
      yc = centre(2)+0*[-1 1]+2*B;
      x1 = [xt,xc(2),fliplr(xb),xc(1),xt(1)];
      y1 = [yt,yc(2),fliplr(yb),yc(1),yt(1)];
      %plot(E*[-1 1],3*B*[1 1],'k','linewidth',lw);
      %plot(D*[-1 1],2*B*[1 1],'r','linewidth',lw);
      plot(x1,y1,'k','linewidth',lw);
      %disp('cell above'),pause;
      %%
      xt = centre2(1)+E*[-1 1];
      yt = centre2(2)+B*[1 1];
      xb = centre2(1)+E*[-1 1];
      yb = centre2(2)-B*[1 1];
      xc = centre2(1)+D*[-1 1];
      yc = centre2(2)+0*[-1 1];
      x1 = [xt,xc(2),fliplr(xb),xc(1),xt(1)];
      y1 = [yt,yc(2),fliplr(yb),yc(1),yt(1)];
      %%
      %plot(xt,yt,'k','linewidth',lw);
      %plot(xb,yb,'k','linewidth',lw);
      %plot(xc,yc,'r','linewidth',lw);
      plot(x1,y1,'k','linewidth',lw);
      disp('cell up & to the right');
      %%
      %daspect([1 1 1]);
      hold off;
      'paused',pause
   end
end
