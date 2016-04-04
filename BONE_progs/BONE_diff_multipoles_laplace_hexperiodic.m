function [dF_dx,dF_dy] = BONE_diff_multipoles_laplace_hexperiodic(...
            X,Y,M,BDrho,centre)
%% CALL: [dF_dx,dF_dy] = BONE_multipoles_laplace_hexperiodic(...
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
      X        = 1;
      Y        = 1;
      M        = 1;
      centre   = [0 0];
   end
end
B  = BDrho(1);
D  = BDrho(2);
E  = D*BDrho(3);
%%
A  = D+E;
AB = [A,B];
[dF_dx,dF_dy]  = BONE_diff_multipoles_laplace_doubly_periodic(...
                  X,Y,M,AB,centre);
%%
%phi   = atan(B/A);
%C     = sqrt(B^2+A^2);
%%
if 1
   centre2(1)  = centre(1)+A;
   centre2(2)  = centre(2)+B;
   [F_x,F_y]   = BONE_diff_multipoles_laplace_doubly_periodic(...
                   X,Y,M,AB,centre2);
   dF_dx = dF_dx+F_x;
   dF_dy = dF_dy+F_y;
end

if do_test
   zp = exp(2i*pi/6*(0:5)');
   xp = real(zp);
   yp = imag(zp);
   plot(xp,yp,'or');
   hold on;
%  plot([-xp(2) xp(2)],B*[1 1],'k');
%  plot([-xp(2) xp(2)],-B*[1 1],'--k');
   plot(E*[-1 1],B*[1 1],'k');
   plot(E*[-1 1],-B*[1 1],'--k');
   plot(D*[-1 1],0*[1 1],'r');
   disp('centre cell'),pause;
   %%
   plot(E*[-1 1],3*B*[1 1],'k');
   plot(D*[-1 1],2*B*[1 1],'r');
   disp('cell above'),pause;
   %%
   plot(centre2(1)+E*[-1 1],centre2(2)+B*[1 1],'k');
   plot(centre2(1)+E*[-1 1],centre2(2)-B*[1 1],'k');
   plot(centre2(1)+D*[-1 1],centre2(2)+0*[1 1],'r');
   disp('cell up & to the right');
   %%
   daspect([1 1 1]);
   hold off;
end
