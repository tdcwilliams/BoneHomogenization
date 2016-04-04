function F=BONE_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre)
%% CALL: F=BONE_multipoles_laplace_doubly_periodic(...
%%          X,Y,M,AB,centre)
%%
%% X,Y are x&y coordinates of points multipoles are to be
%%  eval'd at;
%% M is no of multipoles required;
%% centre=[c1,c2] contains coordinates of point
%%  where multipoles are to be centred about;
%%
%% AB=[A,B] are 2 parameters needed
%%  to describe the rectangular unit cell:
%%  - periodicity is 2A in x dirn
%%  & periodicity is 2B in y dirn;

DO_TEST=0;
if nargin==0%%do test
  DO_TEST   = 1;
  XTEST     = 1;
  AB        = [.5 .5];
  M         = 6;
  centre    = [0 0];
  %%
  Nint      = 200;
  [tq,wq]   = GEN_numint_exp(Nint);
  if XTEST
    X = tq*AB(1);
    Y = .2+0*X;
  else
    Y = tq*AB(2);
    X = .2+0*Y;
  end
end

Z  = X-centre(1)+1i*(Y-centre(2));
np = length(Z);
%%
F     = zeros(np,M);
A     = AB(1);
B     = AB(2);
alp   = pi/(2*A);
bet   = pi/(2*B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0%%do infinite sum by checking the error 1st;
   jj       = (1:M);
   tol      = 1e-26;
   Nsum     = ceil(-log(tol)./(2*pi*jj*B/A));
   Nsum0    = ceil(-log(tol)/(2*pi*A/B));
   Nsum(1)  = max(Nsum(1),Nsum0);
   Nsum     = Nsum+5;

   for k=Nsum(1):-1:1
      Zp     = Z+2i*k*B;
      Zm     = Z-2i*k*B;
      Zi     = -1i*Z;
      F(:,1) = F(:,1)+2*alp*real( sin(2*alp*Z)./...
                    ( cosh(2*k*pi*B/A)-cos(2*alp*Z) ) )+...
                   +2i*bet*real( sin(2*bet*Zi)./...
                    ( cosh(2*k*pi*A/B)-cos(2*bet*Zi) ) );
      for j=2:2:M
         if k<=Nsum(1+j/2)
            F(:,j)   = F(:,j)+alp^j./sin(alp*Zp).^j+...
                         +alp^j./sin(alp*Zm).^j;
         end
      end
      for j=3:2:M
    %  [cos(alp*Zp),sin(alp*Zp).^j],pause
         if k<=Nsum((j+1)/2)
            F(:,j)   = F(:,j)+alp^j.*cos(alp*Zp)./sin(alp*Zp).^j+...
                         +alp^j.*cos(alp*Zm)./sin(alp*Zm).^j;
         end
      end
   end
else%%use while loop
   k        = 0;
   critter  = 1;
   tol      = 1e-8;
   DO_PROD  = 0;
   if DO_PROD==1
      expF  = exp(F);
   end
   while critter
      k     = k+1;
      Zp    = Z+2i*k*B;
      Zm    = Z-2i*k*B;
      Zi    = -1i*Z;
      dF_j  = 2*alp*real( sin(2*alp*Z)./...
                ( cosh(2*k*pi*B/A)-cos(2*alp*Z) ) )+...
              +2i*bet*real( sin(2*bet*Zi)./...
                ( cosh(2*k*pi*A/B)-cos(2*bet*Zi) ) );
      if DO_PROD==0
         F(:,1)   = F(:,1)+dF_j;
      else
         expF(:,1)   = expF(:,1).*exp(dF_j);
      end
      critter  = (max(abs(dF_j))>tol);
      %%
      for j=2:2:M
         dF_j     = alp^j./sin(alp*Zp).^j+...
                     +alp^j./sin(alp*Zm).^j;
         if DO_PROD==0
            F(:,j)   = F(:,j)+dF_j;
         else
            expF(:,j)   = expF(:,j).*exp(dF_j);
         end
         critter  = critter|(max(abs(dF_j))>tol);
      end
      for j=3:2:M
    %  [cos(alp*Zp),sin(alp*Zp).^j],pause
         dF_j     = alp^j.*cos(alp*Zp)./sin(alp*Zp).^j+...
                      +alp^j.*cos(alp*Zm)./sin(alp*Zm).^j;
         if DO_PROD==0
            F(:,j)   = F(:,j)+dF_j;
         else
            expF(:,j)   = expF(:,j).*exp(dF_j);
         end
         critter  = critter|(max(abs(dF_j))>tol);
      end
   end
   if DO_PROD
      F  = log(expF);
   end
end%%of infinite sum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(:,1)   = F(:,1)+alp*real(cot(alp*Z)) +1i*bet*real(cot(bet*Zi));
for j=2:2:M
  F(:,j) = F(:,j)+alp^j./sin(alp*Z).^j;
end
for j=3:2:M
  F(:,j) = F(:,j)+alp^j.*cos(alp*Z)./sin(alp*Z).^j;
end

if DO_TEST
  [IP,hn,CS]      = GEN_inprod_cos_sin(tq,wq,Nint-1);
  [dF_dx,dF_dy]   = ...
    BONE_diff_multipoles_laplace_doubly_periodic(...
             X,Y,M,AB,centre);

  if XTEST%%test x-periodicity & derivatives:
    kn                  = (1:Nint-1)*pi/A;
    DD                  = diag(kn);
    Diff                = zeros(2*Nint-1,2*Nint-1);
    Diff(2:end,2:end)   = [0*DD,DD;-DD,0*DD];
    disp('testing x variations:');

    for m=1:M
      m
      %% test periodicity:
      x     = [X;X+2*A];
      f     = [F(:,m);F(:,m)];
      fn    = IP*F(:,m);
      fap   = CS*fn;
      %%
      subplot(2,2,1);
      plot(x,real(f)), hold on;
      plot(X,real(fap),'--r'), hold off;
      %%
      subplot(2,2,2);
      plot(x,imag(f)), hold on;
      plot(X,imag(fap),'--r'), hold off;
      %%
      df    = [dF_dx(:,m);dF_dx(:,m)];
      dfn   = Diff*fn;
      dfap  = CS*dfn;
      %%
      subplot(2,2,3);
      plot(x,real(df)), hold on;
      plot(X,real(dfap),'--r'), hold off;
      %%
      subplot(2,2,4);
      plot(x,imag(df)), hold on;
      plot(X,imag(dfap),'--r'), hold off;
      pause;
    end
  else
    kn                  = (1:Nint-1)*pi/B;
    DD                  = diag(kn);
    Diff                = zeros(2*Nint-1,2*Nint-1);
    Diff(2:end,2:end)   = [0*DD,DD;-DD,0*DD];
    disp('testing y variations:');

    for m=1:M
      m
      %% test periodicity:
      y     = [Y;Y+2*B];
      f     = [F(:,m);F(:,m)];
      fn    = IP*F(:,m);
      fap   = CS*fn;
      %%
      subplot(2,2,1);
      plot(y,real(f)), hold on;
      plot(Y,real(fap),'--r'), hold off;
      %%
      subplot(2,2,2);
      plot(y,imag(f)), hold on;
      plot(Y,imag(fap),'--r'), hold off;
      %%
      df    = [dF_dy(:,m);dF_dy(:,m)];
      dfn   = Diff*fn;
      dfap  = CS*dfn;
      %%
      subplot(2,2,3);
      plot(y,real(df)), hold on;
      plot(Y,real(dfap),'--r'), hold off;
      %%
      subplot(2,2,4);
      plot(y,imag(df)), hold on;
      plot(Y,imag(dfap),'--r'), hold off;
      pause;
    end
  end
end
