function [ip_Pn,iplv_0diff,...
  iplv_1diff,xtraP]=...
    BONE_Pn_Stuff(soL,ip1_Tn,Nterms);
%% *soL=s/L gives the quadrature points for the T_n chebyshev poly's;
%% - this function aims to do integrals with involving P_n's using these points
%%   instead of using the natural quad pts for the P_n (Legendre) poly's;
%% *if f(s)=a_nT_n(s/L), ip1_Tn*an=(1/L)*\int_{-L}^L[f(s)]ds
%% - ie does integral using chebyshev pts instead of Legendre points;

NgP                  = Nterms-1;
[ip_Pn,hnP,Pn_vals]  = ...
  OP_inprod_legendre(soL,ip1_Tn',NgP);
iplv_0diff  = diag(hnP)*ip_Pn;
%% NB needs to be multiplied by L later;
%%  - if you want to expand a function F=a_n/(h_nL)*P_n(s/L)
%%
tt       = [-1;1];
Pn_ends  = OP_interp_legendre(tt,{NgP})';
xtraP    = {Pn_vals,Pn_ends,hnP};

%% int_{-L}^L.P_m(s/L)f_s(s).ds
%%  =P_m(1).f(L)-P_m(-1).f(-L)
%%   -1/L.int_{-L}^L.P'_m(s/L)f_s(s).ds;
%% P'_m=C_{m-1}^(1.5);
alp         = 1.5;
dPn_vals    = [0*soL,...
   OP_interp_gegenbauer(soL,alp,{NgP-1})];
iplv_1diff  = -dPn_vals'*diag(ip1_Tn);
%% if f(s)=a_nT_n(s/L), iplv_1diff*an=(1/L)*\int_{-L}^L[-P'_m(s/L)f(s)]ds
