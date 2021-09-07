function out = ppmodel
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = @der4;
out{9} = @der5;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,a,b,epsil,d)
dydt=[a*kmrgd(1)*(1-kmrgd(1))-b*kmrgd(1)*kmrgd(2)/(1+epsil*kmrgd(1));
d*kmrgd(1)*kmrgd(2)/(1+epsil*kmrgd(1));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(ppmodel);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,a,b,epsil,d)
jac=[ (b*epsil*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 - a*(kmrgd(1) - 1) - (b*kmrgd(2))/(epsil*kmrgd(1) + 1) - a*kmrgd(1) , -(b*kmrgd(1))/(epsil*kmrgd(1) + 1) ; (d*kmrgd(2))/(epsil*kmrgd(1) + 1) - (d*epsil*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 , (d*kmrgd(1))/(epsil*kmrgd(1) + 1) ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,a,b,epsil,d)
jacp=[ -kmrgd(1)*(kmrgd(1) - 1) , -(kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1) , (b*kmrgd(1)^2*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 , 0 ; 0 , 0 , -(d*kmrgd(1)^2*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 , (kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,a,b,epsil,d)
hess1=[ (2*b*epsil*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 - 2*a - (2*b*epsil^2*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^3 , (b*epsil*kmrgd(1))/(epsil*kmrgd(1) + 1)^2 - b/(epsil*kmrgd(1) + 1) ; (2*d*epsil^2*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^3 - (2*d*epsil*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 , d/(epsil*kmrgd(1) + 1) - (d*epsil*kmrgd(1))/(epsil*kmrgd(1) + 1)^2 ];
hess2=[ (b*epsil*kmrgd(1))/(epsil*kmrgd(1) + 1)^2 - b/(epsil*kmrgd(1) + 1) , 0 ; d/(epsil*kmrgd(1) + 1) - (d*epsil*kmrgd(1))/(epsil*kmrgd(1) + 1)^2 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,a,b,epsil,d)
hessp1=[ 1 - 2*kmrgd(1) , 0 ; 0 , 0 ];
hessp2=[ (epsil*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 - kmrgd(2)/(epsil*kmrgd(1) + 1) , -kmrgd(1)/(epsil*kmrgd(1) + 1) ; 0 , 0 ];
hessp3=[ (2*b*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 - (2*b*epsil*kmrgd(1)^2*kmrgd(2))/(epsil*kmrgd(1) + 1)^3 , (b*kmrgd(1)^2)/(epsil*kmrgd(1) + 1)^2 ; (2*d*epsil*kmrgd(1)^2*kmrgd(2))/(epsil*kmrgd(1) + 1)^3 - (2*d*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 , -(d*kmrgd(1)^2)/(epsil*kmrgd(1) + 1)^2 ];
hessp4=[ 0 , 0 ; kmrgd(2)/(epsil*kmrgd(1) + 1) - (epsil*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^2 , kmrgd(1)/(epsil*kmrgd(1) + 1) ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,a,b,epsil,d)
tens31=[ (6*b*epsil^3*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^4 - (6*b*epsil^2*kmrgd(2))/(epsil*kmrgd(1) + 1)^3 , (2*b*epsil)/(epsil*kmrgd(1) + 1)^2 - (2*b*epsil^2*kmrgd(1))/(epsil*kmrgd(1) + 1)^3 ; (6*d*epsil^2*kmrgd(2))/(epsil*kmrgd(1) + 1)^3 - (6*d*epsil^3*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^4 , (2*d*epsil^2*kmrgd(1))/(epsil*kmrgd(1) + 1)^3 - (2*d*epsil)/(epsil*kmrgd(1) + 1)^2 ];
tens32=[ (2*b*epsil)/(epsil*kmrgd(1) + 1)^2 - (2*b*epsil^2*kmrgd(1))/(epsil*kmrgd(1) + 1)^3 , 0 ; (2*d*epsil^2*kmrgd(1))/(epsil*kmrgd(1) + 1)^3 - (2*d*epsil)/(epsil*kmrgd(1) + 1)^2 , 0 ];
tens33=[ (2*b*epsil)/(epsil*kmrgd(1) + 1)^2 - (2*b*epsil^2*kmrgd(1))/(epsil*kmrgd(1) + 1)^3 , 0 ; (2*d*epsil^2*kmrgd(1))/(epsil*kmrgd(1) + 1)^3 - (2*d*epsil)/(epsil*kmrgd(1) + 1)^2 , 0 ];
tens34=[ 0 , 0 ; 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,a,b,epsil,d)
tens41=[ (24*b*epsil^3*kmrgd(2))/(epsil*kmrgd(1) + 1)^4 - (24*b*epsil^4*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^5 , (6*b*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 - (6*b*epsil^2)/(epsil*kmrgd(1) + 1)^3 ; (24*d*epsil^4*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^5 - (24*d*epsil^3*kmrgd(2))/(epsil*kmrgd(1) + 1)^4 , (6*d*epsil^2)/(epsil*kmrgd(1) + 1)^3 - (6*d*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 ];
tens42=[ (6*b*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 - (6*b*epsil^2)/(epsil*kmrgd(1) + 1)^3 , 0 ; (6*d*epsil^2)/(epsil*kmrgd(1) + 1)^3 - (6*d*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens43=[ (6*b*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 - (6*b*epsil^2)/(epsil*kmrgd(1) + 1)^3 , 0 ; (6*d*epsil^2)/(epsil*kmrgd(1) + 1)^3 - (6*d*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens44=[ 0 , 0 ; 0 , 0 ];
tens45=[ (6*b*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 - (6*b*epsil^2)/(epsil*kmrgd(1) + 1)^3 , 0 ; (6*d*epsil^2)/(epsil*kmrgd(1) + 1)^3 - (6*d*epsil^3*kmrgd(1))/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens46=[ 0 , 0 ; 0 , 0 ];
tens47=[ 0 , 0 ; 0 , 0 ];
tens48=[ 0 , 0 ; 0 , 0 ];
tens4(:,:,1,1,1) =tens41;
tens4(:,:,1,1,2) =tens42;
tens4(:,:,1,2,1) =tens43;
tens4(:,:,1,2,2) =tens44;
tens4(:,:,2,1,1) =tens45;
tens4(:,:,2,1,2) =tens46;
tens4(:,:,2,2,1) =tens47;
tens4(:,:,2,2,2) =tens48;
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,a,b,epsil,d)
tens51=[ (120*b*epsil^5*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^6 - (120*b*epsil^4*kmrgd(2))/(epsil*kmrgd(1) + 1)^5 , (24*b*epsil^3)/(epsil*kmrgd(1) + 1)^4 - (24*b*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 ; (120*d*epsil^4*kmrgd(2))/(epsil*kmrgd(1) + 1)^5 - (120*d*epsil^5*kmrgd(1)*kmrgd(2))/(epsil*kmrgd(1) + 1)^6 , (24*d*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 - (24*d*epsil^3)/(epsil*kmrgd(1) + 1)^4 ];
tens52=[ (24*b*epsil^3)/(epsil*kmrgd(1) + 1)^4 - (24*b*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 , 0 ; (24*d*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 - (24*d*epsil^3)/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens53=[ (24*b*epsil^3)/(epsil*kmrgd(1) + 1)^4 - (24*b*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 , 0 ; (24*d*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 - (24*d*epsil^3)/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens54=[ 0 , 0 ; 0 , 0 ];
tens55=[ (24*b*epsil^3)/(epsil*kmrgd(1) + 1)^4 - (24*b*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 , 0 ; (24*d*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 - (24*d*epsil^3)/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens56=[ 0 , 0 ; 0 , 0 ];
tens57=[ 0 , 0 ; 0 , 0 ];
tens58=[ 0 , 0 ; 0 , 0 ];
tens59=[ (24*b*epsil^3)/(epsil*kmrgd(1) + 1)^4 - (24*b*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 , 0 ; (24*d*epsil^4*kmrgd(1))/(epsil*kmrgd(1) + 1)^5 - (24*d*epsil^3)/(epsil*kmrgd(1) + 1)^4 , 0 ];
tens510=[ 0 , 0 ; 0 , 0 ];
tens511=[ 0 , 0 ; 0 , 0 ];
tens512=[ 0 , 0 ; 0 , 0 ];
tens513=[ 0 , 0 ; 0 , 0 ];
tens514=[ 0 , 0 ; 0 , 0 ];
tens515=[ 0 , 0 ; 0 , 0 ];
tens516=[ 0 , 0 ; 0 , 0 ];
tens5(:,:,1,1,1,1) =tens51;
tens5(:,:,1,1,1,2) =tens52;
tens5(:,:,1,1,2,1) =tens53;
tens5(:,:,1,1,2,2) =tens54;
tens5(:,:,1,2,1,1) =tens55;
tens5(:,:,1,2,1,2) =tens56;
tens5(:,:,1,2,2,1) =tens57;
tens5(:,:,1,2,2,2) =tens58;
tens5(:,:,2,1,1,1) =tens59;
tens5(:,:,2,1,1,2) =tens510;
tens5(:,:,2,1,2,1) =tens511;
tens5(:,:,2,1,2,2) =tens512;
tens5(:,:,2,2,1,1) =tens513;
tens5(:,:,2,2,1,2) =tens514;
tens5(:,:,2,2,2,1) =tens515;
tens5(:,:,2,2,2,2) =tens516;
