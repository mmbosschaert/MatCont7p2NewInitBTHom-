function out = PCTest
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_A,par_B,par_C)
dydt=[-kmrgd(2)-kmrgd(3);
kmrgd(1)+par_A*kmrgd(2);
par_B*kmrgd(1)-par_C*kmrgd(3)+kmrgd(1)*kmrgd(3);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(PCTest);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_A,par_B,par_C)
jac=[ 0 , -1 , -1 ; 1 , par_A , 0 ; kmrgd(3) + par_B , 0 , kmrgd(1) - par_C ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_A,par_B,par_C)
jacp=[ 0 , 0 , 0 ; kmrgd(2) , 0 , 0 ; 0 , kmrgd(1) , -kmrgd(3) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_A,par_B,par_C)
hess1=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 1 ];
hess2=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 1 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_A,par_B,par_C)
hessp1=[ 0 , 0 , 0 ; 0 , 1 , 0 ; 0 , 0 , 0 ];
hessp2=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 1 , 0 , 0 ];
hessp3=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , -1 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_A,par_B,par_C)
tens31=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens32=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,2,1) =tens34;
tens3(:,:,2,2) =tens35;
tens3(:,:,2,3) =tens36;
tens3(:,:,3,1) =tens37;
tens3(:,:,3,2) =tens38;
tens3(:,:,3,3) =tens39;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_A,par_B,par_C)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_A,par_B,par_C)
