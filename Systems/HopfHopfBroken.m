function out = HopfHopfBroken
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
function dydt = fun_eval(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
dydt=[kmrgd(1)*(par_u1+par_p11*(kmrgd(1)^2+kmrgd(2)^2)+par_p12*(kmrgd(3)^2+kmrgd(4)^2)+par_s1*(kmrgd(3)^2+kmrgd(4)^2)^2)-kmrgd(2)*par_w1+3*kmrgd(2)^6;
kmrgd(2)*(par_u1+par_p11*(kmrgd(1)^2+kmrgd(2)^2)+par_p12*(kmrgd(3)^2+kmrgd(4)^2)+par_s1*(kmrgd(3)^2+kmrgd(4)^2)^2)+kmrgd(1)*par_w1-2*kmrgd(3)^6;
kmrgd(3)*(par_u2+par_p21*(kmrgd(1)^2+kmrgd(2)^2)+par_p22*(kmrgd(3)^2+kmrgd(4)^2)+par_s2*(kmrgd(1)^2+kmrgd(2)^2)^2)-kmrgd(4)*par_w2-7*kmrgd(2)^6;
kmrgd(4)*(par_u2+par_p21*(kmrgd(1)^2+kmrgd(2)^2)+par_p22*(kmrgd(3)^2+kmrgd(4)^2)+par_s2*(kmrgd(1)^2+kmrgd(2)^2)^2)+kmrgd(3)*par_w2+kmrgd(1)^6;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(HopfHopfBroken);
y0=[0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
jac=[ par_u1 + par_s1*(kmrgd(3)^2 + kmrgd(4)^2)^2 + 2*kmrgd(1)^2*par_p11 + par_p11*(kmrgd(1)^2 + kmrgd(2)^2) + par_p12*(kmrgd(3)^2 + kmrgd(4)^2) , 18*kmrgd(2)^5 - par_w1 + 2*kmrgd(1)*kmrgd(2)*par_p11 , kmrgd(1)*(2*kmrgd(3)*par_p12 + 4*kmrgd(3)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) , kmrgd(1)*(2*kmrgd(4)*par_p12 + 4*kmrgd(4)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) ; par_w1 + 2*kmrgd(1)*kmrgd(2)*par_p11 , par_u1 + par_s1*(kmrgd(3)^2 + kmrgd(4)^2)^2 + 2*kmrgd(2)^2*par_p11 + par_p11*(kmrgd(1)^2 + kmrgd(2)^2) + par_p12*(kmrgd(3)^2 + kmrgd(4)^2) , kmrgd(2)*(2*kmrgd(3)*par_p12 + 4*kmrgd(3)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) - 12*kmrgd(3)^5 , kmrgd(2)*(2*kmrgd(4)*par_p12 + 4*kmrgd(4)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) ; kmrgd(3)*(2*kmrgd(1)*par_p21 + 4*kmrgd(1)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) , kmrgd(3)*(2*kmrgd(2)*par_p21 + 4*kmrgd(2)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) - 42*kmrgd(2)^5 , par_u2 + par_s2*(kmrgd(1)^2 + kmrgd(2)^2)^2 + 2*kmrgd(3)^2*par_p22 + par_p21*(kmrgd(1)^2 + kmrgd(2)^2) + par_p22*(kmrgd(3)^2 + kmrgd(4)^2) , 2*kmrgd(3)*kmrgd(4)*par_p22 - par_w2 ; kmrgd(4)*(2*kmrgd(1)*par_p21 + 4*kmrgd(1)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) + 6*kmrgd(1)^5 , kmrgd(4)*(2*kmrgd(2)*par_p21 + 4*kmrgd(2)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) , par_w2 + 2*kmrgd(3)*kmrgd(4)*par_p22 , par_u2 + par_s2*(kmrgd(1)^2 + kmrgd(2)^2)^2 + 2*kmrgd(4)^2*par_p22 + par_p21*(kmrgd(1)^2 + kmrgd(2)^2) + par_p22*(kmrgd(3)^2 + kmrgd(4)^2) ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
jacp=[ kmrgd(1) , 0 , kmrgd(1)*(kmrgd(1)^2 + kmrgd(2)^2) , kmrgd(1)*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , 0 , kmrgd(1)*(kmrgd(3)^2 + kmrgd(4)^2)^2 , 0 , -kmrgd(2) , 0 ; kmrgd(2) , 0 , kmrgd(2)*(kmrgd(1)^2 + kmrgd(2)^2) , kmrgd(2)*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , 0 , kmrgd(2)*(kmrgd(3)^2 + kmrgd(4)^2)^2 , 0 , kmrgd(1) , 0 ; 0 , kmrgd(3) , 0 , 0 , kmrgd(3)*(kmrgd(1)^2 + kmrgd(2)^2) , kmrgd(3)*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , kmrgd(3)*(kmrgd(1)^2 + kmrgd(2)^2)^2 , 0 , -kmrgd(4) ; 0 , kmrgd(4) , 0 , 0 , kmrgd(4)*(kmrgd(1)^2 + kmrgd(2)^2) , kmrgd(4)*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , kmrgd(4)*(kmrgd(1)^2 + kmrgd(2)^2)^2 , 0 , kmrgd(3) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
hess1=[ 6*kmrgd(1)*par_p11 , 2*kmrgd(2)*par_p11 , 2*kmrgd(3)*par_p12 + 4*kmrgd(3)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 2*kmrgd(4)*par_p12 + 4*kmrgd(4)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) ; 2*kmrgd(2)*par_p11 , 2*kmrgd(1)*par_p11 , 0 , 0 ; kmrgd(3)*(2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) , 8*kmrgd(1)*kmrgd(3)*kmrgd(2)*par_s2 , 2*kmrgd(1)*par_p21 + 4*kmrgd(1)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 ; kmrgd(4)*(2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) + 30*kmrgd(1)^4 , 8*kmrgd(1)*kmrgd(2)*kmrgd(4)*par_s2 , 0 , 2*kmrgd(1)*par_p21 + 4*kmrgd(1)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) ];
hess2=[ 2*kmrgd(2)*par_p11 , 2*kmrgd(1)*par_p11 + 90*kmrgd(2)^4 , 0 , 0 ; 2*kmrgd(1)*par_p11 , 6*kmrgd(2)*par_p11 , 2*kmrgd(3)*par_p12 + 4*kmrgd(3)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 2*kmrgd(4)*par_p12 + 4*kmrgd(4)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) ; 8*kmrgd(1)*kmrgd(3)*kmrgd(2)*par_s2 , kmrgd(3)*(2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) - 210*kmrgd(2)^4 , 2*kmrgd(2)*par_p21 + 4*kmrgd(2)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 ; 8*kmrgd(1)*kmrgd(2)*kmrgd(4)*par_s2 , kmrgd(4)*(2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2)) , 0 , 2*kmrgd(2)*par_p21 + 4*kmrgd(2)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) ];
hess3=[ 2*kmrgd(3)*par_p12 + 4*kmrgd(3)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , kmrgd(1)*(2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) , 8*kmrgd(1)*kmrgd(3)*kmrgd(4)*par_s1 ; 0 , 2*kmrgd(3)*par_p12 + 4*kmrgd(3)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , kmrgd(2)*(2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) - 60*kmrgd(3)^4 , 8*kmrgd(3)*kmrgd(2)*kmrgd(4)*par_s1 ; 2*kmrgd(1)*par_p21 + 4*kmrgd(1)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 2*kmrgd(2)*par_p21 + 4*kmrgd(2)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 6*kmrgd(3)*par_p22 , 2*kmrgd(4)*par_p22 ; 0 , 0 , 2*kmrgd(4)*par_p22 , 2*kmrgd(3)*par_p22 ];
hess4=[ 2*kmrgd(4)*par_p12 + 4*kmrgd(4)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , 8*kmrgd(1)*kmrgd(3)*kmrgd(4)*par_s1 , kmrgd(1)*(2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) ; 0 , 2*kmrgd(4)*par_p12 + 4*kmrgd(4)*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 8*kmrgd(3)*kmrgd(2)*kmrgd(4)*par_s1 , kmrgd(2)*(2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2)) ; 0 , 0 , 2*kmrgd(4)*par_p22 , 2*kmrgd(3)*par_p22 ; 2*kmrgd(1)*par_p21 + 4*kmrgd(1)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 2*kmrgd(2)*par_p21 + 4*kmrgd(2)*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 2*kmrgd(3)*par_p22 , 6*kmrgd(4)*par_p22 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
hessp1=[ 1 , 0 , 0 , 0 ; 0 , 1 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp2=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 1 , 0 ; 0 , 0 , 0 , 1 ];
hessp3=[ 3*kmrgd(1)^2 + kmrgd(2)^2 , 2*kmrgd(1)*kmrgd(2) , 0 , 0 ; 2*kmrgd(1)*kmrgd(2) , kmrgd(1)^2 + 3*kmrgd(2)^2 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp4=[ kmrgd(3)^2 + kmrgd(4)^2 , 0 , 2*kmrgd(1)*kmrgd(3) , 2*kmrgd(1)*kmrgd(4) ; 0 , kmrgd(3)^2 + kmrgd(4)^2 , 2*kmrgd(3)*kmrgd(2) , 2*kmrgd(2)*kmrgd(4) ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp5=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 2*kmrgd(1)*kmrgd(3) , 2*kmrgd(3)*kmrgd(2) , kmrgd(1)^2 + kmrgd(2)^2 , 0 ; 2*kmrgd(1)*kmrgd(4) , 2*kmrgd(2)*kmrgd(4) , 0 , kmrgd(1)^2 + kmrgd(2)^2 ];
hessp6=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 3*kmrgd(3)^2 + kmrgd(4)^2 , 2*kmrgd(3)*kmrgd(4) ; 0 , 0 , 2*kmrgd(3)*kmrgd(4) , kmrgd(3)^2 + 3*kmrgd(4)^2 ];
hessp7=[ (kmrgd(3)^2 + kmrgd(4)^2)^2 , 0 , 4*kmrgd(1)*kmrgd(3)*(kmrgd(3)^2 + kmrgd(4)^2) , 4*kmrgd(1)*kmrgd(4)*(kmrgd(3)^2 + kmrgd(4)^2) ; 0 , (kmrgd(3)^2 + kmrgd(4)^2)^2 , 4*kmrgd(3)*kmrgd(2)*(kmrgd(3)^2 + kmrgd(4)^2) , 4*kmrgd(2)*kmrgd(4)*(kmrgd(3)^2 + kmrgd(4)^2) ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp8=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 4*kmrgd(1)*kmrgd(3)*(kmrgd(1)^2 + kmrgd(2)^2) , 4*kmrgd(3)*kmrgd(2)*(kmrgd(1)^2 + kmrgd(2)^2) , (kmrgd(1)^2 + kmrgd(2)^2)^2 , 0 ; 4*kmrgd(1)*kmrgd(4)*(kmrgd(1)^2 + kmrgd(2)^2) , 4*kmrgd(2)*kmrgd(4)*(kmrgd(1)^2 + kmrgd(2)^2) , 0 , (kmrgd(1)^2 + kmrgd(2)^2)^2 ];
hessp9=[ 0 , -1 , 0 , 0 ; 1 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp10=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -1 ; 0 , 0 , 1 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
hessp(:,:,10) =hessp10;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
tens31=[ 6*par_p11 , 0 , 0 , 0 ; 0 , 2*par_p11 , 0 , 0 ; 24*kmrgd(1)*kmrgd(3)*par_s2 , 8*kmrgd(3)*kmrgd(2)*par_s2 , 2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 ; 120*kmrgd(1)^3 + 24*kmrgd(1)*kmrgd(4)*par_s2 , 8*kmrgd(2)*kmrgd(4)*par_s2 , 0 , 2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) ];
tens32=[ 0 , 2*par_p11 , 0 , 0 ; 2*par_p11 , 0 , 0 , 0 ; 8*kmrgd(3)*kmrgd(2)*par_s2 , 8*kmrgd(1)*kmrgd(3)*par_s2 , 8*kmrgd(1)*kmrgd(2)*par_s2 , 0 ; 8*kmrgd(2)*kmrgd(4)*par_s2 , 8*kmrgd(1)*kmrgd(4)*par_s2 , 0 , 8*kmrgd(1)*kmrgd(2)*par_s2 ];
tens33=[ 0 , 0 , 2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 8*kmrgd(3)*kmrgd(4)*par_s1 ; 0 , 0 , 0 , 0 ; 2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 8*kmrgd(1)*kmrgd(2)*par_s2 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens34=[ 0 , 0 , 8*kmrgd(3)*kmrgd(4)*par_s1 , 2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 8*kmrgd(1)*kmrgd(2)*par_s2 , 0 , 0 ];
tens35=[ 0 , 2*par_p11 , 0 , 0 ; 2*par_p11 , 0 , 0 , 0 ; 8*kmrgd(3)*kmrgd(2)*par_s2 , 8*kmrgd(1)*kmrgd(3)*par_s2 , 8*kmrgd(1)*kmrgd(2)*par_s2 , 0 ; 8*kmrgd(2)*kmrgd(4)*par_s2 , 8*kmrgd(1)*kmrgd(4)*par_s2 , 0 , 8*kmrgd(1)*kmrgd(2)*par_s2 ];
tens36=[ 2*par_p11 , 360*kmrgd(2)^3 , 0 , 0 ; 0 , 6*par_p11 , 0 , 0 ; 8*kmrgd(1)*kmrgd(3)*par_s2 , 24*kmrgd(3)*kmrgd(2)*par_s2 - 840*kmrgd(2)^3 , 2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 ; 8*kmrgd(1)*kmrgd(4)*par_s2 , 24*kmrgd(2)*kmrgd(4)*par_s2 , 0 , 2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) ];
tens37=[ 0 , 0 , 0 , 0 ; 0 , 0 , 2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 8*kmrgd(3)*kmrgd(4)*par_s1 ; 8*kmrgd(1)*kmrgd(2)*par_s2 , 2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 , 0 ; 0 , 0 , 8*kmrgd(3)*kmrgd(4)*par_s1 , 2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) ; 0 , 0 , 0 , 0 ; 8*kmrgd(1)*kmrgd(2)*par_s2 , 2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 , 0 ];
tens39=[ 0 , 0 , 2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 8*kmrgd(3)*kmrgd(4)*par_s1 ; 0 , 0 , 0 , 0 ; 2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 8*kmrgd(1)*kmrgd(2)*par_s2 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens310=[ 0 , 0 , 0 , 0 ; 0 , 0 , 2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 8*kmrgd(3)*kmrgd(4)*par_s1 ; 8*kmrgd(1)*kmrgd(2)*par_s2 , 2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens311=[ 2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , 24*kmrgd(1)*kmrgd(3)*par_s1 , 8*kmrgd(1)*kmrgd(4)*par_s1 ; 0 , 2*par_p12 + 8*kmrgd(3)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 24*kmrgd(3)*kmrgd(2)*par_s1 - 240*kmrgd(3)^3 , 8*kmrgd(2)*kmrgd(4)*par_s1 ; 0 , 0 , 6*par_p22 , 0 ; 0 , 0 , 0 , 2*par_p22 ];
tens312=[ 8*kmrgd(3)*kmrgd(4)*par_s1 , 0 , 8*kmrgd(1)*kmrgd(4)*par_s1 , 8*kmrgd(1)*kmrgd(3)*par_s1 ; 0 , 8*kmrgd(3)*kmrgd(4)*par_s1 , 8*kmrgd(2)*kmrgd(4)*par_s1 , 8*kmrgd(3)*kmrgd(2)*par_s1 ; 0 , 0 , 0 , 2*par_p22 ; 0 , 0 , 2*par_p22 , 0 ];
tens313=[ 0 , 0 , 8*kmrgd(3)*kmrgd(4)*par_s1 , 2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 2*par_p21 + 8*kmrgd(1)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 8*kmrgd(1)*kmrgd(2)*par_s2 , 0 , 0 ];
tens314=[ 0 , 0 , 0 , 0 ; 0 , 0 , 8*kmrgd(3)*kmrgd(4)*par_s1 , 2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) ; 0 , 0 , 0 , 0 ; 8*kmrgd(1)*kmrgd(2)*par_s2 , 2*par_p21 + 8*kmrgd(2)^2*par_s2 + 4*par_s2*(kmrgd(1)^2 + kmrgd(2)^2) , 0 , 0 ];
tens315=[ 8*kmrgd(3)*kmrgd(4)*par_s1 , 0 , 8*kmrgd(1)*kmrgd(4)*par_s1 , 8*kmrgd(1)*kmrgd(3)*par_s1 ; 0 , 8*kmrgd(3)*kmrgd(4)*par_s1 , 8*kmrgd(2)*kmrgd(4)*par_s1 , 8*kmrgd(3)*kmrgd(2)*par_s1 ; 0 , 0 , 0 , 2*par_p22 ; 0 , 0 , 2*par_p22 , 0 ];
tens316=[ 2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 0 , 8*kmrgd(1)*kmrgd(3)*par_s1 , 24*kmrgd(1)*kmrgd(4)*par_s1 ; 0 , 2*par_p12 + 8*kmrgd(4)^2*par_s1 + 4*par_s1*(kmrgd(3)^2 + kmrgd(4)^2) , 8*kmrgd(3)*kmrgd(2)*par_s1 , 24*kmrgd(2)*kmrgd(4)*par_s1 ; 0 , 0 , 2*par_p22 , 0 ; 0 , 0 , 0 , 6*par_p22 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,1,4) =tens34;
tens3(:,:,2,1) =tens35;
tens3(:,:,2,2) =tens36;
tens3(:,:,2,3) =tens37;
tens3(:,:,2,4) =tens38;
tens3(:,:,3,1) =tens39;
tens3(:,:,3,2) =tens310;
tens3(:,:,3,3) =tens311;
tens3(:,:,3,4) =tens312;
tens3(:,:,4,1) =tens313;
tens3(:,:,4,2) =tens314;
tens3(:,:,4,3) =tens315;
tens3(:,:,4,4) =tens316;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_u1,par_u2,par_p11,par_p12,par_p21,par_p22,par_s1,par_s2,par_w1,par_w2)
