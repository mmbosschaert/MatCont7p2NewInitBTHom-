function out = BogdanovTakensNormalForm
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = @der4;
out{9} = @der5;
out{10} = [];

function dydt = fun_eval(t,in2,beta1,beta2)
dydt = [in2(2,:);beta1+beta2.*in2(2,:)+in2(1,:).*in2(2,:)-in2(1,:).^2];

function [tspan,y0,options] = init
handles = feval(BogdanovTakensNormalForm);
y0=[0;0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

function jac = jacobian(t,in2,beta1,beta2)
jac = reshape([0.0,in2(1,:).*-2.0+in2(2,:),1.0,beta2+in2(1,:)],[2,2]);

function jacp = jacobianp(t,in2,beta1,beta2)
jacp = reshape([0.0,1.0,0.0,in2(2,:)],[2,2]);

function hess = hessians(t,in2,beta1,beta2)
hess(:,:,1) =reshape([0.0,-2.0,0.0,1.0],[2,2]);
hess(:,:,2) =reshape([0.0,1.0,0.0,0.0],[2,2]);

function hessp = hessiansp(t,in2,beta1,beta2)
hessp(:,:,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
hessp(:,:,2) =reshape([0.0,0.0,0.0,1.0],[2,2]);

function tens3 = der3(t,in2,beta1,beta2)
tens3(:,:,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens3(:,:,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens3(:,:,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens3(:,:,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);

function tens4 = der4(t,in2,beta1,beta2)
tens4(:,:,1,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,1,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,1,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,1,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);

function tens5 = der5(t,in2,beta1,beta2)
tens5(:,:,1,1,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,1,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,1,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,1,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
