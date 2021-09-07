function out = Bazykin
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

function dydt = fun_eval(t,in2,alpha,delta)
dydt = [in2(1,:)-in2(1,:).^2./1.0e+2-(in2(1,:).*in2(2,:))./(alpha.*in2(1,:)+1.0);-in2(2,:)-delta.*in2(2,:).^2+(in2(1,:).*in2(2,:))./(alpha.*in2(1,:)+1.0)];

function [tspan,y0,options] = init
handles = feval(Bazykin);
y0=[0;0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

function jac = jacobian(t,in2,alpha,delta)
jac = reshape([in2(1,:).*(-1.0./5.0e+1)-in2(2,:)./(alpha.*in2(1,:)+1.0)+alpha.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2+1.0,in2(2,:)./(alpha.*in2(1,:)+1.0)-alpha.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2,-in2(1,:)./(alpha.*in2(1,:)+1.0),delta.*in2(2,:).*-2.0+in2(1,:)./(alpha.*in2(1,:)+1.0)-1.0],[2,2]);

function jacp = jacobianp(t,in2,alpha,delta)
jacp = reshape([in2(1,:).^2.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2,-in2(1,:).^2.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2,0.0,-in2(2,:).^2],[2,2]);

function hess = hessians(t,in2,alpha,delta)
hess(:,:,1) =reshape([alpha.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2.*2.0-alpha.^2.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0-1.0./5.0e+1,alpha.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2.*-2.0+alpha.^2.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,-1.0./(alpha.*in2(1,:)+1.0)+alpha.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^2,1.0./(alpha.*in2(1,:)+1.0)-alpha.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^2],[2,2]);
hess(:,:,2) =reshape([-1.0./(alpha.*in2(1,:)+1.0)+alpha.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^2,1.0./(alpha.*in2(1,:)+1.0)-alpha.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^2,0.0,delta.*-2.0],[2,2]);

function hessp = hessiansp(t,in2,alpha,delta)
hessp(:,:,1) =reshape([in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2.*2.0-alpha.*in2(1,:).^2.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^2.*-2.0+alpha.*in2(1,:).^2.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,in2(1,:).^2.*1.0./(alpha.*in2(1,:)+1.0).^2,-in2(1,:).^2.*1.0./(alpha.*in2(1,:)+1.0).^2],[2,2]);
hessp(:,:,2) =reshape([0.0,0.0,0.0,in2(2,:).*-2.0],[2,2]);

function tens3 = der3(t,in2,alpha,delta)
tens3(:,:,1,1) =reshape([alpha.^2.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*-6.0+alpha.^3.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,alpha.^2.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*6.0-alpha.^3.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,alpha.*1.0./(alpha.*in2(1,:)+1.0).^2.*2.0-alpha.^2.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,alpha.*1.0./(alpha.*in2(1,:)+1.0).^2.*-2.0+alpha.^2.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0],[2,2]);
tens3(:,:,1,2) =reshape([alpha.*1.0./(alpha.*in2(1,:)+1.0).^2.*2.0-alpha.^2.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,alpha.*1.0./(alpha.*in2(1,:)+1.0).^2.*-2.0+alpha.^2.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,0.0,0.0],[2,2]);
tens3(:,:,2,1) =reshape([alpha.*1.0./(alpha.*in2(1,:)+1.0).^2.*2.0-alpha.^2.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,alpha.*1.0./(alpha.*in2(1,:)+1.0).^2.*-2.0+alpha.^2.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^3.*2.0,0.0,0.0],[2,2]);
tens3(:,:,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);

function tens4 = der4(t,in2,alpha,delta)
tens4(:,:,1,1,1) =reshape([alpha.^3.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*2.4e+1-alpha.^4.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^3.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*-2.4e+1+alpha.^4.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*-6.0+alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*6.0-alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0],[2,2]);
tens4(:,:,1,1,2) =reshape([alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*-6.0+alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*6.0-alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,0.0,0.0],[2,2]);
tens4(:,:,1,2,1) =reshape([alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*-6.0+alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*6.0-alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,0.0,0.0],[2,2]);
tens4(:,:,1,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,1,1) =reshape([alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*-6.0+alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,alpha.^2.*1.0./(alpha.*in2(1,:)+1.0).^3.*6.0-alpha.^3.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^4.*6.0,0.0,0.0],[2,2]);
tens4(:,:,2,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens4(:,:,2,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);

function tens5 = der5(t,in2,alpha,delta)
tens5(:,:,1,1,1,1) =reshape([alpha.^4.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*-1.2e+2+alpha.^5.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^6.*1.2e+2,alpha.^4.*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*1.2e+2-alpha.^5.*in2(1,:).*in2(2,:).*1.0./(alpha.*in2(1,:)+1.0).^6.*1.2e+2,alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*2.4e+1-alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*-2.4e+1+alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1],[2,2]);
tens5(:,:,1,1,1,2) =reshape([alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*2.4e+1-alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*-2.4e+1+alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,0.0,0.0],[2,2]);
tens5(:,:,1,1,2,1) =reshape([alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*2.4e+1-alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*-2.4e+1+alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,0.0,0.0],[2,2]);
tens5(:,:,1,1,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,1,1) =reshape([alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*2.4e+1-alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*-2.4e+1+alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,0.0,0.0],[2,2]);
tens5(:,:,1,2,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,1,2,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,1,1) =reshape([alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*2.4e+1-alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,alpha.^3.*1.0./(alpha.*in2(1,:)+1.0).^4.*-2.4e+1+alpha.^4.*in2(1,:).*1.0./(alpha.*in2(1,:)+1.0).^5.*2.4e+1,0.0,0.0],[2,2]);
tens5(:,:,2,1,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,1,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,1,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,1,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,2,1) =reshape([0.0,0.0,0.0,0.0],[2,2]);
tens5(:,:,2,2,2,2) =reshape([0.0,0.0,0.0,0.0],[2,2]);
