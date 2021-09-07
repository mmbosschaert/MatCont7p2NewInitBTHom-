function bt = BT_nmfm_without_e_b1(odefile, bt, ap, options)

%% Multi-linear forms
F = nmfm_deriv_define(odefile, bt, ap);
%% Compute null vectors q0, q1, p0 and p1
[q0,p1,q1,p0]=nmfm_nullvectors(F.A);
Bord = [F.A p1; q0' 0]; % Setup bordered system (should be moved).
AINV = @(x) bord_inv(Bord, x);
%% Take transpose of the adjoint vectors
p1=p1';
p0=p0';
%% Critical normal form coefficients
a = 1/2*p1*F.B(q0,q0);
b = p1*F.B(q0,q1)+p0*F.B(q0,q0);
H2000hat = -AINV(F.B(q0,q0)-2*a*q1);
H1100hat = -AINV(F.B(q0,q1)-b*q1-H2000hat);
gamma1 = p0*(F.B(q0,q1)-H2000hat)+1/2*p1*F.B(q1,q1);
gamma2 = 1/(6*a)*(p1*(2*F.B(q0,H1100hat)+F.B(H2000hat,q1)+F.C(q0,q0,q1)) ...
				+2*a*p0*F.B(q1,q1)+2*b*p0*(F.B(q0,q1)-H2000hat) ...
				+p0*(3*F.B(H2000hat,q0)+F.C(q0,q0,q0))+gamma1*b ...
				-10*a*p0*H1100hat);
H2000 = H2000hat+gamma1*q0;
H1100 = H1100hat+gamma1*q1+gamma2*q0;
H0200 = -AINV(F.B(q1,q1)-2*H1100);
d = 1/6*p1*(F.C(q0,q0,q0)+3*F.B(q0,H2000)-6*a*H1100);
H3000 = -AINV(3*F.B(H2000,q0)+F.C(q0,q0,q0)-6*d*q1-6*a*H1100);
H2100 = -AINV(-2*a*H0200-2*b*H1100-H3000+2*F.B(H1100,q0) ...
				+ F.B(H2000,q1)+F.C(q0,q0,q1));
%% H0010,K10,H0001,K01
eta = (p1*F.J1)';
K10hat = 1/(eta'*eta)*eta;
H0010hat = q0-AINV(F.J1*K10hat);
H0010hat = AINV(q1-F.J1*K10hat);
K01hat = [[0,-1];[1,0]]*K10hat;
H0001hat = -AINV(F.J1*K01hat);
gamma3 = -(p1*(F.B(H0001hat,q0)+F.A1(q0,K01hat)))/(2*a);
delta1 = 1/(p1*(F.B(H0001hat,q1)+F.A1(q1,K01hat)) ...
								+ p0*(F.B(H0001hat,q0)+F.A1(q0,K01hat))+gamma3*b);
gamma4 = (p1*H1100-p1*(F.B(H0010hat,q0)+F.A1(q0,K10hat)))/(2*a);
delta2 = -p1*(F.B(H0010hat,q1)+F.A1(q1,K10hat))-gamma4*b+p1*H0200 ...
				 -p0*(F.B(H0010hat,q0)+F.A1(q0,K10hat)-H1100);
K01 = delta1*K01hat;
H0001 = delta1*(H0001hat+gamma3*q0);
K10 = K10hat+delta2*K01;
H0010 = H0010hat+delta2*H0001+gamma4*q0;
%% H1010, H0100
H1010 = -AINV(F.B(H0010,q0)+F.A1(q0,K10)-H1100);
H0110 = -AINV(F.B(H0010,q1)+F.A1(q1,K10)-H0200-H1010);
%% H1001, H0101, H2001, H1101
H1001hat = -AINV(F.B(H0001,q0)+F.A1(q0,K01));
H0101hat = -AINV(F.A1(q1,K01)+F.B(H0001,q1)-H1001hat-q1);
zeta2 = p1*(-b*H0101hat-H1100+F.A1(H1100,K01) ...
				+F.B(H0001,H1100)+F.B(H0101hat,q0) ...
				+F.B(H1001hat,q1)+F.B1(q0,q1,K01)+F.C(H0001,q0,q1)) ...
				+p0*(-2*a*H0101hat+F.A1(H2000,K01)+F.B(H0001,H2000) ...
				+2*F.B(H1001hat,q0)+F.B1(q0,q0,K01)+F.C(H0001,q0,q0));
gamma5 = -zeta2/b;
H1001 = H1001hat + gamma5*q0;
H0101 = H0101hat + gamma5*q1;
a1 = 1/2*p1*(-2*a*H0101+F.A1(H2000,K01)+F.B(H0001,H2000)+2*F.B(H1001,q0) ...
				 +F.B1(q0,q0,K01)+F.C(H0001,q0,q0));
H2001 = -AINV(-2*a*H0101+F.A1(H2000,K01)+F.B(H0001,H2000)+2*F.B(H1001,q0) ...
				 -2*a1*q1+F.B1(q0,q0,K01)+F.C(H0001,q0,q0));
H1101 = -AINV(-b*H0101-H1100-H2001+F.A1(H1100,K01) ...
				+F.B(H0001,H1100)+F.B(H0101,q0)+ F.B(H1001,q1) ...
				+F.B1(q0,q1,K01)+F.C(H0001,q0,q1));
%% K11, H0011
K11 = -p1*(F.A1(H0001,K10)+F.A1(H0010,K01)+F.B(H0001,H0010)+F.J2(K01,K10) ...
				-H0101)*K10;
H0011 = -AINV(F.J1*K11+F.A1(H0001,K10)+F.A1(H0010,K01)+F.B(H0001,H0010) ...
				+F.J2(K01,K10)-H0101);
%% H02, H0002, H1002, H0102
K02hat = -(p1*(2*F.A1(H0001,K01) + F.B(H0001,H0001) + F.J2(K01,K01)))*K10;
H0002hat = -AINV( F.J1*K02hat + 2*F.A1(H0001,K01) + F.B(H0001,H0001) ...
								+ F.J2(K01,K01));
gamma6 = -(1/(2*a))*p1*(2*F.A1(H1001,K01)+F.A1(q0,K02hat)+F.A2(q0,K01,K01) ...
				+F.B(q0,H0002hat)+2*F.B(H0001,H1001)+2*F.B1(q0,H0001,K01) ...
				+F.C(q0,H0001,H0001));
delta3 = -p1*(2*F.A1(H0101,K01)+F.A1(q1,K02hat)+F.A2(q1,K01,K01) ...
				+F.B(q1,H0002hat)+2*F.B(H0001,H0101)+2*F.B1(q1,H0001,K01) ...
				+F.C(q1,H0001,H0001)-2*H0101) ...
				-p0*(2*F.A1(H1001,K01)+F.A1(q0,K02hat)+F.A2(q0,K01,K01) ...
				+F.B(q0,H0002hat)+2*F.B(H0001,H1001)+2*F.B1(q0,H0001,K01) ...
				+F.C(q0,H0001,H0001))-gamma6*b;
K02 = K02hat + delta3*K01;
H0002 = H0002hat + delta3*H0001 + gamma6*q0;
H1002 = -AINV(2*F.A1(H1001,K01)+F.A1(q0,K02)+F.A2(q0,K01,K01)+F.B(q0,H0002) ...
				+2*F.B(H0001,H1001)+2*F.B1(q0,H0001,K01)+F.C(q0,H0001,H0001));
H0102 = -AINV(2*F.A1(H0101,K01)+F.A1(q1,K02)+F.A2(q1,K01,K01)+F.B(q1,H0002) ...
				+2*F.B(H0001,H0101)+2*F.B1(q1,H0001,K01)+F.C(q1,H0001,H0001) ...
				-2*H0101-H1002);
%% K03, H0003
K03 = -p1*(F.A1(H0001,K02)+F.A1(H0002,K01)+2*(F.A1(H0001,K02) ...
				+F.A1(H0002,K01))+3*F.B(H0001,H0002)+3*F.J2(K01,K02) ...
				+3*F.A2(H0001,K01,K01)+3*F.B1(H0001,H0001,K01)+F.C(H0001,H0001,H0001) ...
				+F.J3(K01,K01,K01))*K10;
H0003 = -AINV(F.J1*K03+F.A1(H0001,K02)+F.A1(H0002,K01)+2*(F.A1(H0001,K02) ...
				+F.A1(H0002,K01))+3*F.B(H0001,H0002)+3*F.J2(K01,K02) ...
				+3*F.A2(H0001,K01,K01)+3*F.B1(H0001,H0001,K01)+F.C(H0001,H0001,H0001) ...
				+F.J3(K01,K01,K01));

%% calculate residuals
accuracy_cmf_coeff = max(abs( ...
   [p1*(-2*a*q1 + F.A*H2000 + F.B(q0,q0))
    p1*(-H2000-b*q1+F.A*H1100+F.B(q0,q1))
    p1*(-2*H1100 + F.A*H0200 + F.B(q1,q1))
    p1*(-6*a*H1100-6*d*q1+F.A*H3000+3*F.B(H2000,q0)+F.C(q0,q0,q0))
    p1*(-2*a*H0200-2*b*H1100-H3000+F.A*H2100+2*F.B(H1100,q0)+F.B(H2000,q1)+F.C(q0,q0,q1))
    p1*(-q1+F.A*H0010+F.J1*K10)
    p1*(F.A*H0001+F.J1*K01)
    p1*(F.A*H1001+F.A1(q0,K01)+F.B(H0001,q0))
    p1*(-H1001-q1+F.A*H0101+F.A1(q1,K01)+F.B(H0001,q1))
    p1*(-H1100+F.A*H1010+F.A1(q0,K10)+F.B(H0010,q0))
    p1*(-H0200-H1010+F.A*H0110+F.A1(q1,K10)+F.B(H0010,q1))
    p1*(-2*a*H0101-2*a1*q1+F.A*H2001+F.A1(H2000,K01)+F.B(H0001,H2000)+2*F.B(H1001,q0)+F.B1(q0,q0,K01)+F.C(H0001,q0,q0))
    p1*(-b*H0101-H1100-H2001+F.A*H1101+F.A1(H1100,K01)+1/2*(2*F.B(H0001,H1100)+2*F.B(H0101,q0)+2*F.B(H1001,q1))+F.B1(q0,q1,K01)+F.C(H0001,q0,q1))
    p1*(F.A*H0002+F.J1*K02+2*F.A1(H0001,K01)+F.B(H0001,H0001)+F.J2(K01,K01))
    p1*(-H0101+F.A*H0011+F.J1*K11+F.A1(H0001,K10)+F.A1(H0010,K01)+F.B(H0001,H0010)+F.J2(K01,K10))
    p1*(F.A*H1002+2*F.A1(H1001,K01)+F.A1(q0,K02)+1/2*(2*F.B(H0001,H1001)+2*(F.B(H0001,H1001)+F.B(H0002,q0)))+F.A2(q0,K01,K01)+2*F.B1(H0001,q0,K01)+F.C(H0001,H0001,q0))
    p1*(-2*H0101-H1002+F.A*H0102+2*F.A1(H0101,K01)+F.A1(q1,K02)+1/2*(2*F.B(H0001,H0101)+2*(F.B(H0001,H0101)+F.B(H0002,q1)))+F.A2(q1,K01,K01)+2*F.B1(H0001,q1,K01)+F.C(H0001,H0001,q1))
    p1*(F.A*H0003+F.J1*K03+F.A1(H0001,K02)+F.A1(H0002,K01)+2*((F.A1(H0001,K02)+F.A1(H0002,K01)))+3*F.B(H0001,H0002)+3*F.J2(K01,K02)+3*F.A2(H0001,K01,K01)+ 3*F.B1(H0001,H0001,K01)+F.C(H0001,H0001,H0001)+F.J3(K01,K01,K01))]));
if options.messages
    fprintf('Center manifold coefficients'' accuracy: %d\n', accuracy_cmf_coeff);
end

%% return normal form coefficients
bt.nmfm.a = a;
bt.nmfm.b = b;
bt.nmfm.d = d;
bt.nmfm.e = 0;
bt.nmfm.a1 = a1;
bt.nmfm.b1 = 0;
bt.nmfm.K10 = K10;
bt.nmfm.K01 = K01;
bt.nmfm.K02 = K02;
bt.nmfm.K11 = K11;
bt.nmfm.K03 = K03;
bt.nmfm.q0 = q0;
bt.nmfm.q1 = q1;
bt.nmfm.H0010 = H0010;
bt.nmfm.H0001 = H0001;
bt.nmfm.H2000 = H2000;
bt.nmfm.H1100 = H1100;
bt.nmfm.H0200 = H0200;
bt.nmfm.H1010 = H1010;
bt.nmfm.H1001 = H1001;
bt.nmfm.H0110 = H0110;
bt.nmfm.H0101 = H0101;
bt.nmfm.H0002 = H0002;
bt.nmfm.H0011 = H0011;
bt.nmfm.H3000 = H3000;
bt.nmfm.H2100 = H2100;
bt.nmfm.H1101 = H1101;
bt.nmfm.H2001 = H2001;
bt.nmfm.H0003 = H0003;
bt.nmfm.H1002 = H1002;
bt.nmfm.H0102 = H0102;
