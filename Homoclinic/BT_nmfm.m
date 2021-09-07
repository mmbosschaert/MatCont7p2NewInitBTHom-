function bt = BT_nmfm(odefile, bt, ap)
% TODO: change bord_inv to AINV

%% Multi-linear forms
F = nmfm_deriv_define(odefile, bt, ap);

%% 4. Compute null vectors q0, q1, p0 and p1
[q0,p1,q1,p0]=nmfm_nullvectors(F.A);
% p1=-p1;
% p0=-p0;
% q0=-q0;
% q1=-q1;
Bord = [F.A p1; q0' 0]; % Setup bordered system (should be moved).
AINV = @(x) bord_inv(Bord, x);
%% 6. BT NORMAL FORM COEFFICIENTS (a,b)
a = 1/2*p1'*F.B(q0,q0);
b = p1'*F.B(q0,q1)+p0'*F.B(q0,q0);
%% 7.compute H2000, H1100, H0200, H00, K1, K22, H1001, H0002 and H0101
H2000 = bord_inv(Bord,2*a*q1-F.B(q0,q0));
gamma = 1/2*(-2*p0'*H2000 + 2*p0'*F.B(q0,q1) + p1'*F.B(q1,q1));
H2000 = H2000+gamma*q0;
H1100 = bord_inv(Bord,b*q1-F.B(q0,q1)+H2000);
H0200 = bord_inv(Bord,2*H1100-F.B(q1,q1));
%% Solve the 'Big' system
% Id=eye(homds.nphase);
% BigSystem=[
%     A                                J1
%     p1'*F.B(q0,Id)                   p1'*F.A1(q0,eye(2))
%     p0'*F.B(q0,Id)+p1'*F.B(q1,Id)    p0'*F.A1(q0,eye(2))+p1'*F.A1(q1,eye(2))
% ];
% RHS=[ q1                            zeros(size(A,1),1);
%       0.5*p1'*F.B(q1,q1)            0;
%       -p0'*F.B(q1,q1)+3*p0'*H1100   1
% ];
% HK=BigSystem\RHS;
% H0010=HK(1:end-2,1);
% H0001=HK(1:end-2,2);
% K10=HK(end-1:end,1);
% K01=HK(end-1:end,2);
%% solve without 'Big' system
% gamma=p1'*J1;
% s1 = (1/(gamma*gamma')*gamma)';
% s2 = [-gamma(2); gamma(1)];
% K10 = s1;
% K01 = s2;
% H0010 = bord_inv(Bord, q1 - J1*s1);
% H0001 = -bord_inv(Bord, J1*s2);
% xi2 = -(p1'*(F.A1(q0, K01) + F.B(q0, H0001)))/(2*a);
% H0001 = H0001 + xi2*q0;
% delta2 = 1/(p0'*(F.A1(q0, K01) + F.B(q0, H0001)) + p1'*(F.A1(q1, K01) + F.B(q1, H0001)));
% H0001 = delta2*H0001;
% K01 = delta2*K01;
% xi1 = (p0'*H2000 - p0'*F.B(q0, q1) - p1'*F.B(q0, H0010) - p1'*F.A1(q0, K10))/(2*a);
% H0010 = H0010 + xi1*q0;
% delta1 = 3*p0'*H1100 - p0'*F.B(q1, q1) - p0'*F.B(q0, H0010) - p0'*F.A1(q0, K10) - p1'*F.B(q1, H0010) - p1'*F.A1(q1, K10);
% H0010 = H0010 + delta1*H0001;
% K10 = K10 + delta1*K01;

eta = (p1'*F.J1)';
K10hat = 1/(eta'*eta)*eta;
H0010hat = AINV(q1-F.J1*K10hat);
K01hat = [[0,-1];[1,0]]*eta;
H0001hat = -AINV(F.J1*K01hat);
gamma3 = -(p1'*(F.B(H0001hat,q0)+F.A1(q0,K01hat)))/(2*a);
delta1 = 1/(p1'*(F.B(H0001hat,q1)+F.A1(q1,K01hat)) ...
								+ p0'*(F.B(H0001hat,q0)+F.A1(q0,K01hat))+gamma3*b);
gamma4 = (p1'*H1100-p1'*(F.B(H0010hat,q0)+F.A1(q0,K10hat)))/(2*a);
delta2 = -p1'*(F.B(H0010hat,q1)+F.A1(q1,K10hat))-gamma4*b+p1'*H0200 ...
				-p0'*(F.B(H0010hat,q0)+F.A1(q0,K10hat)-H1100);
K01 = delta1*K01hat;
H0001 = delta1*(H0001hat+gamma3*q0);
K10 = K10hat+delta2*K01;
H0010 = H0010hat+delta2*H0001+gamma4*q0;

%% H1001, H0101, H0002, K02, H1002, H0102
H1001=-bord_inv(Bord,F.B(q0,H0001)+F.A1(q0,K01));
H0101=-bord_inv(Bord,F.B(q1,H0001)+F.A1(q1,K01)-H1001-q1);
%
K02=-p1'*(F.B(H0001,H0001)+2*F.A1(H0001,K01)+F.J2(K01,K01))*K10;
H0002=-bord_inv(Bord, F.J1*K02 + F.B(H0001,H0001) + 2*F.A1(H0001,K01) + F.J2(K01,K01));

% save K02 for comparing predictors
K02_old = K02;

xi6 = -1/(2*a)*p1'*(F.A1(q0,K02) + 2*F.A1(H1001,K01) + ...
        F.B(q0,H0002) + 2*F.B(H0001,H1001) + ...
        F.A2(q0,K01,K01) + 2*F.B1(q0,H0001,K01) + F.C(q0,H0001,H0001));
H0002 = H0002 + xi6*q0;
H1002 = -bord_inv(Bord, F.A1(q0,K02) + 2*F.A1(H1001,K01) + ...
        F.B(q0,H0002) + 2*F.B(H0001,H1001) + ...
        F.A2(q0,K01,K01) + 2*F.B1(q0,H0001,K01) + F.C(q0,H0001,H0001));
delta3 =  -p1'*(F.A1(q1,K02) + 2*F.A1(H0101,K01) + ...
        F.B(q1,H0002) + 2*F.B(H0001,H0101) + ...
        F.A2(q1,K01,K01) + 2*F.B1(q1,H0001,K01) + F.C(q1,H0001,H0001) - ...
        2*H0101-H1002);
H0002 = H0002 + delta3*H0001;
K02   = K02 + delta3*K01;
H1002 = H1002 + delta3*H1001;
H0102 = -bord_inv(Bord, F.A1(q0,K02) + 2*F.A1(H1001,K01) + ...
        F.B(q0,H0002) + 2*F.B(H0001,H1001) + ...
        F.A2(q0,K01,K01) + 2*F.B1(q0,H0001,K01) + F.C(q0,H0001,H0001));
%% K11, H0011, H0003, H1010, H0110
K11 = p1'*(H0101 - (F.A1(H0001,K10) + F.A1(H0010,K01) + F.B(H0010,H0001) + F.J2(K10,K01)))*K10;
H0011 = -bord_inv(Bord, -H0101 + F.J1*K11 + F.A1(H0001,K10) + F.A1(H0010,K01) + ...
        F.B(H0001,H0010) + F.J2(K01,K10));
K03 = -p1'*(F.A1(H0002,K01) + F.A1(H0001,K02) + ...
     2*(F.A1(H0002,K01) + F.A1(H0001,K02)) + 3*F.B(H0001,H0002) + ...
     3*F.J2(K01,K02) + 3*F.A2(H0001,K01,K01) + 3*F.B1(H0001,H0001,K01) + ...
     F.C(H0001,H0001,H0001) + F.J3(K01,K01,K01))*K10;
H0003 = -bord_inv(Bord, F.J1*K03 + F.A1(H0002,K01) + F.A1(H0001,K02) + ...
    2*(F.A1(H0002,K01) + F.A1(H0001,K02)) + 3*F.B(H0001,H0002) + 3*F.J2(K01,K02) + ...
    3*F.A2(H0001,K01,K01) + 3*F.B1(H0001,H0001,K01) + F.C(H0001,H0001,H0001) + ...
    F.J3(K01,K01,K01));
H1010 = bord_inv(Bord, H1100 - F.B(H0010,q0) - F.A1(q0,K10));
H0110 = bord_inv(Bord, H0200 + H1010 - F.B(H0010,q1) - F.A1(q1,K10));
%% Compute d, e, a1 AND b1 (The BT normal form coefficients)
d=p1'*(1/6*F.C(q0,q0,q0)+0.5*F.B(q0,H2000)-a*H1100);
H3000=bord_inv(Bord,-6*(1/6*F.C(q0,q0,q0)+0.5*F.B(q0,H2000)-a*H1100-d*q1));
e=p1'*(0.5*F.C(q0,q0,q1)+F.B(q0,H1100)+0.5*F.B(q1,H2000)-b*H1100-a*H0200-0.5*H3000);
a1=p1'*(0.5*F.C(q0,q0,H0001)+0.5*F.B1(q0,q0,K01)+F.B(q0,H1001)+0.5*F.B(H0001,H2000)+0.5*F.A1(H2000,K01)-a*H0101);
H2001=bord_inv(Bord,-2*(0.5*F.C(q0,q0,H0001)+0.5*F.B1(q0,q0,K01)+F.B(q0,H1001)+0.5*F.B(H0001,H2000)+0.5*F.A1(H2000,K01)-a*H0101-a1*q1));
b1=p1'*(F.C(q0,q1,H0001)+F.B1(q0,q1,K01)+F.B(q1,H1001)+F.B(H0001,H1100)+F.B(q0,H0101)+F.A1(H1100,K01)-b*H0101-H1100-H2001);

H2100 = bord_inv(Bord, 2*a*H0200 + 2*b*H1100 + H3000 + 2*e*q1 - ...
    2*F.B(H1100,q0) - F.B(H2000,q1) - F.C(q0,q0,q1));
H1101 = -bord_inv(Bord, -b1*q1 - b*H0101 - H1100 - H2001 + F.B(q0,H0101) + ...
    F.A1(H1100,K01) + F.B(q1,H1001) + F.B(H0001,H1100) + F.B1(q0,q1,K01) + ...
    F.C(q0,q1,H0001));
%% return normal form coefficients
bt.nmfm.a = a;
bt.nmfm.b = b;
bt.nmfm.d = d;
bt.nmfm.e = e;
bt.nmfm.a1 = a1;
bt.nmfm.b1 = b1;
bt.nmfm.K10 = K10;
bt.nmfm.K01 = K01;
bt.nmfm.K02 = K02;
bt.nmfm.K02_old = K02_old;
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
