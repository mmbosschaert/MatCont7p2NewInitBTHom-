function F = SIRmodel_multilinearforms
F.A = @A;
F.J1 = @J1;
F.B = @B;
F.A1 = @A1;
F.J2 = @J2;
F.C = @C;
F.J3 = @J3;
F.B1 = @B1;
F.A2 = @A2;

function A = A(in1,in2)
%A
%    A = A(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:27

I = in1(:,2);
R = in1(:,3);
S = in1(:,1);
b = in2(:,2);
mu1 = in2(:,1);
t2 = I+b;
t3 = I+R+S;
t4 = mu1-1.0e+1;
t5 = 1.0./t2;
t7 = 1.0./t3;
t6 = t5.^2;
t8 = t7.^2;
t9 = I.*t7.*(2.3e+1./2.0);
t10 = S.*t7.*(2.3e+1./2.0);
t11 = b.*t4.*t5;
t12 = I.*b.*t4.*t6;
t13 = I.*S.*t8.*(2.3e+1./2.0);
t14 = -t13;
A = reshape([-t9+t13-1.0./1.0e+1,t9+t14,0.0,-t10+t13,t10-t11+t12+t14-1.11e+2./1.0e+1,t11-t12+1.0e+1,t13,t14,-1.0./1.0e+1],[3,3]);

function J1 = J1(in1,in2)
%J1
%    J1 = J1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
b = in2(:,2);
mu1 = in2(:,1);
t2 = I+b;
t3 = mu1-1.0e+1;
t4 = 1.0./t2;
t5 = t4.^2;
t6 = I.*b.*t4;
t7 = t3.*t4;
t8 = b.*t3.*t5;
t9 = -t8;
t10 = t7+t9;
t11 = I.*t10;
J1 = reshape([0.0,-t6,t6,0.0,-t11,t11],[3,2]);

function B = B(in1,in2,in3,in4)
%B
%    B = B(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
R = in1(:,3);
S = in1(:,1);
b = in2(:,2);
mu1 = in2(:,1);
v11 = in3(1,:);
v12 = in3(2,:);
v13 = in3(3,:);
v21 = in4(1,:);
v22 = in4(2,:);
v23 = in4(3,:);
t2 = I+b;
t3 = I+R+S;
t4 = mu1-1.0e+1;
t5 = 1.0./t2.^2;
t6 = 1.0./t2.^3;
t7 = 1.0./t3;
t8 = t7.^2;
t9 = t7.^3;
t12 = t7.*(2.3e+1./2.0);
t22 = b.*t4.*t5.*2.0;
t23 = I.*b.*t4.*t6.*2.0;
t10 = I.*t8.*2.3e+1;
t11 = S.*t8.*2.3e+1;
t13 = I.*S.*t9.*2.3e+1;
t14 = -t12;
t16 = I.*t8.*(2.3e+1./2.0);
t17 = S.*t8.*(2.3e+1./2.0);
t19 = I.*S.*t9.*v13.*-2.3e+1;
t15 = t13.*v13;
t18 = -t13;
t20 = t16.*v13;
t21 = t17.*v13;
t28 = -v11.*(t13-t16);
t29 = -v12.*(t13-t17);
t31 = -v11.*(t12+t13-t16-t17);
t32 = -v12.*(t12+t13-t16-t17);
t34 = -v23.*(v11.*(t13-t16)+v12.*(t13-t17)+I.*S.*t9.*v13.*2.3e+1);
t24 = t10+t18;
t26 = t16+t18;
t27 = t17+t18;
t33 = t19+t28+t29;
t25 = t24.*v11;
t30 = t14+t17+t26;
t35 = t19+t20+t25+t32;
t36 = v21.*(t19+t25+t32+I.*t8.*v13.*(2.3e+1./2.0));
B = [t34+t36+v22.*(t19+t31+v12.*(t11+t18)+S.*t8.*v13.*(2.3e+1./2.0));-t36+v23.*(v11.*(t13-t16)+v12.*(t13-t17)+I.*S.*t9.*v13.*2.3e+1)+v22.*(t15-v12.*(t11+t18-t22+t23)+v11.*(t12+t13-t16-t17)-S.*t8.*v13.*(2.3e+1./2.0));-v12.*v22.*(t22-t23)];

function J2 = J2(in1,in2,in3,in4)
%J2
%    J2 = J2(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
b = in2(:,2);
mu1 = in2(:,1);
p11 = in3(1,:);
p12 = in3(2,:);
p21 = in4(1,:);
p22 = in4(2,:);
t2 = I+b;
t3 = mu1-1.0e+1;
t4 = 1.0./t2;
t5 = t4.^2;
t6 = t4.^3;
t8 = I.*p11.*t4;
t9 = -t4;
t7 = b.*t5;
t11 = -t8;
t12 = t3.*t5.*2.0;
t13 = b.*t3.*t6.*2.0;
t10 = I.*p11.*t7;
t14 = -t13;
t15 = t7+t9;
t16 = -I.*p12.*p21.*(t4-t7);
t17 = t12+t14;
t18 = I.*p12.*t17;
t19 = t10+t11+t18;
t20 = p22.*t19;
J2 = [0.0;t16+t20;-t20+I.*p12.*p21.*(t4-t7)];

function A1 = A1(in1,in2,in3,in4)
%A1
%    A1 = A1(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
b = in2(:,2);
mu1 = in2(:,1);
p11 = in4(1,:);
p12 = in4(2,:);
v12 = in3(2,:);
t2 = I+b;
t3 = mu1-1.0e+1;
t4 = 1.0./t2;
t5 = t4.^2;
t6 = t4.^3;
t7 = b.*t4;
t10 = t3.*t4;
t8 = I.*b.*t5;
t11 = I.*t3.*t5;
t12 = b.*t3.*t5;
t13 = I.*b.*t3.*t6.*2.0;
t9 = -t8;
t14 = -t11;
t15 = -t12;
t16 = t7+t9;
t18 = t10+t13+t14+t15;
t17 = p11.*t16.*v12;
t19 = p12.*t18.*v12;
A1 = [0.0;-t17-t19;t17+t19];

function C = C(in1,in2,in3,in4,in5)
%C
%    C = C(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
R = in1(:,3);
S = in1(:,1);
b = in2(:,2);
mu1 = in2(:,1);
v11 = in3(1,:);
v12 = in3(2,:);
v13 = in3(3,:);
v21 = in4(1,:);
v22 = in4(2,:);
v23 = in4(3,:);
v31 = in5(1,:);
v32 = in5(2,:);
v33 = in5(3,:);
t2 = I+b;
t3 = I+R+S;
t4 = mu1-1.0e+1;
t5 = 1.0./t2.^3;
t6 = 1.0./t2.^4;
t7 = 1.0./t3.^2;
t8 = 1.0./t3.^3;
t9 = t7.^2;
t10 = t7.*2.3e+1;
t12 = I.*t8.*2.3e+1;
t13 = I.*t8.*4.6e+1;
t14 = S.*t8.*2.3e+1;
t15 = S.*t8.*4.6e+1;
t16 = t7.*(2.3e+1./2.0);
t17 = I.*t8.*6.9e+1;
t18 = S.*t8.*6.9e+1;
t27 = t7.*v13.*(-2.3e+1./2.0);
t29 = b.*t4.*t5.*6.0;
t31 = I.*b.*t4.*t6.*6.0;
t11 = -t10;
t19 = t12.*v13;
t20 = t13.*v13;
t21 = t14.*v13;
t22 = t15.*v13;
t23 = -t16;
t24 = t16.*v13;
t25 = I.*S.*t9.*6.9e+1;
t30 = I.*S.*t9.*v13.*-6.9e+1;
t26 = t25.*v13;
t28 = -t25;
t44 = -v11.*(t10-t12-t15+t25);
t45 = -v11.*(t10-t13-t14+t25);
t46 = -v12.*(t10-t12-t15+t25);
t47 = -v12.*(t10-t13-t14+t25);
t32 = t12+t28;
t33 = t13+t28;
t34 = t14+t28;
t35 = t15+t28;
t36 = t17+t28;
t61 = t19+t21+t27+t30+t45+t46;
t62 = -v21.*(t24+t26+v11.*(t10-t13-t14+t25)+v12.*(t10-t12-t15+t25)-I.*t8.*v13.*2.3e+1-S.*t8.*v13.*2.3e+1);
t63 = -v22.*(t24+t26+v11.*(t10-t13-t14+t25)+v12.*(t10-t12-t15+t25)-I.*t8.*v13.*2.3e+1-S.*t8.*v13.*2.3e+1);
t37 = t32.*v11;
t38 = t33.*v11;
t39 = t34.*v12;
t40 = t35.*v12;
t41 = t36.*v11;
t42 = t11+t15+t32;
t43 = t11+t14+t33;
t48 = t14+t23+t32;
t49 = t48.*v11;
t50 = t48.*v12;
t51 = t30+t37+t39;
t53 = t20+t30+t41+t47;
t56 = v21.*(t30+t41+t47+I.*t8.*v13.*4.6e+1);
t52 = t51.*v23;
t54 = t19+t30+t38+t50;
t55 = t21+t30+t40+t49;
t57 = v22.*(t30+t40+t49+S.*t8.*v13.*2.3e+1);
t58 = v23.*(t30+t40+t49+S.*t8.*v13.*2.3e+1);
t59 = v21.*(t30+t38+t50+I.*t8.*v13.*2.3e+1);
t60 = v23.*(t30+t38+t50+I.*t8.*v13.*2.3e+1);
t64 = t52+t57+t59;
t66 = t56+t60+t63;
t65 = t64.*v33;
t67 = t66.*v31;
C = [-t65-t67-v32.*(t58+t62+v22.*(t30+t44+v12.*(t18+t28)+S.*t8.*v13.*4.6e+1));t65+t67+v32.*(t58+t62+v22.*(t44+v12.*(t18+t28-t29+t31)+S.*t8.*v13.*4.6e+1-I.*S.*t9.*v13.*6.9e+1));v12.*v22.*v32.*(t29-t31)];

function J3 = J3(in1,in2,in3,in4,in5)
%J3
%    J3 = J3(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
b = in2(:,2);
mu1 = in2(:,1);
p11 = in3(1,:);
p12 = in3(2,:);
p21 = in4(1,:);
p22 = in4(2,:);
p31 = in5(1,:);
p32 = in5(2,:);
t2 = I+b;
t3 = mu1-1.0e+1;
t4 = 1.0./t2.^2;
t5 = 1.0./t2.^3;
t6 = t4.^2;
t7 = t4.*2.0;
t8 = b.*t5.*2.0;
t11 = I.*p11.*t4.*-2.0;
t13 = t3.*t5.*6.0;
t9 = -t7;
t10 = I.*p11.*t7;
t12 = I.*p11.*t8;
t14 = b.*t3.*t6.*6.0;
t17 = -I.*p12.*p21.*(t7-t8);
t18 = -I.*p12.*p22.*p31.*(t7-t8);
t15 = -t14;
t16 = t8+t9;
t19 = t13+t15;
t20 = I.*p12.*t19;
t21 = t11+t12+t20;
t22 = p22.*t21;
t23 = t17+t22;
t24 = p32.*t23;
J3 = [0.0;-t24+I.*p12.*p22.*p31.*(t7-t8);t18+t24];

function B1 = B1(in1,in2,in3,in4,in5)
%B1
%    B1 = B1(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
b = in2(:,2);
mu1 = in2(:,1);
p11 = in5(1,:);
p12 = in5(2,:);
v12 = in3(2,:);
v22 = in4(2,:);
t2 = I+b;
t3 = mu1-1.0e+1;
t4 = 1.0./t2.^2;
t5 = 1.0./t2.^3;
t6 = t4.^2;
t7 = b.*t4.*2.0;
t8 = I.*b.*t5.*2.0;
t10 = t3.*t4.*2.0;
t11 = I.*t3.*t5.*2.0;
t12 = b.*t3.*t5.*4.0;
t9 = -t8;
t13 = I.*b.*t3.*t6.*6.0;
t14 = -t11;
t15 = -t12;
t16 = t7+t9;
t18 = t10+t13+t14+t15;
t17 = p11.*t16.*v12.*v22;
t19 = p12.*t18.*v12.*v22;
B1 = [0.0;t17+t19;-t17-t19];

function A2 = A2(in1,in2,in3,in4,in5)
%A2
%    A2 = A2(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Sep-2021 23:45:28

I = in1(:,2);
b = in2(:,2);
mu1 = in2(:,1);
p11 = in4(1,:);
p12 = in4(2,:);
p21 = in5(1,:);
p22 = in5(2,:);
v12 = in3(2,:);
t2 = I+b;
t3 = mu1-1.0e+1;
t4 = 1.0./t2;
t5 = t4.^2;
t6 = t4.^3;
t10 = -t4;
t7 = t5.^2;
t8 = I.*t5;
t9 = b.*t5;
t11 = I.*b.*t6.*2.0;
t13 = t3.*t5.*2.0;
t14 = I.*t3.*t6.*4.0;
t15 = b.*t3.*t6.*2.0;
t12 = -t11;
t16 = I.*b.*t3.*t7.*6.0;
t17 = -t14;
t18 = -t15;
t20 = -p11.*v12.*(t4-t8-t9+t11);
t21 = -p12.*p21.*v12.*(t4-t8-t9+t11);
t19 = t8+t9+t10+t12;
t22 = t13+t16+t17+t18;
t23 = p12.*t22.*v12;
t24 = t20+t23;
t25 = p22.*(t23-p11.*v12.*(t4-t8-t9+t11));
A2 = [0.0;t21+t25;-t25+p12.*p21.*v12.*(t4-t8-t9+t11)];

