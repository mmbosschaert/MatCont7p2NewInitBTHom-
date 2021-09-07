function [q0,p1,q1,p0] = nmfm_nullvectors(A)
% calculuate the eigenvectors and adjoint eigenvectors
% at a codimension 2 Bogdanov-Takens point
[X,D] = eig(A);
nphase = length(A(1,:));
index1 = find(abs(diag(D)) < 1e-3); %If ok, index1 is 1x2 array otherwise
vext = real(X(:,index1(1)));
[X,D] = eig(A');
index1 = find(abs(diag(D)) < 1e-3);
wext = real(X(:,index1(1)));
% correct the eigenvectors using bordered systems
Bord = [ A wext; vext' 0];
bunit=[zeros(nphase,1);1];
q0=Bord\bunit;
q0=q0(1:nphase);          % A q0 = 0 , <vext,q0> = 1
p1=Bord'\bunit;
p1=p1(1:nphase);          % A'p1 = 0 , <wext,p1> = 1
Bord = [ A p1; q0' 0];
q1 = Bord\[q0; 0];
q1 = q1(1:nphase);		% A q1 = q0, <q0, q1>  = 0
p0 = Bord'\[p1; 0];
p0 = p0(1:nphase);		% A'p0 = p1, <p0, p1>  = 0
% normalize such that <p0,q0>=<p1,q1>=1 & <p0,q1>=<p1,q0>=0
mu = sqrt(abs(q0'*q0));
q0 = (1/mu)*q0;
q1 = (1/mu)*q1;
q1 = q1 - (q0'*q1)*q0;
nu = q0'*p0;
p1 = (1/nu)*p1;
p0 = p0 - (p0'*q1)*p1;
p0 = (1/nu)*p0;
end
