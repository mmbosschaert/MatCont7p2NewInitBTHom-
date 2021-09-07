function [x,x0,p,T,eps0,eps1,YS,YU] = bt_rearr(x1)

% Rearranges x1 into all of its components
global homds
x = x1(homds.coords);
x0 = x1(homds.ncoords+1:homds.ncoords+homds.nphase);

p = homds.P0;
p(homds.ActiveParams) = x1(homds.PeriodIdx+1:homds.PeriodIdx+2);
idx = homds.PeriodIdx+3;

if homds.extravec(1)
    T = x1(idx);
    idx = idx+1;
else
    T = homds.T;
end

if homds.extravec(2)
    eps0 = x1(idx);
    idx = idx+1;
else
    eps0 = homds.eps0;
end

if homds.extravec(3)
    eps1 = x1(idx);
    idx = idx+1;
else
    eps1 = homds.eps1;
end

    YU = reshape(x1(idx:idx+homds.npos*homds.nneg-1),homds.nneg,homds.npos);
    idx = idx + homds.npos*homds.nneg;
    YS = reshape(x1(idx:idx+homds.npos*homds.nneg-1),homds.npos,homds.nneg);