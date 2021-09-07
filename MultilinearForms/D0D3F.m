function vec3 = D0D3F(odefile,p1,p2,p3,x0,p,increment,ap)

% This file computes the multilinear function J3(p1,p2,p3) where
% J3 = D0D3(F(x0)), the 3rd derivative of the vectorfield wrt to the parameters

nphase=size(x0,1);
if (p1==p2)
    if (p1==p3)
        vec3 = J3vvv(odefile,p1,x0,p,increment,ap);
    else
        part1 = J3vvv(odefile,p1+p3,x0,p,increment,ap);
        part2 = J3vvv(odefile,p1-p3,x0,p,increment,ap);
        part3 = J3vvv(odefile,p3,x0,p,increment,ap);
        vec3 = (part1 - part2 - 2.0*part3)/6.0;
    end
else
    part1 = J3vvv(odefile,p1+p2+p3,x0,p,increment,ap);
    part2 = J3vvv(odefile,p1+p2-p3,x0,p,increment,ap);
    part3 = J3vvv(odefile,p1-p2+p3,x0,p,increment,ap);
    part4 = J3vvv(odefile,p1-p2-p3,x0,p,increment,ap);
    vec3 = (part1 - part2 - part3 + part4)/24.0;
end

function tempvec = J3vvv(odefile,vp,x0,p,increment,ap)
    par0 = p*0; 
    dp = @(vals) matcont_array_insert(par0, ap, vals);

    if ~norm(vp) == 0
        vs=vp/norm(vp);
    else
        tempvec = 0;
        return
    end
    p1 = p + 3.0*increment*dp(vs);
    p2 = p +     increment*dp(vs);
    p3 = p -     increment*dp(vs);
    p4 = p - 3.0*increment*dp(vs);
    p1 = num2cell(p1);
    p2 = num2cell(p2);
    p3 = num2cell(p3);
    p4 = num2cell(p4);
    p1 = feval(odefile, 0, x0, p1{:});
    p2 = feval(odefile, 0, x0, p2{:});
    p3 = feval(odefile, 0, x0, p3{:});
    p4 = feval(odefile, 0, x0, p4{:});
    tempvec = (p1 - 3.0*p2 + 3.0*p3 - p4)/8.0*(norm(vp)/increment)^3;
