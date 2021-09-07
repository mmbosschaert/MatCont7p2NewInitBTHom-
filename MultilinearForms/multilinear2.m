function vec2 = multilinear2(odefile,tens2,q1,q2,x0,p,increment)
%----------------------------------------------------------
%This file computes the multilinear function B(q1,q2) where
%B = D^2(F(x0)), the second derivative of the vectorfield 
%----------------------------------------------------
nphase=size(x0,1);
if (~isempty(tens2))
    vec2=tensor2op(tens2,q1,q2,nphase);
else
    if (q1==q2)
        vec2 = Bvv(odefile,q1,x0,p,increment);
    else
        part1 = Bvv(odefile,q1+q2,x0,p,increment);
        part2 = Bvv(odefile,q1-q2,x0,p,increment);
        vec2 = (part1-part2)/4.0;
    end
end
%----------------------------------------------------
function tempvec = Bvv(odefile,vq,x0,p,increment)
    if norm(vq) ~= 0
        vs=vq/norm(vq);
    else
        tempvec = 0;
        return
    end
    f0 = x0;
    f1 = x0 + increment*(vs);
    f2 = x0 - increment*(vs);
    f0 = feval(odefile, 0, f0, p{:});
    f1 = feval(odefile, 0, f1, p{:});
    f2 = feval(odefile, 0, f2, p{:});
    tempvec = (f1+f2-2.0*f0)*(norm(vq)/increment)^2;
