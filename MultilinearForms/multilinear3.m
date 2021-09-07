function vec3 = multilinear3(odefile,tens3,q1,q2,q3,x0,p,increment)
%--------------------------------------------------------------
% This file computes the multilinear function C(q1,q2,q3) where
% C = D^3(F(x0)), the 3rd derivative of the vectorfield wrt to phase
% variables only. 
%--------------------------------------------------------------
nphase=size(x0,1);
if (~isempty(tens3))
    vec3=tensor3op(tens3,q1,q2,q3,nphase);
else
    if (q1==q2)
        if (q1==q3)
            vec3 = Cvvv(odefile,q1,x0,p,increment);
        else
            part1 = Cvvv(odefile,q1+q3,x0,p,increment);
            part2 = Cvvv(odefile,q1-q3,x0,p,increment);
            part3 = Cvvv(odefile,q3,x0,p,increment);
            vec3 = (part1 - part2 - 2.0*part3)/6.0;
        end
    else
        part1 = Cvvv(odefile,q1+q2+q3,x0,p,increment);
        part2 = Cvvv(odefile,q1+q2-q3,x0,p,increment);
        part3 = Cvvv(odefile,q1-q2+q3,x0,p,increment);
        part4 = Cvvv(odefile,q1-q2-q3,x0,p,increment);
        vec3 = (part1 - part2 - part3 + part4)/24.0;
    end
end
%----------------------------------------------------
function tempvec = Cvvv(odefile,vq,x0,p,increment)
    if norm(vq) ~= 0
        vs=vq/norm(vq);
    else
        tempvec = 0;
        return
    end
    f1 = x0 + 3.0*increment*vs;
    f2 = x0 +     increment*vs;
    f3 = x0 -     increment*vs;
    f4 = x0 - 3.0*increment*vs;
    f1 = feval(odefile, 0, f1, p{:});
    f2 = feval(odefile, 0, f2, p{:});
    f3 = feval(odefile, 0, f3, p{:});
    f4 = feval(odefile, 0, f4, p{:});
    tempvec = (f1 - 3.0*f2 + 3.0*f3 - f4)/8.0*(norm(vq)/increment)^3;

