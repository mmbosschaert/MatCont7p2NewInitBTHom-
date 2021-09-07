function H = bord_inv(Bord,RHS)

H=Bord\[RHS; 0];
H=H(1:end-1);

end