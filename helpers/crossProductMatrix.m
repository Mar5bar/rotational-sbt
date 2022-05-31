function mat = crossProductMatrix(v)
%% Compute the matrix representation of doing cross(a,v), returning the
% matrix operator acting on [a1;a2;a3].
    mat = [0,v(3),-v(2);-v(3),0,v(1);v(2),-v(1),0];
end