function R = genRRotletAnsatzRTT(in)
% Form the RHS of the linear system for the rotlet ansatz with RTT, 
% which is simply Omega(s).
	R = in.Omega(:);
end