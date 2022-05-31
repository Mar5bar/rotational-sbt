function R = genRRotletAnsatzRTTI1Approx(in)
% Form the RHS of the linear system for the rotlet ansatz with RTT and I1
% approximated, which is simply Omega(s).
	R = in.Omega(:);
end