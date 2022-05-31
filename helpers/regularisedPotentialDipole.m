function res = regularisedPotentialDipole(x,y,chi)
%% Returns the tensor corresponding to a regularised potential dipole at y
% evaluated at x for local regularisation parameter chi.
	r = x - y;
	norm_r_squared = r'*r;
	res = (-(norm_r_squared - 2*chi) * eye(3) + 3*r*r') / (norm_r_squared + chi)^(5/2);
end
