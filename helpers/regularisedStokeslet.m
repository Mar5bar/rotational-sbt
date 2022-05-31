function res = regularisedStokeslet(x,y,chi)
%% Returns the tensor corresponding to a regularised Stokeslet at y evaluated
% at x for local regularisation parameter chi.
	r = x - y;
	norm_r_squared = r'*r;
	res = ((norm_r_squared + 2*chi) * eye(3) + r*r') / (norm_r_squared + chi)^(3/2);
end
