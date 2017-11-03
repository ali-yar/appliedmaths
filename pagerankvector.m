function [v, iteration] = eigenvector_func (P, w0)
	% if matrix P not given, create one
	if ~exist('P')
		P = [0 0.5 0 0 0.2;
		0.5 0 0 0 0.2;
		0 0.5 0 0.5 0.2;
		0.5 0 0 0 0.2;
		0 0 1 0.5 0.2];
	end
	[M,N] = size(P);
	% check if P is a square matrix
	if M != N
		error('The input matrix is not a square matrix.')
	end
	% check if each column of P sums to 1
	if ~isequal(sum(P),ones(1,N))
		error('The input matrix is not a Markov matrix.')
	end
	% if initial vector w0 not given, create one
	if ~exist('w0')
		w0 = (1/N) * [ones(N,1)];
	end
	r = .85; % porbability parameter
	s = 1 - r; % compelment of probability parameter
	u = ones(N,1); % vector of 1's
	rP = r * P; % first part of equation of w
	n = (s / N) * u; % second part of equation of w
	count = 0; % keep track of repetitions
	max = 100; % max number of iterations
	w = w0; % keep track of latest computed w
	while (count < max)
		w_new = rP * w + n; % compute new w
		% check if reached a stable state
		if isequal(w, w_new)
			break; % stop the loop if we found the final w
		end
		w = w_new; % update w for the next iteration
		count++;
	end
	v = w; % return eigenvector
	iteration = count; % return number of iterations to find eigenvector
endfunction