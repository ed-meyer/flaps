
% function to find a Hopf bifurcation point using Eqn 7.24 Seydel
function  g = hopfpt(t,z)
	% f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
	% send x to flaps to compute new RFA, send back error
	enviar = memmapfile('z', 'Writable',true,'Format', 'double');
	recibir = memmapfile('g', 'Writable',true,'Format', 'double');
	% m.Data(1) = 0.0;	% message received
	len = length(x);
	enviar.Data(2:len+1) = z;
	enviar.Data(1) = len;		% message sent
	fprintf(1,'sent z len %d\n',len)
	% wait for the client to read the data and set the first element to zero
	while enviar.Data(1) ~= 0.0
		pause(1);
	end
	'wait for reply'
	% now wait for a reply
	% wait until first element is non-zero
	while recibir.Data(1) == 0.0
		pause(1);
	end
	% the reply will be g: nt doubles
	len = recibir.Data(1)
	g = recibir.Data(2:len+1)
	fprintf(1,'got reply length %d, f %9.5f',len,g)
	recibir.Data(1) = 0.0;	% message received

