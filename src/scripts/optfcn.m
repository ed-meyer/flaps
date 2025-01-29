
% optimization function evaluation with fminsearch
function  f = optfcn(x)
	% f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
	% send x to flaps to compute new RFA, send back error
	enviar = memmapfile('betas', 'Writable',true,'Format', 'double');
	recibir = memmapfile('function', 'Writable',true,'Format', 'double');
	% m.Data(1) = 0.0;	% message received
	len = length(x);
	enviar.Data(2:len+1) = x;
	enviar.Data(1) = len;		% message sent
	fprintf(1,'sent x len %d\n',len)
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
	% the reply will be only 1 double
	len = recibir.Data(1)
	f = recibir.Data(2)
	fprintf(1,'got reply length %d, f %9.5f',len,f)
	recibir.Data(1) = 0.0;	% message received

