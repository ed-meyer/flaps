
function [V,D] = gyroeig(m,g,k)

% create I
[n,nc] = size(m);
z = zeros(n,n);
I = [m z;z k];
% create K = [ g'*mi*g + k  g'*mi*k ]
%            [ k*mi*g         k*mi*k  ]
mi = inv(m)
gmigk = g'*mi*g + k
gmik = g'*mi*k
kmik = k*mi*k
K = [ gmigk gmik;gmik' kmik]
[V,D] = eig(K,I)

