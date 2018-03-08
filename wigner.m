function W = wigner(Ex,X)
%MYWIGNER: Calculates the Wigner distribution from a column vector
%
%	W  = mywigner(Ex)
%
%	W  = output Wigner distribution
%	Ex = Input electric field (MUST be a column vector
%
%	Notes:
%		W = Int(-inf..inf){E(x+y)E(x-y)exp[2ixy]}
%
%		E(x+y) & E(x-y) are calculated via a FFT (fast Fourier transform) using the
%		shift theorem. The integration is performed via a FFT. Thus it is important
%		for the data to satisfy the sampling theorem:
%		dy = 2*pi/X			X = span of all x-values	dy = y resolution
%		dx = 2*pi/Y			Y = span of all y-values	dx = x resolution
%		The data must be completely contained within the range x(0)..x(N-1) &
%		y(0)..y(N-1) (i.e. the function must fall to zero within this range).
%
%	v1.0
%
%	Currently waiting for update:
%		Remove the fft/ifft by performing this inside the last function calls
%		Allow an arbitrary output resolution
%		Allow an input vector for x (and possibly y).

if (size(Ex, 2)-1)
    error('E(x) must be a column vector');
end

N = length(Ex);														%   Get length of vector
nn = -(N-1)/2:(N-1)/2;                    							%   Generate linear vector
x = (2*nn*pi/(N*abs(X(2)-X(1))));
EX1 = ifft( (fft(Ex)*ones(1,N)).*exp(+1i*x*X/2 ));					%   +ve shift
EX2 = ifft( (fft(Ex)*ones(1,N)).*exp(-1i*x*X/2 ));					%   -ve shift

figure
hold on
plot(abs(Ex).^2)
plot(abs(EX1).^2)
plot(abs(EX2).^2)


W = real(fftshift(fft(fftshift(EX1.*conj(EX2), 2), [], 2), 2));		%   Wigner function


% function [tfr, t, f] = wigner(x, t, N)
% %WVCH Wigner-Ville distribution
% %	[TFR,T,F]=WVCH(X,T,N) computes the Wigner-Ville distribution
% %	of a discrete-time signal X, which is given in the form (1xN), 
% %   or the cross Wigner-Ville representation between two signals. 
% % 
% %	X     : signal if auto-WV, or [X1,X2] if cross-WV.
% %	T     : time instant(s)          (default : 1:length(X)).
% %	N     : number of frequency bins (default : length(X)).
% %	TFR   : time-frequency representation.
% %	F     : vector of normalized frequencies.
% 
% if (nargin == 0),
%  error('At least one parameter required');
% end;
% 
% [xrow,xcol] = size(x);
% if xrow ~= 1,
%     error('X must have one row, (1xN) form');
% end;
% 
% if (nargin == 1),
%  t=1:xcol; N=xcol ;
% elseif (nargin == 2),
%  N=xcol ;
% end;
% 
% if (N<0),
%  error('N must be greater than zero');
% end;
% 
% [trow,tcol] = size(t);
% if (trow~=1),
%   error('T must only have one row'); 
% end; 
% 
% N1 = length(x)+rem(N,2);
% length_FFT = N1; % Take an even value of N
% R = zeros(N,N);  
% for n = 0:N-1
%     M = min(n, N-1-n);
%     for k = 0:M
%         R(n+1, k+1) = x(n+k+1) * conj(x(n-k+1));
%     end
%     for k = N-1 : -1 : N-M
%         R(n+1, k+1) = conj(R(n+1, N-k+1));
%     end
% end
% tfr = zeros(N1, N1);
% for n = 0: N-1
%     temp = ifft(R(n+1,:), length_FFT)/(2*length_FFT);
%     temp = ifftshift(temp);
%     tfr(n+1, :) = temp; 
% end 
% 
% %f = -pi/(N*abs(t(2)-t(1)))*(N-1):2*pi/(N*abs(t(2)-t(1))):pi/(N*abs(t(2)-t(1)))*(N-1);%0 : N1-1;
% nn = -(N-1)/2:(N-1)/2;
% f = 2*nn*pi/(N*abs(t(2)-t(1)));
% 
% 
% %t = 0 : N1-1;