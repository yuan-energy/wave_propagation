function [ times, ampli ] = ricker_wavelet( peak_frequency, time_leng,  N )
%RICKER_WAVELET return time series of ricker wavelet
%   Detailed explanation goes here
% f, length=0.128, dt=0.001
	% default arguments
	if nargin < 3
	   peak_frequency =   2 ;
	   time_leng =   2 ; 
	   N =   1000 ; 
	end

	f= peak_frequency ; % peak frequency
	times = linspace( -(time_leng/2), time_leng/2 , N+1) ;
	% times 
	tmp = pi ^2 * f^2 * times.^2 ; 
	ampli = (1 - 2 * tmp) .* exp(- tmp) ;

	times = times + time_leng/2 ;
end

