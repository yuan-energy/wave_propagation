
%% time2freq: function description
function [freqs, amplitudes] = time2freq(time_step, disp_time_series)


% time = northridge(:,1);
% displ = northridge(:,4);
% fig = plot(time, displ, 'Linewidth',3);
% grid;
% xlabel('Time (s)');
% ylabel('Displacement (cm)');
% title('Northridge Time Series')
% saveas(fig,'Northridge_time_disp.jpg');

	Y=fft(disp_time_series);
	L=length(disp_time_series);
	P2= abs(Y/L) ;
	P1=P2(1:floor(L/2)+1);
	P1(2:end-1) = 2*P1(2:end-1);
	Fs=1/time_step ;
	f=Fs*(0:(L/2))/L;

	freqs = f;
	amplitudes = P1 ;

	% fig = semilogx(freqs,amplitudes,'Linewidth',3);
	% grid;
	% xlabel('Frequency (Hz)');
	% ylabel('Displacement (cm)');
	% title('Northridge Frequency Analysis')
	% saveas(fig,'Northridge_frequency_disp.jpg');