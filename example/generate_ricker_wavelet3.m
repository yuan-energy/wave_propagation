peak_freq = 5;
time_leng =   1 ; 
N = 1000 ; 
[t1, acc1 ] = ricker_wavelet(peak_freq, time_leng, N) ; 

t2 = t1 + time_leng ;

% t = [t1 t2];
% acc = [acc1 acc1] ;

t3 = t1 + 2*time_leng ;
t = [t1 t2 t3];
acc_zero = zeros(length(acc1), 1) ;
acc = [acc1 acc_zero acc_zero] ;



plot(t, acc); 

time_step = t(2) - t(1)
[freq, amp] = time2freq(time_step, acc) ;

semilogx(freq, amp); 


fileID = fopen('ricker_acc.txt','w');


for	i = 1:length(t)
	fprintf(fileID,'%f %f\n',t(i) , acc(i));	
end

vel = integrate_vel(t, acc); 
dis = integrate_dis(t, acc, vel) ;


fileID = fopen('ricker_vel.txt','w');
for	i = 1:length(t)
	fprintf(fileID,'%f %f\n',t(i) , vel(i));	
end



fileID = fopen('ricker_dis.txt','w');
for	i = 1:length(t)
	fprintf(fileID,'%f %f\n',t(i) , dis(i));	
end
plot(t, vel)
plot(t, dis)

