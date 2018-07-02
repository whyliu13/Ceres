clear all;
load output2.dat;
N = 32;
step = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1/N;
x = 1/N/2:1/N:1-1/N/2;
y=x;
A = zeros(N,N);
for i  =1:step
A = output2(1+(i-1)*N:N+(i-1)*N,2:33);
pcolor(x,y,A);

% view(2);
 eval(['print -djpeg step' num2str(i) '.jpg']);
end