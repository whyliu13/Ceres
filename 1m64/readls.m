clear;
lo = 0;
hi = 1;
N = 32;
h = (hi-lo)/N;

set(gcf, 'Position', [1000, 1000, 1000, 1000]);
% % % % % % % %sketch the grid
for i = 0:N
 A = zeros(2,2);
 A(1,2) = lo;
 A(2,2) = hi;
 A(:,1) = h*i;
 plot(A(:,1),A(:,2),'k');   
 hold on;   
end

for i = 0:N
 A = zeros(2,2);
 A(1,1) = lo;
 A(2,1) = hi;
 A(:,2) = h*i;
 plot(A(:,1),A(:,2),'k');   
 hold on;   
end
% % % % % % % % % % % % % % % % 


% % % % % % % read and sketch the linear recon
load ls1.dat;
load ls2.dat;
m1=size(ls1,1);
m2=size(ls2,1);

z=zeros(2,2);
for i=1:m1/2
  z(1,:)=ls1(1+(i-1)*2,:);
  z(2,:)=ls1(i*2,:);
 plot(z(:,1),z(:,2),'b');
 hold on;
end
for i=1:m2/2
  z(1,:)=ls2(1+(i-1)*2,:);
  z(2,:)=ls2(i*2,:);
 plot(z(:,1),z(:,2),'b');
 hold on;
end

% % % % % % % % % % % % plot the centroids
  b =  load('cen.dat');
  y1 = b(:,1);
  y2 = b(:,2);
  scatter (y1,y2,'.');
  
  hold on;
  c =  load('cen1.dat');
  y1 = c(:,1);
  y2 = c(:,2);
  scatter (y1,y2,'.r');
  
  hold on;
  d =  load('cen2.dat');
  y1 = d(:,1);
  y2 = d(:,2);
  scatter (y1,y2,'.g');






