clear;
lo = 0;
hi = 1;
N = 32;
h = (hi-lo)/N;

mm=[];

set(gcf, 'Position', [1000, 1000, 1000, 1000]);
% % % % % % % %sketch the grid
for i = 0:N
 A = zeros(2,2);
 A(1,2) = lo;
 A(2,2) = hi;
 A(:,1) = h*i;
 plot(A(:,1),A(:,2),'k','HandleVisibility','off');   
 hold on;   
end

for i = 0:N
 A = zeros(2,2);
 A(1,1) = lo;
 A(2,1) = hi;
 A(:,2) = h*i;
 plot(A(:,1),A(:,2),'k','HandleVisibility','off');   
 hold on;   
end
% % % % % % % % % % % % % % % % 


% % % % % % % read and sketch the linear recon
load ls1.dat;
load ls2.dat;
load ls3.dat;
load ls4.dat;
m1=size(ls1,1);
m2=size(ls2,1);
m3=size(ls3,1);
m4=size(ls4,1);

z=zeros(2,2);
for i=1:m1/2
  z(1,:)=ls1(1+(i-1)*2,:);
  z(2,:)=ls1(i*2,:);
 mm(6)=plot(z(:,1),z(:,2),'b');
 
 hold on;
end
for i=1:m2/2
  z(1,:)=ls2(1+(i-1)*2,:);
  z(2,:)=ls2(i*2,:);
 plot(z(:,1),z(:,2),'b');
 hold on;
end

for i=1:m3/2
  z(1,:)=ls3(1+(i-1)*2,:);
  z(2,:)=ls3(i*2,:);
%  plot(z(:,1),z(:,2),'b');
 hold on;
end

for i=1:m4/2
  z(1,:)=ls4(1+(i-1)*2,:);
  z(2,:)=ls4(i*2,:);
 plot(z(:,1),z(:,2),'b');
 hold on;
end

% % % % % % % % % % % % plot the centroids
  b =  load('cen.dat');
  y1 = b(:,1);
  y2 = b(:,2);
  mm(1)=scatter (y1,y2,'.');
  
  hold on;
  c =  load('cen1.dat');
  y1 = c(:,1);
  y2 = c(:,2);
  mm(2)=scatter (y1,y2,'.r');
  
  hold on;
  d =  load('cen2.dat');
  y1 = d(:,1);
  y2 = d(:,2);
  mm(3)=scatter (y1,y2,'.g');
  
  
    hold on;
  e =  load('cen3.dat');
  y1 = e(:,1);
  y2 = e(:,2);
  mm(4)=scatter (y1,y2,'.k');
%   
     hold on;
  h =  load('cen4.dat');
  y1 = h(:,1);
  y2 = h(:,2);
  mm(5)=scatter (y1,y2,'.m');
hold on;

 legend(mm,'Material1','Material2','Material3','Material4','Material5','Reconstructed Interface');







