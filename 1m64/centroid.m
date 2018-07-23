clear;

lo = 0;
hi = 1;
N = 32;
h = (hi-lo)/N;

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
% 
%   a =  load('lineseg10.dat');
%   for j= 1:size(a,1)
%    x1(1)= a(j,1);
%    x2(1)= a(j,2);
%    x1(2)= a(j,3);
%    x2(2)= a(j,4);
%    plot(x1,x2,'r');
%    axis([0,100,0,100]);
%    hold on;
%   end
  
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
  
  
    hold on;
  e =  load('cen3.dat');
  y1 = e(:,1);
  y2 = e(:,2);
  scatter (y1,y2,'.k');
%   
     hold on;
  h =  load('cen4.dat');
  y1 = h(:,1);
  y2 = h(:,2);
  scatter (y1,y2,'.m');

hold on;



% x0 = 0+1/(2*N):1/N:1-1/(2*N);
% y0 = 0+1/(2*N):1/N:1-1/(2*N);
% z = zeros(N,N);


 x0 = 0:1/N:1;
 y0 = 0:1/N:1;
z = zeros(N+1,N+1);

[X,Y] = meshgrid(x0,y0);




pi = 4.0*atan(1.0);

% x = 2.0*x0-1.0;
% y = 2.0*y0-1.0;
% 
% % c1 = 0.0d0;
% % c2 = 0.0d0;
%  c1 = 0.02*sqrt(5.0);
%  c2 = 0.02*sqrt(5.0);
% for i = 1:N
%     for j = 1:N
%      beta = atan((y(j)-c1)/(x(i)-c2));  
%   if(x(i)== c1 && y(j)== c2) 
%      beta = 0;
%   end    
%   if((x(i)-c1)>= 0 && (y(j)-c2) >=0 )
%     
%   elseif((x(i)-c1)<=0 && (y(j)-c2)>0 )
%     beta = beta+pi;
%   elseif((x(i)-c1)<0 && (y(j)-c2)<0 )
%     beta = beta +pi;
%   else
%     beta = 2.0*pi + beta;
%   end 
%  z(i,j) = sqrt((x(i)-c1)^2 + (y(j)-c2)^2) - ...
%           (0.5 + 0.2*sin(5.0*beta));
%     end
% end



% radcen = 0.25;
% radeps = 0.005;
% 
% 
% x= x0;
% y =y0;
% 
% for i = 1:64
%    for j=1:64
%  center1 = 0.5;
%  center2 = 0.5;
%  r1=radcen-radeps;
%  r2=radcen+radeps;
%  dist1= sqrt((x(i)-center1)^2.0+(y(j)-center2)^2.0);
%  z(i,j) = r1-dist1;
% % 
% %    if(dist1 <= r2 &&  dist1 >= r1)
% %       dist2 = r2 - dist1;
% %       dist3 = dist1 - r1;
% %       z(i,j) = min(dist2,dist3);
% %    elseif(dist1 >= r2)
% %     z(i,j)=r2-dist1;
% %    elseif(dist1 <= r1)
% %     z(i,j)=dist1-r1;
% %    end      
%    end
% end
% 
% 
% % figure;
% v = [1,0];
% 
% contour(X,Y,transpose(z),v);
% 
% hold on;
% for i = 1:64
%    for j=1:64
%  center1 = 0.5;
%  center2 = 0.5;
%  r1=radcen-radeps;
%  r2=radcen+radeps;
%  dist1= sqrt((x(i)-center1)^2.0+(y(j)-center2)^2.0);
%  z(i,j) = dist1-r2;
%    end
% end


% figure;
v = [100,0];
load 'levelset.dat';
%contour(X,Y,transpose(z),v);
 contour(X,Y,levelset,v);
hold on; 
load 'levelset1.dat';
contour(X,Y,levelset1,v);

hold on;
 load 'levelset2.dat';
 contour(X,Y,levelset2,v);





figure;
load 'vf.dat';
% for i = 1:N
%  for j= 1:N
%    if(vf(i,j) < 1.0 && vf(i,j)>0.0)
%       vf(i,j) = 1.0;
%    else
%       vf(i,j) = 0.0;
%    end
%  end
% end
hold on;
surf(X,Y,vf);
view(2);
%contour(X,Y,vf);



  
 