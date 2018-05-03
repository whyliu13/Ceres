lo = 0;
hi = 1;
N = 32;
h = (hi-lo)/N;

x = 0:1/N:1;
y = 0:1/N:1;
z = zeros(N+1,N+1);

[X,Y] = meshgrid(x,y);

pi = 4.0*atan(1.0);

% x = 2.0*x0-1.0;
% y = 2.0*y0-1.0;
% 
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
%           sqrt(10.0 + 6.0*cos(4.0*beta))*0.2;
%     end
% end


num =2000;
for i = 1:num
 theta(i) = (i-1)*2.0d0*pi/num;
end

for i = 1:num
 xt(1,i) = 0.02d0*sqrt(5.0d0) + ...
           0.6d0*cos(theta(i)) + 0.2d0*cos(3.0d0*theta(i));
 xt(2,i) = 0.02d0*sqrt(5.0d0) + ...
           0.6d0*sin(theta(i)) - 0.2d0*sin(3.0d0*theta(i));
end

%plot(xt(1,:),xt(2,:),'.');

plot(0.5*(xt(1,:)+1.0d0),0.5*(xt(2,:)+1.0),'.r');





for j = 1:N
 for k = 1:N
 crit_dist = 0.0d0;
 dist1 = 1e+5;
 for i = 1:num
  dist2 = sqrt((x(j)-(xt(1,i)+1.0d0)/2.0d0)^2.0d0 + ...
               (y(k)-(xt(2,i)+1.0d0)/2.0d0)^2.0d0); 
  if(dist2 < dist1)
   crit_dist(1) = (xt(1,i)+1.0d0)/2.0d0;
   crit_dist(2) = (xt(2,i)+1.0d0)/2.0d0;
   dist1 = dist2;
  end
 end

 x0 = 2.0d0*x(j)-1.0d0;
 y0 = 2.0d0*y(k)-1.0d0;

 c1 = (0.02d0*sqrt(5.0d0)+1.0d0)/2.0d0;
 c2 = (0.02d0*sqrt(5.0d0)+1.0d0)/2.0d0;
 c3 = 0.02d0*sqrt(5.0d0);
 c4 = 0.02d0*sqrt(5.0d0);
 
 tt = atan((y(k)-c2)/(x(j)-c1));  
 if((x(j)-c1) >= 0.0d0 && (y(k)-c2) >= 0.0d0)
     
 elseif((x(j)-c1) <= 0.0d0 && (y(k)-c2) > 0.0d0)
    tt = tt + pi;
 elseif((x(j)-c1) < 0.0d0 && (y(k)-c2) < 0.0d0)
    tt = tt +pi;
 else
    tt = 2.0d0*pi + tt;
 end
 
%  dist3 =sqrt((0.5d0*(0.6d0*cos(tt) + 0.2d0*cos(3.0d0*tt)+c3+1.0d0)-c1)^2.0d0 + ...
%          (0.5d0*(0.6d0*sin(tt)-0.2d0*sin(3.0d0*tt)+c4+1.0d0)-c2)^2.0d0);
xxx(1) = 0.5d0*(0.6d0*cos(tt) + 0.2d0*cos(3.0d0*tt)+c3+1.0d0);
yyy(1) = 0.5d0*(0.6d0*sin(tt) - 0.2d0*sin(3.0d0*tt)+c4+1.0d0);
xxx(2) = x(j);
yyy(2) = y(k);
hold on;
plot(xxx,yyy,'bl');
hold on;

dist3 = sqrt((0.5d0*(0.6d0*cos(tt) + 0.2d0*cos(3.0d0*tt)))^2.0d0 + ...
         (0.5d0*(0.6d0*sin(tt)-0.2d0*sin(3.0d0*tt)))^2.0d0);
dist2 = sqrt((x(j)-c1)^2.0d0+(y(k)-c2)^2.0d0);
  
% ! c1 = (0.02d0*sqrt(5.0d0) + 1.0d0)/2.0d0
% ! c2 = (0.02d0*sqrt(5.0d0) + 1.0d0)/2.0d0
%  
% ! dist3= (crit_dist(1)-c1)*(x-crit_dist(1)) +&
% !        (crit_dist(2)-c2)*(y-crit_dist(2))

 ss = sign(dist3-dist2);

 z(j,k) = ss*dist1;
  
  end
 end
 
figure; 
%  v = [20,0];
% contour(X,Y,z,v);

contour(X,Y,z);

