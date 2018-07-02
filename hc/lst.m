clear;

theta=0:pi/2/1000:pi/2;

x=0.6*cos(theta)+0.2*cos(3.0*theta);
y=0.6*sin(theta)-0.2*sin(3.0*theta);

xd=-0.6*sin(theta)-0.6*sin(3.0*theta);
yd=0.6*cos(theta)-0.6*cos(3.0*theta);

x0(1:1001)=1.0;
y0(1:1001)=0.9;

z=(x-x0).*xd+(y-y0).*yd;

plot(theta,z);