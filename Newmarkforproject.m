clear all;
clc;close all;
load('Project1_2_GMAM_GMTA_Matlab.mat')

#Initial condition
u0=zeros(3782,1);
v0=zeros(3782,1);
C=zeros(size(M));

#Initial acceleratio
a0=inv(M)*(-C*v0-K*u0);
tspan=[0:0.00001:0.0008];
dt=0.00001;

u=[u0];
v=[v0];
a=[a0];
#forcing function
p=[10^5*sin(25000*tspan);10^5*sin(10000*tspan)];
p=[zeros(3780,size(p,2));p];

for i=2:81
  upred=u(:,i-1)+dt*v(:,i-1)+dt^2/2*a(:,i-1)/2;
  vpred=v(:,i-1)+dt*a(:,i-1)/2;
  RHS=p(:,i)-C*vpred-K*upred;
  anp=inv(M+1/2*dt*C+1/4*dt^2*K)*RHS;
  a=[a,anp];
  unp=upred+1/4*dt^2*anp;
  u=[u,unp];
  vnp=vpred+1/2*dt*anp;
  v=[v,vnp];
endfor

figure;
plot(tspan,u(3781,:),"linewidth",2);
legend('Newmark beta')
xlabel('time')
ylabel('Reponse at 1891st node in xdir')
title('Response at 1891st node in xdir, Newmark Beta')
