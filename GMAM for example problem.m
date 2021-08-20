clc;clear all;
format short;
wx=1.5;   #forcing frequency
m=2;   # no of modes considered for approximation
K=[3,-2,0;-2,3,-1;0,-1,1];
M=diag([1,2,2]);
[V,D] = eigs(K,M,3,0);

tspan=[0:0.01:30];
fo=[5;0;0];    # vector of force amplitude
f=V'*fo;       # modal participation
x=zeros(m,size(tspan,2));
#Calculation of modal coordinates
for i=1:3
  x(i,:)=[f(i)/(D(i,i)-wx^2)*(sin(wx*tspan)-wx/sqrt(D(i,i))*sin(sqrt(D(i,i))*tspan))];
endfor
zo=V(:,1:m)*x(1:2,:);   #MDM
uexact=V*x;             #Exact solution 
Kinv=inv(K);
co=zeros(size(K));
for i=1:m
  co=co+V(:,i)*V(:,i)'/(D(i,i));
endfor
c1=Kinv-co;
y2=c1*fo*sin(wx(1)*tspan);   #first order correction
y3=c1*(M*Kinv)*fo*wx(1)^2*sin(wx(1)*tspan);  #Second order correction
y4=c1*(M*Kinv)^2*fo*wx(1)^4*sin(wx(1)*tspan);  #Third order correction
y5=c1*(M*Kinv)^3*fo*wx(1)^6*sin(wx(1)*tspan);  #Fourth order correction
##y6=c1*(M*Kinv)^4*fo*wx(1)^8*sin(wx(1)*tspan);
##y7=c1*(M*Kinv)^5*fo*wx(1)^10*sin(wx(1)*tspan);
##y8=c1*(M*Kinv)^6*fo*wx(1)^12*sin(wx(1)*tspan);
##y9=c1*(M*Kinv)^7*fo*wx(1)^14*sin(wx(1)*tspan);

## Newmark Beta method

#Initial conditions
u0=[0;0;0]; 
v0=[0;0;0];
C=zeros(size(M));
a0=inv(M)*(-C*v0-K*u0);  #Initial acceleration
tspan=[0:0.01:30];
dt=0.01;
u=[u0];
v=[v0];
a=[a0];
p=[5*sin(1.5*tspan)];
p=[p;zeros(2,size(p,2))];
for i=2:3001
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
plot(tspan,u(2,:),"linewidth",2);
hold on;
plot(tspan,zo(2,:),"linewidth",2);
err1=u-(zo);
err1=sum(diag(err1*err1'))
legend('Newmark beta','GMAM 0th order(MDM)')
xlabel('time')
ylabel('Reponse at 2nd node')
title('Response at 2nd node, Newmark Beta vs 0th order GMAM(MDM)')

figure;
plot(tspan,u(2,:),"linewidth",2);
hold on;
plot(tspan,(zo+y2)(2,:),"linewidth",2);
err2=uexact-(zo+y2);
err2=sum(diag(err2*err2'))
legend('Newmark beta','GMAM 1st order(MAM)')
xlabel('time')
ylabel('Reponse at 2nd node')
title('Response at 2nd node, Newmark Beta vs 1st order GMAM(MAM)')

figure;
plot(tspan,u(2,:),"linewidth",2);
hold on;
plot(tspan,(zo+y2+y3)(2,:),"linewidth",2);
err3=uexact-(zo+y2+y3);
err3=sum(diag(err3*err3'))
legend('Newmark beta','GMAM 2nd order')
xlabel('time')
ylabel('Reponse at 2nd node')
title('Response at 2nd node, Newmark Beta vs 2nd order GMAM')

figure;
plot(tspan,u(2,:),"linewidth",2);
hold on;
plot(tspan,(zo+y2+y3+y4)(2,:),"linewidth",2);
err4=uexact-(zo+y2+y3+y4);
err4=sum(diag(err4*err4'))
legend('Newmark beta','GMAM 3rd order(MAM)')
xlabel('time')
ylabel('Reponse at 2nd node')
title('Response at 2nd node, Newmark Beta vs 3rd order GMAM')
