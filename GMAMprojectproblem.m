close all;
format short;
close all;
##First forcing function

wx=25000;   #Forcing frequency
m=19;       #Modes considered for approximation
load('Project1_2_GMAM_GMTA_Matlab.mat')
[V,D] = eigs(K,M,m,0);

#Vector of force amplitude
fo=zeros(size(K,1),1);
fo(1891*2-1)=10^5;

tspan=[0:0.00001:0.0008];
f=V'*fo;      #Modal participation
x=zeros(m,size(tspan,2));
for i=1:m
  x(i,:)=[f(i)/(D(i,i)-wx^2)*(sin(wx*tspan)-wx/sqrt(D(i,i))*sin(sqrt(D(i,i))*tspan))];
endfor
zo=V(:,1:m)*x;   #MDM

Kinv=inv(K);
co=zeros(size(K));
for i=1:m
  co=co+V(:,i)*V(:,i)'/(D(i,i));
endfor
c1=Kinv-co;
y2=c1*fo*sin(wx(1)*tspan);    #first order correction
y3=c1*(M*Kinv)*fo*wx(1)^2*sin(wx(1)*tspan);  #Second order correction
y4=c1*(M*Kinv)^2*fo*wx(1)^4*sin(wx(1)*tspan); #Third order correction
##y5=c1*(M*Kinv)^3*fo*wx(1)^6*sin(wx(1)*tspan);
##y6=c1*(M*Kinv)^4*fo*wx(1)^8*sin(wx(1)*tspan);
##y7=c1*(M*Kinv)^5*fo*wx(1)^10*sin(wx(1)*tspan);
##y8=c1*(M*Kinv)^6*fo*wx(1)^12*sin(wx(1)*tspan);
##y9=c1*(M*Kinv)^7*fo*wx(1)^14*sin(wx(1)*tspan);

##Second forcing function

wy=10000;   #forcing frequency
fo1=zeros(size(K,1),1);
fo1(1891*2)=10^5;  #vector of force amplitude
f1=V'*fo1;         #Modal participation
m=4;               #No of modes considered for approximation 

x=zeros(m,size(tspan,2));
for i=1:m
  x(i,:)=[f1(i)/(D(i,i)-wy^2)*(sin(wy*tspan)-wy/sqrt(D(i,i))*sin(sqrt(D(i,i))*tspan))];
endfor
zo1=V(:,1:m)*x;     #MDM

co1=zeros(size(K));
for i=1:m
  co1=co1+V(:,i)*V(:,i)'/(D(i,i));
endfor
c11=Kinv-co1;
y21=c11*fo1*sin(wy(1)*tspan);    #first order correction
y31=c11*(M*Kinv)*fo1*wy(1)^2*sin(wy(1)*tspan);   #Second order correction
y41=c11*(M*Kinv)^2*fo1*wy(1)^4*sin(wy(1)*tspan); #Third order correction
##y51=c11*(M*Kinv)^3*fo1*wy(1)^6*sin(wy(1)*tspan);
##y61=c11*(M*Kinv)^4*fo1*wy(1)^8*sin(wy(1)*tspan);
##y71=c11*(M*Kinv)^5*fo1*wy(1)^10*sin(wy(1)*tspan);
##y81=c11*(M*Kinv)^6*fo1*wy(1)^12*sin(wy(1)*tspan);
##y91=c11*(M*Kinv)^7*fo1*wy(1)^14*sin(wy(1)*tspan);

figure;
plot(tspan,(zo+zo1)(3781,:),"linewidth",2);
legend('GMAM 0th order(MDM)')
xlabel('time')
ylabel('Reponse at 1891st node in xdir')
title('Response at 1891st node in xdir, 0th order GMAM(MDM)')

figure;
plot(tspan,(zo+y2+ zo1+y21)(3781,:),"linewidth",2);
legend('GMAM 1st order(MAM)')
xlabel('time')
ylabel('Reponse at 1891st node in xdir')
title('Response at 1891st node in xdir, 1st order GMAM(MAM)')

figure;
plot(tspan,(zo+y2+y3 +zo1+y21+y31)(3781,:),"linewidth",2);
legend('GMAM 2nd order')
xlabel('time')
ylabel('Reponse at 1891st node in xdir')
title('Response at 1891st node in xdir, 2nd order GMAM')

figure;
plot(tspan,(zo+y2+y3+y4 +zo1+y21+y31+y41)(3781,:),"linewidth",2);
legend('GMAM 3rd order(MAM)')
xlabel('time')
ylabel('Reponse at 1891st node in xdir')
title('Response at 1891st node in xdir, 3rd order GMAM')
