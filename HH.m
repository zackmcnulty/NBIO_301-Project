% Hodgkin/Huxley Equations (V_HH = -V-65), with current definition 
% of membrane potential (V=Vin-Vout)

vna=50;  %set the constants
vk=-77;
vl=-54.4;
gna=120;
gk=36;
gl=.3;
c=1;
I=6.35;

spikecount = 0;

v_init=-65;  %the initial conditions
m_init=.052;
h_init=.596;
n_init=.317;

npoints=50000;  %number of timesteps to integrate
dt=0.01;        %timestep

m=zeros(npoints,1); %initialize everything to zero
n=zeros(npoints,1);
h=zeros(npoints,1);
v=zeros(npoints,1);
time=zeros(npoints,1);

m(1)=m_init; %set the initial conditions to be the first entry in the vectors
n(1)=n_init;
h(1)=h_init;
v(1)=v_init;
time(1)=0.0;

tic
for step=1:npoints-1,
    v(step+1)=v(step)+((I - gna*h(step)*(v(step)-vna)*m(step)^3 ...
               -gk*(v(step)-vk)*n(step)^4-gl*(v(step)-vl))/c)*dt;
    m(step+1)=m(step)+ (alpha_m(v(step))*(1-m(step))-beta_m(v(step))*m(step))*dt;
    h(step+1)=h(step)+ (alpha_h(v(step))*(1-h(step))-beta_h(v(step))*h(step))*dt;
    n(step+1)=n(step)+ (alpha_n(v(step))*(1-n(step))-beta_n(v(step))*n(step))*dt;
    time(step+1)=time(step)+dt;
    
    if (v(step) > 1 & v(step+1) < v(step) & v(step-1) < v(step) & v(step)>20)
        spikecount = spikecount+1;
    end
end
toc

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
plot(time,v);
xlabel('t');
ylabel('V');

% figure(2)
% subplot(2,1,1);
% plot(time(10001:11500),v(10001:11500));
% xlabel('t');
% ylabel('V');
% subplot(2,1,2);
% hold on;
% plot(time(10001:11500),n(10001:11500));
% plot(time(10001:11500),m(10001:11500),'--');
% plot(time(10001:11500),h(10001:11500),':');
% xlabel('t');
% ylabel('n,m,h');
