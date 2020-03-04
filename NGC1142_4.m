c=299792.458;       % ligth speed km/s
Ho=70;              % Hubble const. km/s/Mpc
G = 6.67*10^(-11);  %Grav. const.
fun = @(z)(0.27*(1+z).^3+1-0.27).^(-0.5);
% v = c*((1+z)^2-1)/((1+z)^2+1)

ra0 = 43.795935*pi/180;
dec0 = -0.180749*pi/180;
r0 = 124.38;

load 'neigh_5Mpc.txt' %nearest galaxies 2x2 sq. deg.
env(1,:) = neigh_5Mpc(:,1)*pi/180; %ra
env(2,:) = neigh_5Mpc(:,2)*pi/180; %dec
env(3,:) = neigh_5Mpc(:,4);        %redshift
[~,sz] = size(env);
for i=1:sz
    x(i) = r0 *tan(env(2,i)-dec0);
    y(i) = r0 *tan(env(1,i)-ra0);
end

env(4,:) = neigh_5Mpc(:,5);
env(5,:) = neigh_5Mpc(:,3);         %velocity

for i=1:sz
    env(6,i) = sqrt(x(i)^2+y(i)^2); %projection dist.
end
[~, N] = size(env);
m = 0;
for i=1:sz 
   m = m + (env(5,i)*1000)^2* env(6,i)*3.086*10 ^ 19; %sum for eq.
end

M = 32*m/(pi*G*(N-1.5)*2*10^30); %final result

%systematic err.
d_M = 0;
for i=1:sz
    d_M = d_M + (2*env(5,i)*env(6,i)*env(4,i)*3.086*10 ^ 19)^2;   
end

d_M = 32*d_M^0.5/(pi*G*(N-1.5)*2*10^30);