% Galerkin Spectral Method
% Bryan Kaiser
% 11/6/14

close all
clear all
clc

% Solution to:
% du/dt = a^2*(d^2u/dx^2)+sin(5x)
% by the Galerkin method

%---------------------------------------------------------------

Nc = 50; % Number of coeffs
Nx = 60; % Number of spatial points
Nt = 40; % Timesteps
T = 10; % Final time
t = 0; % Initial time
dt = T/Nt;

%---------------------------------------------------------------

a = sqrt(0.5); % A diffusion constant
x = linspace(0,pi,Nx);
As = zeros(Nc,Nx); % An(t)*sin(n*x)
As5 = zeros(1,Nc); % A5(t)*sin(5*x)
m = 1; % Movie index
IC = x.*(pi-x);

while t < T+dt
 
% Coefficients * sin(n*x) = An(t)*sin(n*x):
for n = 1:Nc
    if n ~= 5
        for i = 1:Nx 
        As(n,i) = ((4*(1-cos(pi*n)))/(pi*n^3))*exp(-t*(n*a)^2)*sin(n*x(i));
        end
    elseif n == 5
        for i = 1:Nx
        As5(i) = (((4*(1-cos(pi*5)))/(pi*5^3)-pi/(50*a^2))*exp(-t*(5*a)^2)+pi/(50*a^2))*sin(5*x(i));
        end
    end
end

% Solution:
u = sum(As)+As5;

% Movie
hfig = figure(1);
set(hfig,'Position',[300,300,1400,1000]);
plot(x,IC,'k--','LineWidth',5);
hold on
plot(x,u)
xlabel('u')
ylabel('x')
axis([0,pi,-0.5,pi^2/4+0.1])
hold on
M(m) = getframe;

% Advancement
m = m+1;
t = t+dt;
end
