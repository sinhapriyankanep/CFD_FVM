clc
close all
clear all

%% discretizing
N = 25;
L = 1;
dx = L / (N - 1);

dt = 0.1 * dx^2; % stable  using Fourier stability
tf = 5;
Nt = round(tf / dt);
t = linspace(0, tf, Nt);

%% Initializing
T_old = zeros(1, N);

% BCs
T_old(1) = 0;
T_old(end) = 1;
T_new = T_old;

%% initial condition
p1 = (N - 1) / L;
p2 = T_old(end) - T_old(1);
for i = 2:N-1
    T_old(i) = ((i-1)*dx)^2;
end

%%  for length of L
x = linspace(0, L, N);
tol = 1e-7;  
figure;
hold on

conv = 0;  
%looping for time
for n = 1:Nt
    for i = 2:N-1
        T_new(i) = T_old(i) + (dt / dx^2) * (T_old(i+1) + T_old(i-1) - 2*T_old(i));
    end
    T_old = T_new;
    residual = sum(abs(T_new - T_old)); 
    %checking convergencene
    if residual < tol && conv == 0
        conv = n * dt; 
    end

  %plotting
color = {'r-', 'b--', 'g-.'}; 

mid = round(Nt / 2); 

%plotting
if n == 1
    plot(x, T_new, color{1}, 'DisplayName', sprintf('t = %.2f s', n * dt));
    elseif n == 100
    plot(x, T_new, 'k--','LineWidth',2, 'DisplayName', sprintf('t = %.2f s', n * dt));
elseif n == mid
    plot(x, T_new, color{2},'LineWidth',2, 'DisplayName', sprintf('t = %.2f s', n * dt));
elseif n == Nt
    plot(x, T_new, color{3}, 'DisplayName', sprintf('t = %.2f s', n * dt));
end
ylabel('Normalized Temperature')
xlabel('Normalized length ')
legend('show'); 
end
