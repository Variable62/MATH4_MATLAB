%% 
% *66010892 อดิศร สมมาตร*
% 
% 

clc;clf;clear;
close all;
syms t w;
f = piecewise(t < -1, -t-1, t >= -1 & t < 1, 2*t, t >= 1, -t+1);
%% 
% *ก) Plot the original signal*

figure;
fplot(f, [-2, 2], 'LineWidth', 1.5);
grid on;
title('Original Signal');
xlabel('t');
ylabel('f(t)');
legend('f(t)');
%% 
% *Fourier Series Approximation*

L = 4;  % Period
n = 10; % Number of terms 
a0 = (1/L) * int(f, t, -2, 2);
fourier_series = a0 / 2;
for k = 1:n
    ak = (2/L) * int(f * cos(k * pi * t / 2), t, -2, 2);
    bk = (2/L) * int(f * sin(k * pi * t / 2), t, -2, 2);
    fourier_series = fourier_series + ak * cos(k * pi * t / 2) + bk * sin(k * pi * t / 2);
end
F_w = fourier(f, t, w);
%% 
% *Plot Fourier Series Approximation* 

figure;
fplot(fourier_series, [-6, 6], 'r', 'LineWidth', 1.5);
grid on;
title('Fourier Series Approximation (Periodic Extension)');
xlabel('t');
ylabel('f(t)');
legend('Fourier Series Approximation');
%% 
% *Show result Fourier Series Expression*

disp('Fourier Series Expression:');
disp(vpa(fourier_series, 4)); 
%% 
% *ข)*
% 
% *Show result Fourier Transform Expression*

disp('Fourier Transform Expression:');
disp(vpa(F_w, 4));  
%% 
% 

clearvars; clf; clc;
L = 2;               
N = 5;               
dx = L / (N-1);       
x = linspace(0, L, N); 
T = 1;                
dt = 0.01;           
t = 0:dt:T;        
u = zeros(N, length(t));
u(:,1) = sin(x/2) + 10*x + 50;
alpha = dt / (dx^2);
for n = 1:length(t)-1
    for i = 2:N-1
        u(i,n+1) = u(i,n) + 4*alpha * (u(i+1,n) - 2*u(i,n) + u(i-1,n)) + ...
                   dt * sin(x(i)/2); 
    end   
    u(1,n+1) = 50; 
    u(N,n+1) = u(N-1,n+1) + dx * 10;
end
%% 
% *Plot 3D graph*

figure;
mesh(t, x, u);
xlabel('Time (t)');
ylabel('Space (x)');
zlabel('u(x,t)');
title('Graph Equation Solution');
%% 
% 

clearvars; clf; clc;
L = 3;          
T = 10;             
N = 20;             
M = 1000;           
dx = L / (N-1);    
dt = T / M;        
x = linspace(0, L, N); 
t = linspace(0, T, M); 
u = zeros(N, M);    
u(:,1) = 3;        
alpha = dt / (9 * dx^2); 
for n = 1:M-1
    for i = 2:N-1
        u(i, n+1) = u(i, n) + alpha * (u(i+1, n) - 2*u(i, n) + u(i-1, n));
    end
    u(1, n+1) = u(1, n) + alpha * (u(2, n) - u(1, n));  
    u(N, n+1) = u(N-1, n+1); 
end
%% 
% *Plot graph 3D*

figure;
mesh(t, x, u);
xlabel('Time (t)');
ylabel('Position (x)');
zlabel('u(x,t)');
title('Solution of the Heat Equation');
%% 
% 

clearvars;clf;clc;
syms u(x,t) k
%% 
% *initial condition u(x,0) as piecewise*

u_initial = piecewise(0 <= x & x <= 10, x - 10, x > 10, 0);
u_hat_k = fourier(u_initial, x, k);
%% 
% *solve ODE*

u_hat_k_t = u_hat_k * exp(-k^2 * t / 4);
u_solution = ifourier(u_hat_k_t, k, x);
%% 
% *Display the solution*

disp('The solution u(x,t) is:')
disp(u_solution)
%% 
% 

clearvars; clf; clc;
L = pi;               
T = 10;               
n = 5;             
m = 500;             
dx = L / (n-1);      
dt = T / m;            
x = linspace(0, L, n);
t = linspace(0, T, m);
u = zeros(n, m);       
u(:,1) = sin(x) + (3*pi/2 - x);  
u(:,2) = sin(x) + (pi - x);      
for k = 2:m-1
    for i = 2:n-1
        u(i, k+1) = 2*(1 - 2*dt^2/dx^2)*u(i,k) - u(i,k-1) + (dt^2 * sin(x(i))) + (dt^2 / dx^2) * (u(i+1, k) - 2*u(i,k) + u(i-1, k));
    end
    u(1, k+1) = u(1, k);  
    u(n, k+1) = 0;        
end
%% 
% *Plot 3D*

figure;
mesh(t, x, u);
xlabel('Time (t)');
ylabel('Position (x)');
zlabel('u(x,t)');
title('Solution of the Wave Equation');
%% 
% 

clearvars;clc;clf;
%% 
% 

clearvars;
syms u(x,t) k h lambda;
%% 
% 

clearvars;clf;clc;
L = 1; 
n = 5;
u0 = @(y) cos(y); 
u1 = @(y) 1 + cos(2 * y); 
du0 = @(x) 0; 
du1 = @(x) 0; 
[x, y] = meshgrid(linspace(0, L, n), linspace(0, L, n));
u = zeros(n, n);
for i = 1:n
    for j = 1:n
        u(i, j) = u0(y(j)) + (u1(y(j)) - u0(y(j))) * (x(i) / L);
    end
end
%% 
% *Plot graph*

figure;
surf(x, y, u);
title('ผลเฉลยของสมการ U_{xx} + U_{yy} = 0');
xlabel('x');
ylabel('y');
zlabel('u(x, y)');
%% 
% 

clearvars;clc;clf
syms x y k A B
u_x0 = 0;  % u(x,0) = 0
u_x1 = exp(-x);  % u(x,1) = e^(-x)
%% 
% *Fourier transform of u(x,1)*

U_k = fourier(u_x1, x, k);
U_k_y = A*exp(-k*y) + B*exp(k*y);
eq1 = A + B == 0;  
eq2 = subs(U_k_y, y, 1) == U_k;  
sol = solve([eq1, eq2], [A, B]);
U_k_y_solution = subs(U_k_y, {A, B}, {sol.A, sol.B});
%% 
% *Inverse Fourier transform to get u(x, y)*

u_xy = ifourier(U_k_y_solution, k, x);
disp(u_xy)
%% 
% 

clearvars;clc;
syms k lambda

eq = sin(k) - 2*k*cos(k);
k_values = [];
for n = 1:5  
    k_guess = n*pi; 
    k_sol = vpasolve(eq, k, k_guess);  
    if ~isempty(k_sol) && all(~ismember(k_sol, k_values)) 
        k_values = [k_values, k_sol];  
    end
end
%% 
% *eigenvalues lambda*

lambda_values = (k_values.^2) / 9; 
%% 
% *Display eigenvalues*

disp('Eigenvalues (lambda):');
disp(lambda_values);
eigenfunctions = [];
for i = 1:length(lambda_values)
    eigenfunction = sin(sqrt(9*lambda_values(i))*sym('x'));
    eigenfunctions = [eigenfunctions; simplify(eigenfunction)];
end
%% 
% *Display eigenfunctions*

disp('Eigenfunctions:');
disp(eigenfunctions);
%% 
%