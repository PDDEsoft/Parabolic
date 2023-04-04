function test1
%%%%%%%%%%%%%%%%%%%%%
% Пример 1
%%%%%%%%%%%%%%%%%%%%%

a = 1;
x0 = 0;
xf = pi;
t0 = 0;
tf = 3;
Nx = 32;
Nt = 32;

% Задаем функцию начальных условий
pdep1ic = @(x,t) exp(t)*sin(x);
% Задаем функцию граничных условий
pdep1bc = @(t) [exp(t)*sin(x0); exp(t)*sin(xf)];
% Задаем функцию запаздывания
tau = @(t) 1;
% Задаем функцию неоднородности
 pdep1het = @(x,t,ut) (2-exp(-1))*exp(t)*sin(x)+ut;
% точное решение для сравнения:
u_exact = @(x,t) exp(t).*sin(x); 

x = linspace(x0,xf,Nx+1);
t = linspace(t0,tf,Nt+1);

[X, T, U] = pdep_fd_vd_CN(a, pdep1ic, pdep1bc, pdep1het, tau, x, t);
%[X, T, U] = pdep_fd_vd_R(a, pdep1ic, pdep1bc, pdep1het, tau, x, t);

% Строим график трехмерной поверхности
% mesh(X,T,U)
% %surf(X,T,U)
% title('u(x,t)')
% xlabel('Distance x')
% ylabel('Time t')
% shading interp

% Находим постолбцовую норму разности точного и приближенного решений
% err = norm(u_exact(X,T) - U,1)
% Находим максимум разности точного и приближенного решений
err = max(max(abs(u_exact(X,T) - U)))
%hx = x(2)-x(1)
%ht = t(2)-t(1)
