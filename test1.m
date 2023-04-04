function test1
%%%%%%%%%%%%%%%%%%%%%
% ������ 1
%%%%%%%%%%%%%%%%%%%%%

a = 1;
x0 = 0;
xf = pi;
t0 = 0;
tf = 3;
Nx = 32;
Nt = 32;

% ������ ������� ��������� �������
pdep1ic = @(x,t) exp(t)*sin(x);
% ������ ������� ��������� �������
pdep1bc = @(t) [exp(t)*sin(x0); exp(t)*sin(xf)];
% ������ ������� ������������
tau = @(t) 1;
% ������ ������� ��������������
 pdep1het = @(x,t,ut) (2-exp(-1))*exp(t)*sin(x)+ut;
% ������ ������� ��� ���������:
u_exact = @(x,t) exp(t).*sin(x); 

x = linspace(x0,xf,Nx+1);
t = linspace(t0,tf,Nt+1);

[X, T, U] = pdep_fd_vd_CN(a, pdep1ic, pdep1bc, pdep1het, tau, x, t);
%[X, T, U] = pdep_fd_vd_R(a, pdep1ic, pdep1bc, pdep1het, tau, x, t);

% ������ ������ ���������� �����������
% mesh(X,T,U)
% %surf(X,T,U)
% title('u(x,t)')
% xlabel('Distance x')
% ylabel('Time t')
% shading interp

% ������� ������������ ����� �������� ������� � ������������� �������
% err = norm(u_exact(X,T) - U,1)
% ������� �������� �������� ������� � ������������� �������
err = max(max(abs(u_exact(X,T) - U)))
%hx = x(2)-x(1)
%ht = t(2)-t(1)
