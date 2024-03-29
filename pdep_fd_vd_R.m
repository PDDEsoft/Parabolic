function [X, T, U] = pdep_fd_vd_R(a, pdep1ic, pdep1bc, pdep1het, tau, x, t)
% ������� � ����������� ������������� ����� ������ ����� � ������ � ���������� �� ������ ����������.
% ����� ��� ��������� ��� ��������:
%
% u_t(x,t) = a^2 u_xx(x,t) + f(x,t,u(x,t-tau(t)))
%
% �������� ������� ����������:
% a - ����������� ��� 2-� ������� ����������� �� x;
% pdep1ic - ������� �� ���� ���������� (x � t), �������� ��������� �������;
% pdep1bc - ������-������� �� ����� ���������� (t), �������� ��������� 
%           �������: ������ - �� ������ ������� �� ���������� x, ������ - �� �����;
% pdep1het - ��������� ������� ���� ���������� f(x,t,u(x,t-tau(t))),
%            ���������� ������������ ������ ���������;
% tau - ��������������� ������� ������������ ����� ����������;  
% x - ������, �������� ���������� �������������� ���������������� ���������� �����;
% t - ������, �������� ���������� �������������� ��������� ���������� �����.

% ������ ������� ������-�������� � ��������� ������
[X, T, U] = pdep_fd_vd_CN(a, pdep1ic, pdep1bc, pdep1het, tau, x, t);
% ������ ������� ������-�������� � ���������� ������
Nx = length(x)-1;
Nt = length(t)-1;
x2 = linspace(x(1),x(end),2*Nx+1);
t2 = linspace(t(1),t(end),2*Nt+1);
[~, ~, U2] = pdep_fd_vd_CN(a, pdep1ic, pdep1bc, pdep1het, tau, x2, t2);
U = 4/3 * U2(1:2:end,1:2:end) - U/3;