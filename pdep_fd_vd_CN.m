function [X, T, U] = pdep_fd_vd_CN(a, pdep1ic, pdep1bc, pdep1het, tau, x, t)
% ������� � ����������� ������������� ����� ������ ������-��������.
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

% �������� ������������ ������� ������
if feval(tau,t(1)) < 0
    error('������� ������������ ������ ���� ���������������')
end
temp_diff = diff(x);
if any(abs(temp_diff-temp_diff(1)) > 10*eps)
    error('������ x ������ ��������� ���������� �������������� ��������')
end
if x(end) - x(1) < 0
    error('������ x ������ ��������� ������������� �� ����������� ��������')
end
if length(x) < 3
    error('������ x ������ ��������� �� ����� 3-� ���������')
end
temp_diff = diff(t);
if any(abs(temp_diff-temp_diff(1)) > 10*eps)
    error('������ t ������ ��������� ���������� �������������� ��������')
end
if t(end) - t(1) < 0
    error('������ t ������ ��������� ������������� �� ����������� ��������')
end
if length(t) < 2
    error('������ t ������ ��������� �� ����� 2-� ���������')
end
if any(feval(pdep1bc,t(1)) - [feval(pdep1ic,x(1),t(1)); feval(pdep1ic,x(end),t(1))] > 10*eps)
    error('����� ���������� � �������� ��������� ��������� ����������� �������')
end
    
nx = length(x); nt = length(t);
hx = (x(end)-x(1))/(nx-1); ht=(t(end)-t(1))/(nt-1); 
k = ht*a*a/(hx*hx);
s = 0.5; % �������� ������ ����� (�� ������!)
U = zeros(nx,nt); % ������� �������
for i=1:nt
    U([1,nx],i) = feval(pdep1bc,t(i));
end
for j=1:nx
    U(j,1) = feval(pdep1ic,x(j),t(1));
end
temp_u = U; % ��������������� �������
A = zeros(1,nx-2); % ������������ ������ ��������
B = zeros(1,nx-2); % ������������ ������ ��������

% ����� ����� � �������������� ��� f, � ������������� � ��������������
for i=2:nt
    for j=2:(nx-1)
        tcur = t(1)+(i-1.5)*ht;
        % ��������� ������� ������������ u(x,t-tau(t)) � ������� ����� tcur
        % �� �������� ��� ������ x-��������� ����������� � ���������� del_cur
        if (tcur - feval(tau,tcur) <= t(1))  % ������, ����� t-tau(t) �������� � ������� ��������� �������
            %del_cur1 = feval(pdep1ic,x(1)+(j-2)*hx,tcur-feval(tau,tcur));
            del_cur2 = feval(pdep1ic,x(1)+(j-1)*hx,tcur-feval(tau,tcur));
            %del_cur3 = feval(pdep1ic,x(1)+j*hx,tcur-feval(tau,tcur));
        elseif (0.5*ht - feval(tau,tcur) >= 0)  % �������� �������������
            error('not corrected yet for cubic polynoms')
            if (i == 2)
                %del_cur1 = U(j-1,1)+(0.5 - feval(tau,tcur)/ht)*(U(j-1,1)-feval(pdep1ic,x(1)+(j-2)*hx,t(1)-ht));
                del_cur2 = U(j,1)+(0.5 - feval(tau,tcur)/ht)*(U(j,1)-feval(pdep1ic,x(1)+(j-1)*hx,t(1)-ht));
                %del_cur3 = U(j+1,1)+(0.5 - feval(tau,tcur)/ht)*(U(j+1,1)-feval(pdep1ic,x(1)+j*hx,t(1)-ht));
            else
                %del_cur1 = U(j-1,i-1)+(0.5 - feval(tau,tcur)/ht)*(U(j-1,i-1)-U(j-1,i-2));
                del_cur2 = U(j,i-1)+(0.5 - feval(tau,tcur)/ht)*(U(j,i-1)-U(j,i-2));
                %del_cur3 = U(j+1,i-1)+(0.5 - feval(tau,tcur)/ht)*(U(j+1,i-1)-U(j+1,i-2));
            end
        else
            ind = floor((tcur-feval(tau,tcur)-t(1))/ht)+1;
            % �������� ������������
            %del_cur2 = U(j,ind)+((i-ind-0.5)-feval(tau,tcur)/ht)*(U(j,ind+1)-U(j,ind));
            % ���������� ������������ (������ ��������� �� ����� � ��������� ind-1, ind, ind+1, ind+2)
            ti = tcur-feval(tau,tcur);
            % ���������� 4 ����� ��� ����������� ����������������� ��������
            if (ind+2<i)&&(ind-1>=1)
                t_int = (ind-2:ind+1)*ht+t(1);
                u_int = U(j,ind-1:ind+2);
            elseif (i>4)&&(ind-1<1)
                t_int = (0:3)*ht+t(1);
                u_int = U(j,1:4);
            elseif ind+2>=i
                error('�������, ����� ������������ ������ ����, ���� �� ����������')
            elseif i<4
                error('In current version of program time step must be decreased')
            end
            tm = ti-t_int; 
            ht3 = ht^3;
            del_cur2 = u_int(4)*tm(1)*tm(2)*tm(3)/(6*ht3)-u_int(3)*tm(1)*tm(2)*tm(4)/(2*ht3)+...
                u_int(2)*tm(1)*tm(3)*tm(4)/(2*ht3)-u_int(1)*tm(2)*tm(3)*tm(4)/(6*ht3);
        end
        
        % ��������� ���������� ���������� �����, ��������� � ���������������
        f2 = feval(pdep1het,x(1)+(j-1)*hx,tcur,del_cur2);
        temp_u(j,i) = U(j,i-1) + (1-s)*k*(U(j+1,i-1)-2*U(j,i-1)+U(j-1,i-1))+ht*f2; 
    end

    % ���������� ������ �������� ��� ���������������� �������
    di = 2+1/(k*s); % ������������ ������� ������� �������
    d = temp_u(2:nx-1,i)/(k*s); % ������ ��������� ������
    d(1) = d(1) + temp_u(1,i);
    d(end) = d(end) + temp_u(nx,i);
    
    % ����������� �������� ����������� ������������� A, B
    A(1) = 1/di; B(1) = d(1)/di;
    for j=2:nx-2
        A(j) = 1/(di-A(j-1));
        B(j) = (d(j)+B(j-1))/(di-A(j-1));
    end

    % ���������� ������� � ����� � ���������� ��������� t
    U(nx-1,i) = B(end);
    for j=nx-2:-1:2
        U(j,i) = A(j-1)*U(j+1,i)+B(j-1);
    end

end

% ��������� ������ X, T
[T, X] = meshgrid(t, x);