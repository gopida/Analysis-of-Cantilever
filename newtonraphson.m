clear all;
clc;

% Grid Settings
gd = 40;
max_itr = 50;
N = gd - 4;
lhat = 1;
steps = lhat/(gd-1);
tolerance = 10^-60;


Z = sym('z', [1 gd]);

%voltage
volt = 300;


% Parameter values
w = 10^-6;
tb = 3*10^-6;
length = 250e-6;
td = 10^-3;
massp = 2329.6;
g = 3*10^-6;
E = 169e9;
e0 = 8.85e-12;
er = 11.68;
I = w*(tb^3)/12;

%Normalised values
Vhat = 2*volt*sqrt((e0*w*length^4)/(2*E*I*g^3));
h = td/(g*er);

val = Vhat^2*steps^4;

% Declate the finite difference function
[fn] = getfunction(gd, val, Z, h);
% Compute the jacobian function
jac = getjacobian(fn, Z, gd);
%initial value
init = zeros(size(Z));
for ii=1:gd
    init(ii) = g;
end

% Perform newton raphson
for i = 1:max_itr
    i
    J = subs(jac, Z, init);
    Fn = subs(fn, Z, init);
    Fn_db = double(Fn);
    J_db = evaljacobian(J, init, val, h)
    X = transpose(init) - (J_db\Fn_db);
    err(:,i) = abs(X-transpose(init));
    init = transpose(X);
    error(i) = sum(err(:,i));
    if(err(:,i) < tolerance)
        break;
    end
end

% Error vs itertion
figure, semilogy(1:i,error(1:i));
title('Newton Raphson');
xlabel('Iteration');
ylabel('Error');
grid minor;grid;

% plot x vs z
X = X*g;
figure, plot(1:N, X(3:gd-2));
title('Deflection of Cantilever beam');
xlabel('x');
ylabel('z');
grid;

