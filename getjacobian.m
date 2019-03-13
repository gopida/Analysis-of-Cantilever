function jac = getjacobian(fn, Z, N)

for m=1:N
    for n=1:N
        jac(m,n) = diff(fn(m,1),Z(n));
    end
end