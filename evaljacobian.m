function J = evaljacobian(Jac, Z, val, h)

J = double(Jac);
[row, col] = size(J);

for i = 3:row-2
    for j = 3:col-2
        if i==j
            value = 6 - (2*val)/(Z(i)+h)^3;
            J(i,j) = double(value);
        end
    end
end