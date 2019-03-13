function [fn, Z] = getfunction(N,val,Z, h)

%Governing equation
%Z(i-2) -4*Z(i-1) +6*Z(i) -4*Z(i+1) +Z(i+2) + V.^2.steps^4/(Z(i)+h).^2
fn(1,1) = Z(3) -1;
fn(2,1) = Z(1) -8*Z(2) + 8*Z(4) -Z(5);
index = 2;

for i=3:N-2
    index = index + 1;
    fn(index, 1) = Z(i-2) -4*Z(i-1) +6*Z(i) -4*Z(i+1) +Z(i+2) + val/(Z(i)+h)^2;
end

fn(N-1,1) = -Z(N-4) + 16*Z(N-3) - 30*Z(N-2) + 16*Z(N-1) -Z(N) ;
fn(N,1) = -Z(N-4) + 2*Z(N-3) -2*Z(N-1) + Z(N);


