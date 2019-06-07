
function ljtest(N)

% N = 50;
x = linspace(0, N, N+1);
[x, y] = meshgrid(x, x);
x = [x(:), y(:)]';


tic();
disp('Pretty: ');
[E2, dE2] = energy(x, @lj);
toc();
tic();
disp('Pretty Again: ');
[E2, dE2] = energy(x, @lj);
toc();


tic();
disp('C-Style: ');
[E1, dE1] = energy_inline(x);
toc();
tic();
disp('C-Style Again: ');
[E1, dE1] = energy_inline(x);
toc();


end


function [J, dJ] = lj(r)
    s = sum(r.^2, 1);
    J = s^(-6) - 2 * s^(-3);
    dJ = -12 * (s^(-7) - s^(-4)) * r;
end


function [E, dE] = energy_inline(x)

d = size(x,1);
E = 0.0;
dE = zeros(size(x));

for n = 1:(size(x, 2)-1)
    for k = (n+1):size(x,2)
      s = 0.0;
      for i = 1:d
        r(i) = x(i, k) - x(i,n);
        s = s+r(i)*r(i);
      end
      %J = s^(-6) - 2 * s^(-3);
      %dJ = 12 * (s^(-7) - s^(-4)) * r;
      t= 1./(s*s*s);
      J = t*t-2.*t;
      dJ = 12 * (t*t-t)*s*r;
      
      E = E + J;
      for i = 1:d
        dE(i, k) = dE(i,k) + dJ(i);
        dE(i, n) = dE(i,n) - dJ(i);
      end
    end
end

end


function [E, dE] = energy(x, phi)

E = 0.0;
dE = zeros(size(x));

for n = 1:(size(x, 2)-1)
    for k = (n+1):size(x,2)
        [J, dJ] = phi(x(:,k)-x(:,n));
        E = E + J;
        dE(:, k) = dE(:,k) + dJ;
        dE(:, n) = dE(:,n) - dJ;
    end
end

end







