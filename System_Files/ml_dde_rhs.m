%% function written in the form required by dde23 solver
function f = sys_rhs_dde23(t,xx,yy,par)
f(1,1) = 1/par(12)*(par(13)-par(3)*(xx(1)-par(9))-par(4)*xx(2,1)*(xx(1)-par(11))-par(2)*(1/2+1/2*tanh((xx(1)-par(5))/par(6)))*(xx(1)-par(10))+1/2*par(14)*(1+tanh((yy(1,1)-par(15))/par(16))));
f(2,1) = par(1)*cosh(1/2*(xx(1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2));
end
