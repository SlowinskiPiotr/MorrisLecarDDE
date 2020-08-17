function J = sys_deri(xx,par,nx,np,v)

J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [1/par(12)*(-par(3)-par(4)*xx(2,1)-1/2*par(2)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))-par(2)*(1/2+1/2*tanh((xx(1,1)-par(5))/par(6)))), -1/par(12)*par(4)*(xx(1,1)-par(11)); 1/2*par(1)/par(8)*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))+1/2*par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))/par(8)*(1-tanh((xx(1,1)-par(7))/par(8))^2), -par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))];
	case 1
		J = [1/2/par(12)*par(14)/par(16)*(1-tanh((xx(1,2)-par(15))/par(16))^2), 0; 0, 0];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [0; cosh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))];
	case 2
		J = [-1/par(12)*(1/2+1/2*tanh((xx(1,1)-par(5))/par(6)))*(xx(1,1)-par(10)); 0];
	case 3
		J = [1/par(12)*(-xx(1,1)+par(9)); 0];
	case 4
		J = [-1/par(12)*xx(2,1)*(xx(1,1)-par(11)); 0];
	case 5
		J = [1/2/par(12)*par(2)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10)); 0];
	case 6
		J = [1/2/par(12)*par(2)*(xx(1,1)-par(5))/par(6)^2*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10)); 0];
	case 7
		J = [0; -1/2*par(1)/par(8)*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))-1/2*par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))/par(8)*(1-tanh((xx(1,1)-par(7))/par(8))^2)];
	case 8
		J = [0; -1/2*par(1)*(xx(1,1)-par(7))/par(8)^2*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))-1/2*par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))*(xx(1,1)-par(7))/par(8)^2*(1-tanh((xx(1,1)-par(7))/par(8))^2)];
	case 9
		J = [1/par(12)*par(3); 0];
	case 10
		J = [1/par(12)*par(2)*(1/2+1/2*tanh((xx(1,1)-par(5))/par(6))); 0];
	case 11
		J = [1/par(12)*par(4)*xx(2,1); 0];
	case 12
		J = [-1/par(12)^2*(par(13)-par(3)*(xx(1,1)-par(9))-par(4)*xx(2,1)*(xx(1,1)-par(11))-par(2)*(1/2+1/2*tanh((xx(1,1)-par(5))/par(6)))*(xx(1,1)-par(10))+1/2*par(14)*(1+tanh((xx(1,2)-par(15))/par(16)))); 0];
	case 13
		J = [1/par(12); 0];
	case 14
		J = [1/par(12)*(1/2+1/2*tanh((xx(1,2)-par(15))/par(16))); 0];
	case 15
		J = [-1/2/par(12)*par(14)/par(16)*(1-tanh((xx(1,2)-par(15))/par(16))^2); 0];
	case 16
		J = [-1/2/par(12)*par(14)*(xx(1,2)-par(15))/par(16)^2*(1-tanh((xx(1,2)-par(15))/par(16))^2); 0];
	case 17
		J = [0; 0];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [0, 0; 1/2/par(8)*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))+1/2*cosh(1/2*(xx(1,1)-par(7))/par(8))/par(8)*(1-tanh((xx(1,1)-par(7))/par(8))^2), -cosh(1/2*(xx(1,1)-par(7))/par(8))];
		case 2
			J = [-1/2/par(12)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))-1/par(12)*(1/2+1/2*tanh((xx(1,1)-par(5))/par(6))), 0; 0, 0];
		case 3
			J = [-1/par(12), 0; 0, 0];
		case 4
			J = [-1/par(12)*xx(2,1), -1/par(12)*(xx(1,1)-par(11)); 0, 0];
		case 5
			J = [-1/par(12)*par(2)/par(6)^2*tanh((xx(1,1)-par(5))/par(6))*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))+1/2/par(12)*par(2)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2), 0; 0, 0];
		case 6
			J = [1/2/par(12)*par(2)/par(6)^2*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))-1/par(12)*par(2)*(xx(1,1)-par(5))/par(6)^3*tanh((xx(1,1)-par(5))/par(6))*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))+1/2/par(12)*par(2)*(xx(1,1)-par(5))/par(6)^2*(1-tanh((xx(1,1)-par(5))/par(6))^2), 0; 0, 0];
		case 7
			J = [0, 0; -1/4*par(1)/par(8)^2*cosh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))-1/2*par(1)/par(8)^2*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1-tanh((xx(1,1)-par(7))/par(8))^2)+par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))/par(8)^2*tanh((xx(1,1)-par(7))/par(8))*(1-tanh((xx(1,1)-par(7))/par(8))^2), 1/2*par(1)/par(8)*sinh(1/2*(xx(1,1)-par(7))/par(8))];
		case 8
			J = [0, 0; -1/2*par(1)/par(8)^2*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))-1/4*par(1)*(xx(1,1)-par(7))/par(8)^3*cosh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))-1/2*par(1)*(xx(1,1)-par(7))/par(8)^3*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1-tanh((xx(1,1)-par(7))/par(8))^2)-1/2*par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))/par(8)^2*(1-tanh((xx(1,1)-par(7))/par(8))^2)+par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))*(xx(1,1)-par(7))/par(8)^3*tanh((xx(1,1)-par(7))/par(8))*(1-tanh((xx(1,1)-par(7))/par(8))^2), 1/2*par(1)*(xx(1,1)-par(7))/par(8)^2*sinh(1/2*(xx(1,1)-par(7))/par(8))];
		case 9
			J = [0, 0; 0, 0];
		case 10
			J = [1/2/par(12)*par(2)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2), 0; 0, 0];
		case 11
			J = [0, 1/par(12)*par(4); 0, 0];
		case 12
			J = [-1/par(12)^2*(-par(3)-par(4)*xx(2,1)-1/2*par(2)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))-par(2)*(1/2+1/2*tanh((xx(1,1)-par(5))/par(6)))), 1/par(12)^2*par(4)*(xx(1,1)-par(11)); 0, 0];
		case 13
			J = [0, 0; 0, 0];
		case 14
			J = [0, 0; 0, 0];
		case 15
			J = [0, 0; 0, 0];
		case 16
			J = [0, 0; 0, 0];
		case 17
			J = [0, 0; 0, 0];
		end
	case 1
		switch np
		case 1
			J = [0, 0; 0, 0];
		case 2
			J = [0, 0; 0, 0];
		case 3
			J = [0, 0; 0, 0];
		case 4
			J = [0, 0; 0, 0];
		case 5
			J = [0, 0; 0, 0];
		case 6
			J = [0, 0; 0, 0];
		case 7
			J = [0, 0; 0, 0];
		case 8
			J = [0, 0; 0, 0];
		case 9
			J = [0, 0; 0, 0];
		case 10
			J = [0, 0; 0, 0];
		case 11
			J = [0, 0; 0, 0];
		case 12
			J = [-1/2/par(12)^2*par(14)/par(16)*(1-tanh((xx(1,2)-par(15))/par(16))^2), 0; 0, 0];
		case 13
			J = [0, 0; 0, 0];
		case 14
			J = [1/2/par(12)/par(16)*(1-tanh((xx(1,2)-par(15))/par(16))^2), 0; 0, 0];
		case 15
			J = [1/par(12)*par(14)/par(16)^2*tanh((xx(1,2)-par(15))/par(16))*(1-tanh((xx(1,2)-par(15))/par(16))^2), 0; 0, 0];
		case 16
			J = [-1/2/par(12)*par(14)/par(16)^2*(1-tanh((xx(1,2)-par(15))/par(16))^2)+1/par(12)*par(14)*(xx(1,2)-par(15))/par(16)^3*tanh((xx(1,2)-par(15))/par(16))*(1-tanh((xx(1,2)-par(15))/par(16))^2), 0; 0, 0];
		case 17
			J = [0, 0; 0, 0];
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [1/par(12)*(par(2)/par(6)^2*tanh((xx(1,1)-par(5))/par(6))*(1-tanh((xx(1,1)-par(5))/par(6))^2)*(xx(1,1)-par(10))-par(2)/par(6)*(1-tanh((xx(1,1)-par(5))/par(6))^2))*v(1)-1/par(12)*par(4)*v(2), -1/par(12)*par(4)*v(1); (1/4*par(1)/par(8)^2*cosh(1/2*(xx(1,1)-par(7))/par(8))*(1/2+1/2*tanh((xx(1,1)-par(7))/par(8))-xx(2,1))+1/2*par(1)/par(8)^2*sinh(1/2*(xx(1,1)-par(7))/par(8))*(1-tanh((xx(1,1)-par(7))/par(8))^2)-par(1)*cosh(1/2*(xx(1,1)-par(7))/par(8))/par(8)^2*tanh((xx(1,1)-par(7))/par(8))*(1-tanh((xx(1,1)-par(7))/par(8))^2))*v(1)-1/2*par(1)/par(8)*sinh(1/2*(xx(1,1)-par(7))/par(8))*v(2), -1/2*par(1)/par(8)*sinh(1/2*(xx(1,1)-par(7))/par(8))*v(1)];
		case 1
			J = [0, 0; 0, 0];
		end
	case 1
		switch nx2
		case 0
			J = [0, 0; 0, 0];
		case 1
			J = [-1/par(12)*par(14)/par(16)^2*tanh((xx(1,2)-par(15))/par(16))*(1-tanh((xx(1,2)-par(15))/par(16))^2)*v(1), 0; 0, 0];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end