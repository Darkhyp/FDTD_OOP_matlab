function out = snapshotE(obj,type)
%#codegen

if nargin==1
	type = '';
end

switch type
	case 'x'
		out = obj.Ex;
	case 'y'
		out = obj.Ey;
	case 'z'
		out = obj.Ez;
	% fit to space grid
	case 'x0'
		out = (obj.Ex(:,1:end-1,1:end-1)+obj.Ex(:,1:end-1,2:end)+obj.Ex(:,2:end,1:end-1)+obj.Ex(:,2:end,2:end))/4;
	case 'y0'
		out = (obj.Ey(1:end-1,:,1:end-1)+obj.Ey(1:end-1,:,2:end)+obj.Ey(2:end,:,1:end-1)+obj.Ey(2:end,:,2:end))/4;
	case 'z0'
		out = (obj.Ez(1:end-1,1:end-1,:)+obj.Ez(1:end-1,2:end,:)+obj.Ez(2:end,1:end-1,:)+obj.Ez(2:end,2:end,:))/4;
	case 'intensity'
		out = obj.snapshotE('x0').^2 + obj.snapshotE('y0').^2 + obj.snapshotE('z0').^2;
	otherwise
		% return as it
		out = obj.Ez;
end

end
