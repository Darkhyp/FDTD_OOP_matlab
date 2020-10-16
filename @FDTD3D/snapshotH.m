function out = snapshotH(obj,type)
%#codegen

if nargin==1
	type = '';
end

switch type
	% return as it
	case 'x'
		out = obj.Hx;
	case 'y'
		out = obj.Hy;
	case 'z'
		out = obj.Hz;
	% fit to space grid
	case 'x0'
		out = (obj.Hx(1:end-1,:,:)+obj.Hx(2:end,:,:))/2;
	case 'y0'
		out = (obj.Hy(:,1:end-1,:)+obj.Hy(:,2:end,:))/2;
	case 'z0'
		out = (obj.Hz(:,:,1:end-1)+obj.Hz(:,:,2:end))/2;
	case 'intensity'
		out = obj.snapshotH('x0').^2 + obj.snapshotH('y0').^2 + obj.snapshotH('z0').^2;
	otherwise
		out = obj.Hz;
end

end
