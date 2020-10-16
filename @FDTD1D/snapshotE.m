function out = snapshotE(obj,type)
%#codegen

if nargin==1
	type = '';
end

switch type
	case 'z0'
		out = (obj.Ez(1:end-1)+obj.Ez(2:end))/2;
	case 'intensity'
		out = obj.snapshotE('z0').^2;
    otherwise
    	% return as it
		out = obj.Ez;
end

end
