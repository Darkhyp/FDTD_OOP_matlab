function out = snapshotE(obj,type)
%#codegen

if nargin==1
	type = '';
end

if obj.isTE
    if strcmpi(type,'intensity')
        out = obj.Ez.^2;
    else
        out = obj.Ez;
    end
else
    switch type
        % return as it
        case 'x'
            out = obj.Ex;
        case 'y'
            out = obj.Ey;
        % fit to space grid
        case 'x0'
            out = (obj.Ex(:,1:end-1)+obj.Ex(:,2:end))/2;
        case 'y0'
            out = (obj.Ey(1:end-1,:)+obj.Ey(2:end,:))/2;
        case 'intensity'
            out = obj.snapshotE('x0').^2 + obj.snapshotE('y0').^2;
        otherwise
            out = obj.Ex;
    end
end
    
end
