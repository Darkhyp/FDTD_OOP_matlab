function out = snapshotH(obj,type)
    if ~exist('type','var')
        type = '';
    end
    if strcmpi(type,'intensity')
        switch obj.Dimensionality
            case 1
                out = obj.Hy.^2;
            case 2
                if obj.isTE
                    out = obj.snapshotH('x0').^2 + obj.snapshotH('y0').^2;
                else
                    out = obj.Hz.^2;
                end
            case 3
                out = obj.snapshotH('x0').^2 + obj.snapshotH('y0').^2 + obj.snapshotH('z0').^2;
        end
        return
    end

    switch obj.Dimensionality
        case 1
            out = obj.Hy;
        case 2
            if obj.isTE
                switch type
                    % return as it
                    case 'x'
                        out = obj.Hx;
                    case 'y'
                        out = obj.Hy;
                    % fit to space grid
                    case 'x0'
                        out = (obj.Hx(1:end-1,:)+obj.Hx(2:end,:))/2;
                    case 'y0'
                        out = (obj.Hy(:,1:end-1)+obj.Hy(:,2:end))/2;
                    otherwise
                        out = obj.Hx;
                end
            else
                out = obj.Hz;
            end
        case 3
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
                otherwise
                    out = obj.Ez;
            end
    end
end
