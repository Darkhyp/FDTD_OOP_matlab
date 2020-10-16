%% initialize Absorbing boundary conditioms
function obj = initABC(obj,type)
    switch obj.Dimensionality
        case 1
            switch type
                case '1abc' % first order ABC
                    obj.ABC.type = 1;
                    obj.ABC.oldEz_x	= [0, 0];
                otherwise % second order ABC
                    obj.ABC.type = 2;
                    obj.ABC.oldEz_x = zeros(3,2,2);
            end                   
%                     obj.ABC.coef = (obj.CourantNumber-1)/(obj.CourantNumber+1);
            % calculate coefficient on left-end and right-end of grid
%                     obj.ABC.coef = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end)
            obj.ABC.coef = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end)
        case 2
            % allocate memory for ABC arrays (second order)
            switch type
                case '1abc' % first order ABC
                    obj.ABC.type = 1;
                    if obj.isTE
                        obj.ABC.oldEz_x	= zeros(2, obj.Ny+1); % (m(start:end),n(all)) <- along x-axis
                        obj.ABC.oldEz_y	= zeros(obj.Nx+1, 2); % (m(all),n(start:end)) <- along y-axis
                        % calculate coefficient on left-end and right-end of grid
                        %{
                        obj.ABC.coef(:,2) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % ((start:end),y) <- along y-axis
                        obj.ABC.coef(:,1) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % ((start:end),x) <- along x-axis
                        %}
                        obj.ABC.coef(:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % ((start:end),y) <- along y-axis
                        obj.ABC.coef(:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % ((start:end),x) <- along x-axis
                    else
                        obj.ABC.oldEy_x	= zeros(2, obj.Ny); % (m(start:end),n(all)) <- along x-axis
                        obj.ABC.oldEx_y	= zeros(obj.Nx, 2); % (m(all),n(start:end)) <- along y-axis
                        % calculate coefficient on left-end and right-end of grid
                        %{
                        obj.ABC.coef(:,2) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % ((start:end),y) <- along y-axis
                        obj.ABC.coef(:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % ((start:end),x) <- along x-axis
                        %}
                        obj.ABC.coef(:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % ((start:end),y) <- along y-axis
                        obj.ABC.coef(:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % ((start:end),x) <- along x-axis
                    end
                otherwise % second order ABC
                    obj.ABC.type = 2;
                    if obj.isTE
                        obj.ABC.oldEz_x	= zeros(3,obj.Ny+1,2,2); % (m(1,2,3),n(all),prev:next,start:end) <- along x-axis
                        obj.ABC.oldEz_y	= zeros(obj.Nx+1,3,2,2); % (m(all),n(1,2,3),prev:next,start:end) <- along x-axis
                        % calculate coefficient on left-end and right-end of grid
                        %{
                        obj.ABC.coef(:,:,2) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % (3,start:end,y) <- along y-axis
                        obj.ABC.coef(:,:,1) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end,x) <- along x-axis
                        %}
                        obj.ABC.coef(:,:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end,y) <- along y-axis
                        obj.ABC.coef(:,:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end,x) <- along x-axis
                    else
                        obj.ABC.oldEy_x	= zeros(3,obj.Ny,2,2); % (m(1,2,3),n(all),prev:next,start:end) <- along x-axis
                        obj.ABC.oldEx_y	= zeros(obj.Nx,3,2,2); % (m(all),n(1,2,3),prev:next,start:end) <- along x-axis
                        % calculate coefficient on left-end and right-end of grid
                        %{
                        obj.ABC.coef(:,:,2) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (3,start:end) <- along y-axis
                        obj.ABC.coef(:,:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (3,start:end) <- along x-axis
                        %}
                        obj.ABC.coef(:,:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end) <- along y-axis
                        obj.ABC.coef(:,:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end) <- along x-axis
                    end
            end
        case 3
            switch type
                case '1abc' % first order ABC
                    obj.ABC.type = 1;
                    % allocate memory for ABC arrays
                    obj.ABC.oldEy_x = zeros(2,obj.Ny, obj.Nz+1); % (m(start:end,n,p)
                    obj.ABC.oldEz_x = zeros(2,obj.Ny+1,   obj.Nz);

                    obj.ABC.oldEx_y = zeros(obj.Nx, 2,obj.Nz+1);
                    obj.ABC.oldEz_y = zeros(obj.Nx+1,   2,obj.Nz);

                    obj.ABC.oldEx_z = zeros(obj.Nx, obj.Ny+1,  2);
                    obj.ABC.oldEy_z = zeros(obj.Nx+1,   obj.Ny,2);
                    % calculate ABC coefficients (first order)
                    %{
                    obj.ABC.coef(:,6) = coefABC(obj,'coefEy_x','coefHx_y',obj.ABC.type); % (2,start:end,z=3) <- Ey along z-axis
                    obj.ABC.coef(:,5) = coefABC(obj,'coefEx_y','coefHy_x',obj.ABC.type); % (2,start:end,z=3) <- Ex along z-axis
                    obj.ABC.coef(:,4) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % (2,start:end,y=2) <- Ez along y-axis
                    obj.ABC.coef(:,3) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (2,start:end,y=2) <- Ex along y-axis
                    obj.ABC.coef(:,2) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (2,start:end,x=1) <- Ez along x-axis
                    obj.ABC.coef(:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (2,start:end,x=1) <- Ey along x-axis
                    %}
                    obj.ABC.coef(:,6) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (2,start:end,z=3) <- Ey along z-axis
                    obj.ABC.coef(:,5) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (2,start:end,z=3) <- Ex along z-axis
                    obj.ABC.coef(:,4) = coefABC2(obj.CourantNumsber(2),obj.ABC.type); % (2,start:end,y=2) <- Ez along y-axis
                    obj.ABC.coef(:,3) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (2,start:end,y=2) <- Ex along y-axis
                    obj.ABC.coef(:,2) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (2,start:end,x=1) <- Ez along x-axis
                    obj.ABC.coef(:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (2,start:end,x=1) <- Ey along x-axis
%                             obj.ABC.coef(2,1:2) = (obj.CourantNumber(1:2)-1)./(obj.CourantNumber(1:2)+1);
%                             obj.ABC.coef(1,1:2) = obj.ABC.coef(2,1:3); % (start:end),(x,y,z)
                otherwise % second order ABC
                    obj.ABC.type = 2;
                    % allocate memory for ABC arrays
                    obj.ABC.oldEy_x = zeros(3,obj.Ny,obj.Nz+1,  2,2); % (m,n,p,prev:next,start:end)
                    obj.ABC.oldEz_x = zeros(3,obj.Ny+1,  obj.Nz,2,2);

                    obj.ABC.oldEx_y = zeros(obj.Nx, 3,obj.Nz+1,  2,2);
                    obj.ABC.oldEz_y = zeros(obj.Nx+1,   3,obj.Nz,2,2);

                    obj.ABC.oldEx_z = zeros(obj.Nx, obj.Ny+1,  3,2,2);
                    obj.ABC.oldEy_z = zeros(obj.Nx+1,   obj.Ny,3,2,2);
                    % calculate ABC coefficients (second order)
                    %{
                    obj.ABC.coef(:,:,6) = coefABC(obj,'coefEy_x','coefHx_y',obj.ABC.type); % (3,start:end,z=3) <- Ey along z-axis
                    obj.ABC.coef(:,:,5) = coefABC(obj,'coefEx_y','coefHy_x',obj.ABC.type); % (3,start:end,z=3) <- Ex along z-axis
                    obj.ABC.coef(:,:,4) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % (3,start:end,y=2) <- Ez along y-axis
                    obj.ABC.coef(:,:,3) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (3,start:end,y=2) <- Ex along y-axis
                    obj.ABC.coef(:,:,2) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end,x=1) <- Ez along x-axis
                    obj.ABC.coef(:,:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (3,start:end,x=1) <- Ey along x-axis
                    %}

                    obj.ABC.coef(:,:,6) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (3,start:end,z=3) <- Ey along z-axis
                    obj.ABC.coef(:,:,5) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (3,start:end,z=3) <- Ex along z-axis
                    obj.ABC.coef(:,:,4) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end,y=2) <- Ez along y-axis
                    obj.ABC.coef(:,:,3) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end,y=2) <- Ex along y-axis
                    obj.ABC.coef(:,:,2) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end,x=1) <- Ez along x-axis
                    obj.ABC.coef(:,:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end,x=1) <- Ey along x-axis

%                             obj.ABC.coef(1:3,2,1:3) = [-(1./obj.CourantNumber-2+obj.CourantNumber);
%                                                          -2*(obj.CourantNumber-1./obj.CourantNumber);
%                                                           4*(obj.CourantNumber+1./obj.CourantNumber)]./(1./obj.CourantNumber([1 1 1],:)+2+obj.CourantNumber([1 1 1],:));
%                             obj.ABC.coef(:,1,1:3) = obj.ABC.coef(:,2,1:3); % m,(start:end),(x,y,z)
            end
    end
end
