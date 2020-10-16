%% initialize Grid for FDTD calculations
function obj = initGrid(obj, Objects, CourantNumber)
	global imp0
% 	global c
	if ~isempty(CourantNumber)
		obj.CourantNumber = CourantNumber;
	end
	
	switch obj.Dimensionality
		case 1 % 1D system
            % set electric-field update coefficients
            obj.coefEz_y(:) = obj.CourantNumber(1) * imp0;
            % set magnetic-field update coefficients */
            obj.coefHy_z(:) = obj.CourantNumber(1) / imp0;

            %grid coordinates
            obj.xEz = obj.getGridXYZ('Ez');
            obj.xHy = obj.getGridXYZ('Hy');
		case 2 % 2D system
			if obj.isTE
				% set electric-field update coefficients
				obj.coefEz_y(:,:) = obj.CourantNumber(1) * imp0;
				obj.coefEz_x(:,:) = obj.CourantNumber(2) * imp0;
				% set magnetic-field update coefficients */
				obj.coefHy_z(:,:) = obj.CourantNumber(1) / imp0;
				obj.coefHx_z(:,:) = obj.CourantNumber(2) / imp0;
				
				%grid coordinates
                [X,Y]               = obj.getGridXYZ('Ez');
				[obj.xEz, obj.yEz]  = ndgrid(X,Y);
                [X,Y]               = obj.getGridXYZ('Hx');
				[obj.xHx, obj.yHx]	= ndgrid(X,Y);
                [X,Y]               = obj.getGridXYZ('Hy');
				[obj.xHy, obj.yHy]	= ndgrid(X,Y);
			else
				% set electric-field update coefficients
				obj.coefEy_z(:,:) = obj.CourantNumber(1) * imp0;
				obj.coefEx_z(:,:) = obj.CourantNumber(2) * imp0;
				% set magnetic-field update coefficients */
				obj.coefHz_y(:,:) = obj.CourantNumber(1) / imp0;
				obj.coefHz_x(:,:) = obj.CourantNumber(2) / imp0;
				
				%grid coordinates
                [X,Y]               = obj.getGridXYZ('Ex');
				[obj.xEx, obj.yEx]	= ndgrid(X,Y);
                [X,Y]               = obj.getGridXYZ('Ey');
				[obj.xEy, obj.yEy]	= ndgrid(X,Y);
                [X,Y]               = obj.getGridXYZ('Hz');
				[obj.xHz, obj.yHz]	= ndgrid(X,Y);
			end
            [X,Y]           = obj.getGridXYZ('');
            [obj.x, obj.y]	= ndgrid(X,Y);
	case 3 % 3D system
		% set electric-field update coefficients
		obj.coefEx_y(:,:,:) = obj.CourantNumber(3) * imp0;
		obj.coefEx_z(:,:,:) = obj.CourantNumber(2) * imp0;
		obj.coefEy_x(:,:,:) = obj.CourantNumber(3) * imp0;
		obj.coefEy_z(:,:,:) = obj.CourantNumber(1) * imp0;
		obj.coefEz_x(:,:,:) = obj.CourantNumber(2) * imp0;
		obj.coefEz_y(:,:,:) = obj.CourantNumber(1) * imp0;
		
		% set magnetic-field update coefficients */
		obj.coefHx_y(:,:,:) = obj.CourantNumber(3) / imp0;
		obj.coefHx_z(:,:,:) = obj.CourantNumber(2) / imp0;
		obj.coefHy_x(:,:,:) = obj.CourantNumber(3) / imp0;
		obj.coefHy_z(:,:,:) = obj.CourantNumber(1) / imp0;
		obj.coefHz_x(:,:,:) = obj.CourantNumber(2) / imp0;
		obj.coefHz_y(:,:,:) = obj.CourantNumber(1) / imp0;

		%grid coordinates
		[X,Y,Z]   = obj.getGridXYZ('');
		[obj.x,   obj.y,   obj.z] = ndgrid(X,Y,Z);
		
		[X,Y,Z]   = obj.getGridXYZ('Ex');
		[obj.xEx, obj.yEx, obj.zEx] = ndgrid(X,Y,Z);
		[X,Y,Z]   = obj.getGridXYZ('Ey');
		[obj.xEy, obj.yEy, obj.zEy] = ndgrid(X,Y,Z);
		[X,Y,Z]   = obj.getGridXYZ('Ez');
		[obj.xEz, obj.yEz, obj.zEz] = ndgrid(X,Y,Z);
		[X,Y,Z]   = obj.getGridXYZ('Hx');
		[obj.xHx, obj.yHx, obj.zHx] = ndgrid(X,Y,Z);
		[X,Y,Z]   = obj.getGridXYZ('Hy');
		[obj.xHy, obj.yHy, obj.zHy] = ndgrid(X,Y,Z);
		[X,Y,Z]   = obj.getGridXYZ('Hz');
		[obj.xHz, obj.yHz, obj.zHz] = ndgrid(X,Y,Z);
	end
	
	if ~isempty(Objects)
		for n_obj=1:length(Objects)
			otmp = Objects{n_obj};
			switch otmp.type
				case 'sphere'
					switch obj.Dimensionality
						case 2
							test_func = @(X,Y) ((X-otmp.position(1))/otmp.radius(1)).^2 ...
								             + ((Y-otmp.position(2))/otmp.radius(2)).^2 <= 1;

							if obj.isTE 
								% surrounding Ez nodes
                                condition = test_func(obj.xEz(2:end-1,2:end-1),obj.yEz(2:end-1,2:end-1));
								% lossy dielectric
								obj.coefEz(condition)   = obj.coefEz(condition)*(1-otmp.loss)/(1+otmp.loss);
								obj.coefEz_y(condition) = obj.coefEz_y(condition)/(otmp.RI)^2/(1+otmp.loss);
								if obj.isPML
									obj.epsilon_z(condition) = otmp.RI^2;
								end
							else
								% surrounding Ex nodes
                                condition = test_func(obj.xEx(:,2:end-1),obj.yEx(:,2:end-1));
								% lossy dielectric
								obj.coefEx(condition)   = obj.coefEx(condition)*(1-otmp.loss)/(1+otmp.loss);
								obj.coefEx_y(condition) = obj.coefEx_y(condition)/(otmp.RI)^2/(1+otmp.loss);
								if obj.isPML
									obj.epsilon_x(condition) = otmp.RI^2;
								end

								% surrounding Ey nodes
                                condition = test_func(obj.xEy(2:end-1,:),obj.yEy(2:end-1,:));
								% lossy dielectric
								obj.coefEy(condition)   = obj.coefEy(condition)*(1-otmp.loss)/(1+otmp.loss);
								obj.coefEy_x(condition) = obj.coefEy_x(condition)/(otmp.RI)^2/(1+otmp.loss);
								if obj.isPML
									obj.epsilon_y(condition) = otmp.RI^2;
								end
							end
						case 3
							test_func = @(X,Y,Z) ((X-otmp.position(1))/otmp.radius(1)).^2 ...
								+ ((Y-otmp.position(2))/otmp.radius(2)).^2 ...
								+ ((Z-otmp.position(3))/otmp.radius(3)).^2 <= 1;

                            % surrounding Ex nodes
                            condition = test_func(obj.xEx(:,2:end-1,2:end-1),obj.yEx(:,2:end-1,2:end-1),obj.zEx(:,2:end-1,2:end-1));
							% lossy dielectric
							obj.coefEx(condition)   = obj.coefEx(condition)*(1-otmp.loss)/(1+otmp.loss);
							obj.coefEx_y(condition) = obj.coefEx_y(condition)/(otmp.RI)^2/(1+otmp.loss);
							obj.coefEx_z(condition) = obj.coefEx_z(condition)/(otmp.RI)^2/(1+otmp.loss);
							if obj.isPML
								obj.epsilon_x(condition) = otmp.RI^2;
							end

							% surrounding Ey nodes
							condition = test_func(obj.xEy(2:end-1,:,2:end-1),obj.yEy(2:end-1,:,2:end-1),obj.zEy(2:end-1,:,2:end-1));
							% lossy dielectric
							obj.coefEy(condition)   = obj.coefEy(condition)*(1-otmp.loss)/(1+otmp.loss);
							obj.coefEy_x(condition) = obj.coefEy_x(condition)/(otmp.RI)^2/(1+otmp.loss);
							obj.coefEy_z(condition) = obj.coefEy_z(condition)/(otmp.RI)^2/(1+otmp.loss);
							if obj.isPML
								obj.epsilon_y(condition) = otmp.RI^2;
							end

							% surrounding Ez nodes
                            condition = test_func(obj.xEz(2:end-1,2:end-1,:),obj.yEz(2:end-1,2:end-1,:),obj.zEz(2:end-1,2:end-1,:));
							% lossy dielectric
							obj.coefEz(condition)   = obj.coefEz(condition)*(1-otmp.loss)/(1+otmp.loss);
							obj.coefEz_x(condition) = obj.coefEz_x(condition)/(otmp.RI)^2/(1+otmp.loss);
							obj.coefEz_y(condition) = obj.coefEz_y(condition)/(otmp.RI)^2/(1+otmp.loss);
							if obj.isPML
								obj.epsilon_z(condition) = otmp.RI^2;
							end
					end
			end
		end
	end
	
	% initializations for PML
	if obj.isPML
		switch obj.Dimensionality
			case 1
				% Hy grid
				sigmaxHy = zeros(size(obj.Hy));
				condition = obj.PML.x1>obj.xHy;
				sigmaxHy(condition) = ((obj.PML.x1-obj.xHy(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
				condition = obj.PML.x2<obj.xHy;
				sigmaxHy(condition) = ((obj.xHy(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
				obj.PML.bxHy = exp(-sigmaxHy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

				% Ez grid (2:end-1)
				sigmaxEz = zeros(obj.Nx-1,1);
				x_tmp = obj.xEz(2:end-1); condition = obj.PML.x1>x_tmp;
				sigmaxEz(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
				condition = obj.PML.x2<x_tmp;
				sigmaxEz(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
				obj.PML.bxEz = exp(-sigmaxEz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx
			case 2
				if obj.isTE
					% Hx grid
					sigmayHx = zeros(size(obj.Hx));
					condition = obj.PML.y1>obj.yHx;
					sigmayHx(condition) = ((obj.PML.y1-obj.yHx(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
					condition = obj.PML.y2<obj.yHx;
					sigmayHx(condition) = ((obj.yHx(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
					obj.PML.byHx = exp(-sigmayHx*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy

					% Hy grid
					sigmaxHy = zeros(size(obj.Hy));
					condition = obj.PML.x1>obj.xHy;
					sigmaxHy(condition) = ((obj.PML.x1-obj.xHy(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
					condition = obj.PML.x2<obj.xHy;
					sigmaxHy(condition) = ((obj.xHy(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
					obj.PML.bxHy = exp(-sigmaxHy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1));  % c*obj.dt/obj.dx

					% Ez grid (2:end-1,2:end-1)
					sigmaxEz = zeros(obj.Nx-1,obj.Ny-1);
					x_tmp = obj.xEz(2:end-1,2:end-1); condition = obj.PML.x1>x_tmp;
					sigmaxEz(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
					condition = obj.PML.x2<x_tmp;
					sigmaxEz(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
					obj.PML.bxEz = exp(-sigmaxEz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

					sigmayEz = zeros(obj.Nx-1,obj.Ny-1);
					y_tmp = obj.yEz(2:end-1,2:end-1); condition = obj.PML.y1>y_tmp;
					sigmayEz(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
					condition = obj.PML.y2<y_tmp;
					sigmayEz(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
					obj.PML.byEz = exp(-sigmayEz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy
				else
					% Hz grid
					sigmaxHz = zeros(size(obj.Hz));
					condition = obj.PML.x1>obj.xHz;
					sigmaxHz(condition) = ((obj.PML.x1-obj.xHz(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
					condition = obj.PML.x2<obj.xHz;
					sigmaxHz(condition) = ((obj.xHz(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
					obj.PML.bxHz = exp(-sigmaxHz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

					sigmayHz = zeros(size(obj.Hz));
					condition = obj.PML.y1>obj.yHz;
					sigmayHz(condition) = ((obj.PML.y1-obj.yHz(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
					condition = obj.PML.y2<obj.yHz;
					sigmayHz(condition) = ((obj.yHz(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
					obj.PML.byHz = exp(-sigmayHz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy

					% Ex grid (:,2:end-1)
					sigmayEx = zeros(obj.Nx,obj.Ny-1);
					y_tmp = obj.yEx(:,2:end-1); condition = obj.PML.y1>y_tmp;
					sigmayEx(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
					condition = obj.PML.y2<y_tmp;
					sigmayEx(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
					obj.PML.byEx = exp(-sigmayEx*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy
					
					% Ey grid (2:end-1,:)
					sigmaxEy = zeros(obj.Nx-1,obj.Ny);
					x_tmp = obj.xEy(2:end-1,:); condition = obj.PML.x1>x_tmp;
					sigmaxEy(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
					condition = obj.PML.x2<x_tmp;
					sigmaxEy(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
					obj.PML.bxEy = exp(-sigmaxEy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx
				end
			case 3
				% Hx grid
				sigmayHx = zeros(size(obj.Hx));
				condition = obj.PML.y1>obj.yHx;
				sigmayHx(condition) = ((obj.PML.y1-obj.yHx(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
				condition = obj.PML.y2<obj.yHx;
				sigmayHx(condition) = ((obj.yHx(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
				obj.PML.byHx = exp(-sigmayHx*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy

				sigmazHx = zeros(size(obj.Hx));
				condition = obj.PML.z1>obj.zHx;
				sigmazHx(condition) = ((obj.PML.z1-obj.zHx(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
				condition = obj.PML.z2<obj.zHx;
				sigmazHx(condition) = ((obj.zHx(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
				obj.PML.bzHx = exp(-sigmazHx*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(3)); % c*obj.dt/obj.dz

				% Hy grid
				sigmaxHy = zeros(size(obj.Hy));
				condition = obj.PML.x1>obj.xHy;
				sigmaxHy(condition) = ((obj.PML.x1-obj.xHy(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
				condition = obj.PML.x2<obj.xHy;
				sigmaxHy(condition) = ((obj.xHy(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
				obj.PML.bxHy = exp(-sigmaxHy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

				sigmazHy = zeros(size(obj.Hy));
				condition = obj.PML.z1>obj.zHy;
				sigmazHy(condition) = ((obj.PML.z1-obj.zHy(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
				condition = obj.PML.z2<obj.zHy;
				sigmazHy(condition) = ((obj.zHy(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
				obj.PML.bzHy = exp(-sigmazHy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(3)); % c*obj.dt/obj.dz

				% Hz grid
				sigmaxHz = zeros(size(obj.Hz));
				condition = obj.PML.x1>obj.xHz;
				sigmaxHz(condition) = ((obj.PML.x1-obj.xHz(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
				condition = obj.PML.x2<obj.xHz;
				sigmaxHz(condition) = ((obj.xHz(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
				obj.PML.bxHz = exp(-sigmaxHz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

				sigmayHz = zeros(size(obj.Hz));
				condition = obj.PML.y1>obj.yHz;
				sigmayHz(condition) = ((obj.PML.y1-obj.yHz(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
				condition = obj.PML.y2<obj.yHz;
				sigmayHz(condition) = ((obj.yHz(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
				obj.PML.byHz = exp(-sigmayHz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy

				% Ex grid (:,2:end-1,2:end-1)
				sigmayEx = zeros(obj.Nx,obj.Ny-1,obj.Nz-1);
				y_tmp = obj.yEx(:,2:end-1,2:end-1); condition = obj.PML.y1>y_tmp;
				sigmayEx(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
				condition = obj.PML.y2<y_tmp;
				sigmayEx(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
				obj.PML.byEx = exp(-sigmayEx*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy

				sigmazEx = zeros(obj.Nx,obj.Ny-1,obj.Nz-1);
				z_tmp = obj.zEx(:,2:end-1,2:end-1); condition = obj.PML.z1>z_tmp;
				sigmazEx(condition) = ((obj.PML.z1-z_tmp(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
				condition = obj.PML.z2<z_tmp;
				sigmazEx(condition) = ((z_tmp(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
				obj.PML.bzEx = exp(-sigmazEx*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(3)); % c*obj.dt/obj.dz

				% Ey grid (2:end-1,:,2:end-1)
				sigmaxEy = zeros(obj.Nx-1,obj.Ny,obj.Nz-1);
				x_tmp = obj.xEy(2:end-1,:,2:end-1); condition = obj.PML.x1>x_tmp;
				sigmaxEy(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
				condition = obj.PML.x2<x_tmp;
				sigmaxEy(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
				obj.PML.bxEy = exp(-sigmaxEy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

				sigmazEy = zeros(obj.Nx-1,obj.Ny,obj.Nz-1);
				z_tmp = obj.xEy(2:end-1,:,2:end-1); condition = obj.PML.z1>z_tmp;
				sigmazEy(condition) = ((obj.PML.z1-z_tmp(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
				condition = obj.PML.z2<z_tmp;
				sigmazEy(condition) = ((z_tmp(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
				obj.PML.bzEy = exp(-sigmazEy*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(3)); % c*obj.dt/obj.dz

				% Ez grid (2:end-1,2:end-1,:)
				sigmaxEz = zeros(obj.Nx-1,obj.Ny-1,obj.Nz);
				x_tmp = obj.xEz(2:end-1,2:end-1,:); condition = obj.PML.x1>x_tmp;
				sigmaxEz(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
				condition = obj.PML.x2<x_tmp;
				sigmaxEz(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
				obj.PML.bxEz = exp(-sigmaxEz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(1)); % c*obj.dt/obj.dx

				sigmayEz = zeros(obj.Nx-1,obj.Ny-1,obj.Nz);
				y_tmp = obj.yEz(2:end-1,2:end-1,:); condition = obj.PML.y1>y_tmp;
				sigmayEz(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
				condition = obj.PML.y2<y_tmp;
				sigmayEz(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
				obj.PML.byEz = exp(-sigmayEz*(obj.PML.m+1)*obj.PML.coef*obj.CourantNumber(2)); % c*obj.dt/obj.dy
		end
	end

end
