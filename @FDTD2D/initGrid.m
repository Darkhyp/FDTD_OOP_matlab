%% initialize Grid for FDTD calculations
function obj = initGrid(obj, Objects, Courant)
%#codegen

global imp0 c

if ~isempty(Courant)
	obj.Courant = Courant;
end

if obj.isTE
	%grid coordinates
	[X,Y]			   = obj.getGridXYZ('Ez');
	[obj.xEz, obj.yEz]  = ndgrid(X,Y);
	[X,Y]			   = obj.getGridXYZ('Hx');
	[obj.xHx, obj.yHx]	= ndgrid(X,Y);
	[X,Y]			   = obj.getGridXYZ('Hy');
	[obj.xHy, obj.yHy]	= ndgrid(X,Y);
else
	%grid coordinates
	[X,Y]			   = obj.getGridXYZ('Ex');
	[obj.xEx, obj.yEx]	= ndgrid(X,Y);
	[X,Y]			   = obj.getGridXYZ('Ey');
	[obj.xEy, obj.yEy]	= ndgrid(X,Y);
	[X,Y]			   = obj.getGridXYZ('Hz');
	[obj.xHz, obj.yHz]	= ndgrid(X,Y);
end
[X,Y]		   = obj.getGridXYZ('');
[obj.x, obj.y]	= ndgrid(X,Y);
	

if ~isempty(Objects)
	for n_obj=1:length(Objects)
		otmp = Objects{n_obj};
		switch otmp.type
			case 'sphere'
				test_func = @(X,Y) ((X-otmp.position(1))/otmp.radius(1)).^2 ...
								 + ((Y-otmp.position(2))/otmp.radius(2)).^2 <= 1;

				if obj.isTE 
					% surrounding Ez nodes
					condition = test_func(obj.xEz(2:end-1,2:end-1),obj.yEz(2:end-1,2:end-1));
					% lossy dielectric

					obj.epsilon_zz(condition) = otmp.RI^2;
					obj.sigma_zz(condition)   = otmp.loss;
				else
					% surrounding Ex nodes
					condition = test_func(obj.xEx(:,2:end-1),obj.yEx(:,2:end-1));
					% lossy dielectric
					obj.epsilon_xx(condition) = otmp.RI^2;
					obj.sigma_xx(condition) = otmp.loss;

					% surrounding Ey nodes
					condition = test_func(obj.xEy(2:end-1,:),obj.yEy(2:end-1,:));
					% lossy dielectric
					obj.epsilon_yy(condition) = otmp.RI^2;
					obj.sigma_yy(condition)   = otmp.loss;
				end
		end
	end
end
					
if obj.isTE
	% set electric-field update coefficients
	obj.coefEz   = (1-obj.sigma_zz)./(1+obj.sigma_zz);

	obj.coefEz_y = obj.Courant(1)*imp0 ./obj.epsilon_zz./(1+obj.sigma_zz);
	obj.coefEz_x = obj.Courant(2)*imp0 ./obj.epsilon_zz./(1+obj.sigma_zz);
	% set magnetic-field update coefficients
	obj.coefHy_z = obj.Courant(1)/imp0 ./obj.mu_yy;
	obj.coefHx_z = obj.Courant(2)/imp0 ./obj.mu_xx;
else
	% set electric-field update coefficients
	obj.coefEx   = (1-obj.sigma_xx)./(1+obj.sigma_xx);
	obj.coefEx_z = obj.Courant(2)*imp0 ./obj.epsilon_xx./(1+obj.sigma_xx);
	
	obj.coefEy   = (1-obj.sigma_yy)./(1+obj.sigma_yy);
	obj.coefEy_z = obj.Courant(1)*imp0 ./obj.epsilon_yy./(1+obj.sigma_yy);
	% set magnetic-field update coefficients
	obj.coefHz_y = obj.Courant(1)/imp0 ./obj.mu_zz;
	obj.coefHz_x = obj.Courant(2)/imp0 ./obj.mu_zz;
end

	
% initializations for PML
switch obj.PMLtype
	case 1
		% PML = struct('type','split PML', 'layers',10, 'm',3.5, 'coef',0.8);

		if obj.isTE
			% Hx grid
			f = funcPML1(obj.yHx, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
			obj.PML.byHx = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(2)*f); % c*obj.dt/obj.dy

			% Hy grid
			f = funcPML1(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
			obj.PML.bxHy = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(1)*f); % c*obj.dt/obj.dx

			% Ez grid (2:end-1,2:end-1,:)
			f = funcPML1(obj.xEz(2:end-1,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
			obj.PML.bxEz = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx

			f = funcPML1(obj.yEz(2:end-1,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
			obj.PML.byEz = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(2)*f); % c*obj.dt/obj.dy
		else
			% Hz grid
			f = funcPML1(obj.xHz, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
			obj.PML.bxHz = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(1)*f); % c*obj.dt/obj.dx

			f = funcPML1(obj.yHz, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
			obj.PML.byHz = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(2)*f); % c*obj.dt/obj.dy

			% Ex grid (:,2:end-1,2:end-1)
			f = funcPML1(obj.yEx(:,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
			obj.PML.byEx = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(2)*f); % c*obj.dt/obj.dy

			% Ey grid (2:end-1,:,2:end-1)
			f = funcPML1(obj.xEy(2:end-1,:), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
			obj.PML.bxEy = exp(-obj.PML.coef*(obj.PML.m+1)*obj.Courant(1)*f); % c*obj.dt/obj.dx
		end
	case 2
		% PML = struct('type','UPML', 'layers',10, 'm',4, 'R_err',1e-16, 'ka_max',1);

		sigma_max_x = -log(obj.PML.R_err)*(obj.PML.m + 1)/(2*obj.PML.Wx) *c*obj.dt/2; %	imp0<-eta = sqrt(mu0/epsilon0*Material(1,1)/Material(1,2));
		sigma_max_y = -log(obj.PML.R_err)*(obj.PML.m + 1)/(2*obj.PML.Wy) *c*obj.dt/2;

		if obj.isTE
			% Ez grid (2:end-1,2:end-1,:)
			[fx,kx] = funcPML2(obj.xEz(2:end-1,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
			[fy,ky] = funcPML2(obj.yEz(2:end-1,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
			obj.PML.cFz_1 = (obj.epsilon_zz - obj.sigma_zz*obj.dt/2)./(obj.epsilon_zz + obj.sigma_zz*obj.dt/2);
			obj.PML.cFz_2 =									 imp0./(obj.epsilon_zz + obj.sigma_zz*obj.dt/2);
			obj.PML.cGz_1 = (kx - sigma_max_x*fx)./(kx + sigma_max_x*fx);
			obj.PML.cGz_2 =					 1./(kx + sigma_max_x*fx);
			obj.PML.cEz_1 = (ky - sigma_max_y*fy)./(ky + sigma_max_y*fy);
			obj.PML.cEz_2 =					 1./(ky + sigma_max_y*fy);

			% Hx grid
			[fx,kx] = funcPML2(obj.xHx, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
			[fy,ky] = funcPML2(obj.yHx, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
			obj.PML.cMx_1 = (ky - sigma_max_y*fy)./(ky + sigma_max_y*fy);
			obj.PML.cMx_2 = obj.Courant(2)/imp0./obj.mu_xx./(ky + sigma_max_y*fy);
			obj.PML.cHx_2 = (kx + sigma_max_x*fx);
			obj.PML.cHx_3 = (kx - sigma_max_x*fx);

			% Hy grid
			[fx,kx] = funcPML2(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
			[fy,ky] = funcPML2(obj.yHy, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
			obj.PML.cMy_2 =	 obj.Courant(1)/imp0./obj.mu_yy;
			obj.PML.cHy_1 = (kx - sigma_max_x*fx)./(kx + sigma_max_x*fx);
			obj.PML.cHy_2 = (ky + sigma_max_y*fy)./(kx + sigma_max_x*fx);
			obj.PML.cHy_3 = (ky - sigma_max_y*fy)./(kx + sigma_max_x*fx);
		else
			% Ex grid (:,2:end-1,2:end-1)
			[fx,kx] = funcPML2(obj.xEx(:,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
			[fy,ky] = funcPML2(obj.yEx(:,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
			obj.PML.cFx_1 = (obj.epsilon_xx - obj.sigma_xx*obj.dt/2)./(obj.epsilon_xx + obj.sigma_xx*obj.dt/2);
			obj.PML.cFx_2 =					  obj.Courant(2)*imp0./(obj.epsilon_xx + obj.sigma_xx*obj.dt/2);
			obj.PML.cGx_1 = (ky - sigma_max_y*fy)./(ky + sigma_max_y*fy);
			obj.PML.cGx_2 =					 1./(ky + sigma_max_y*fy);
			obj.PML.cEx_2 = (kx + sigma_max_x*fx);
			obj.PML.cEx_3 = (kx - sigma_max_x*fx);

			% Ey grid (2:end-1,:,2:end-1)
			[fx,kx] = funcPML2(obj.xEy(2:end-1,:), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
			[fy,ky] = funcPML2(obj.yEy(2:end-1,:), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
			obj.PML.cFy_1 = (obj.epsilon_yy - obj.sigma_yy*obj.dt/2)./(obj.epsilon_yy + obj.sigma_yy*obj.dt/2);
			obj.PML.cFy_2 =					  obj.Courant(1)*imp0./(obj.epsilon_yy + obj.sigma_yy*obj.dt/2);
			obj.PML.cEy_1 = (kx - sigma_max_x*fx)./(kx + sigma_max_x*fx);
			obj.PML.cEy_2 = (ky + sigma_max_y*fy)./(kx + sigma_max_x*fx);
			obj.PML.cEy_3 = (ky - sigma_max_y*fy)./(kx + sigma_max_x*fx);

			% Hz grid
			[fx,kx] = funcPML2(obj.xHz, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
			[fy,ky] = funcPML2(obj.yHz, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
			obj.PML.cMz_1 = (kx - sigma_max_x*fx)./(kx + sigma_max_x*fx);
			obj.PML.cMz_2 =	 1/imp0./obj.mu_zz./(kx + sigma_max_x*fx);
			obj.PML.cHz_1 = (ky - sigma_max_y*fy)./(ky + sigma_max_y*fy);
			obj.PML.cHz_2 =					 1./(ky + sigma_max_y*fy);
		end
end

end


%% Addition functions
function f = funcPML1(x,x1,x2,Wx,m)
	f = zeros(size(x));
	condition = x1>x;
	f(condition) = ((x1-x(condition))/Wx).^m;
	condition = x2<x;
	f(condition) = ((x(condition)-x2)/Wx).^m;
end

function [f,k] = funcPML2(x,x1,x2,Wx,m,ka_max)
	f = zeros(size(x));
	condition = x1>x;
	f(condition) = ((x1-x(condition))/Wx).^m;
	condition = x2<x;
	f(condition) = ((x(condition)-x2)/Wx).^m;
	k = 1 + (ka_max - 1)*f;
end
