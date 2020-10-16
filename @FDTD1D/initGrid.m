%% initialize Grid for FDTD calculations 1D
function obj = initGrid(obj, Objects, Courant)
%#codegen

global imp0 c
if ~isempty(Courant)
	obj.Courant = Courant;
end
	

%grid coordinates
obj.xEz = obj.getGridXYZ('Ez');
obj.xHy = obj.getGridXYZ('Hy');
obj.x   = obj.getGridXYZ('');


if ~isempty(Objects)
	for n_obj=1:length(Objects)
		otmp = Objects{n_obj};
		switch otmp.type
			case 'sphere' % line in 1D
				test_func = @(X) ((X-otmp.position(1))/otmp.radius(1)).^2 <= 1;

                % surrounding Ez nodes
                condition = test_func(obj.xEz(2:end-1));
                % lossy dielectric
                obj.epsilon_zz(condition) = otmp.RI^2;
                obj.sigma_zz(condition) = otmp.RI^2;
		end
	end
end

% set electric-field update coefficients
obj.coefEz   = (1-obj.sigma_zz)./(1+obj.sigma_zz);

obj.coefEz_y = obj.Courant(1)*imp0 ./obj.epsilon_zz./(1+obj.sigma_zz);
% set magnetic-field update coefficients
obj.coefHy_z = obj.Courant(1)/imp0 ./obj.mu_yy;

% initializations for PML
switch obj.PMLtype
	case 1
		% PML = struct('type','split PML', 'layers',10, 'm',3.5, 'coef',0.8);

		% Hy grid
		f = funcPML1(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
		obj.PML.bxHy = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx

		% Ez grid (2:end-1)
		f = funcPML1(obj.xEz(2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
		obj.PML.bxEz = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx
	case 2
        % PML = struct('type','UPML', 'layers',10, 'm',4, 'R_err',1e-16, 'ka_max',1);
		
		sigma_max_x = -log(obj.PML.R_err)*(obj.PML.m + 1)/(2*obj.PML.Wx) *c*obj.dt/2; %    imp0<-eta = sqrt(mu0/epsilon0*Material(1,1)/Material(1,2));

        % Ez grid (2:end-1,2:end-1,:)
		[fx,kx] = funcPML2(obj.xEz(2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
        obj.PML.cFz_1 = (obj.epsilon_zz - obj.sigma_zz*obj.dt/2)./(obj.epsilon_zz + obj.sigma_zz*obj.dt/2);
        obj.PML.cFz_2 =                                     imp0./(obj.epsilon_zz + obj.sigma_zz*obj.dt/2);
        obj.PML.cGz_1 = (kx - sigma_max_x*fx)./(kx + sigma_max_x*fx);
        obj.PML.cGz_2 =                     1./(kx + sigma_max_x*fx);
        
		% Hy grid
		[fx,kx] = funcPML2(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
        obj.PML.cHy_1 = (kx - sigma_max_x*fx)./(kx + sigma_max_x*fx);
        obj.PML.cHy_2 =     1/imp0./obj.mu_yy./(kx + sigma_max_x*fx);
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
