%% initialize Grid for FDTD calculations
function obj = initGrid(obj, Objects, Courant)
%#codegen

global imp0 c eps0 h eV

if ~isempty(Courant)
	obj.Courant = Courant;
end

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

if ~isempty(Objects)
	for n_obj=1:length(Objects)
		otmp = Objects{n_obj};
		switch otmp.type
			case 'sphere'
				test_func = @(X,Y,Z) ((X-otmp.position(1))/otmp.radius(1)).^2 ...
					+ ((Y-otmp.position(2))/otmp.radius(2)).^2 ...
					+ ((Z-otmp.position(3))/otmp.radius(3)).^2 <= 1;

                % surrounding Ex nodes
                conditionx = test_func(obj.xEx(:,2:end-1,2:end-1),obj.yEx(:,2:end-1,2:end-1),obj.zEx(:,2:end-1,2:end-1));
                % surrounding Ey nodes
                conditiony = test_func(obj.xEy(2:end-1,:,2:end-1),obj.yEy(2:end-1,:,2:end-1),obj.zEy(2:end-1,:,2:end-1));
                % surrounding Ez nodes
                conditionz = test_func(obj.xEz(2:end-1,2:end-1,:),obj.yEz(2:end-1,2:end-1,:),obj.zEz(2:end-1,2:end-1,:));

                if isfield(otmp,'dispersive')
                    switch otmp.dispersive
                        case 'Au Johnson'
                            eps_Inf = 1;
                            % in eV
                            c_CCPR = [1.53+420i;
                                     0.387+0.0314i;
                                     7.246+0.1796i];

                            a_CCPR = [0.0235+0.0865i;
                                     0.233+2.52i;
                                     1.186+2.390i;];
                        case 'Ag Palik'
                            eps_Inf = 1;
                            % in eV
                            c_CCPR = [0.5987+4195i;
                                    -0.2211+0.268i;
                                    -4.24+732.4i;
                                     0.6391-0.07186i;
                                     1.806+4.563i;
                                     1.443-82.19i];
                            a_CCPR = [0.02502+0.008626i;
                                     0.2021+0.9407i;
                                    14.67+1.338i;
                                     0.2997+4.034i;
                                     1.896+4.808i;
                                     9.396+6.477i];
                    end
                    
                    a_CCPR = a_CCPR*(eV/(h/2/pi)) *obj.dt/2;
                    c_CCPR = c_CCPR*(eV/(h/2/pi)) *obj.dt;
                    
                    obj.CCPR.REx = 0;
                    obj.CCPR.REy = 0;
                    obj.CCPR.REz = 0;

                    tmp = (1-a_CCPR)./(1+a_CCPR);
                    obj.CCPR.alphaEx = bsxfun(@times, conditionx, permute(tmp, [4 3 2 1]));
                    obj.CCPR.alphaEy = bsxfun(@times, conditiony, permute(tmp, [4 3 2 1]));
                    obj.CCPR.alphaEz = bsxfun(@times, conditionz, permute(tmp, [4 3 2 1]));

                    tmp = c_CCPR./(1+a_CCPR);
                    obj.CCPR.betaEx = bsxfun(@times, conditionx, permute(tmp, [4 3 2 1]));
                    obj.CCPR.betaEy = bsxfun(@times, conditiony, permute(tmp, [4 3 2 1]));
                    obj.CCPR.betaEz = bsxfun(@times, conditionz, permute(tmp, [4 3 2 1]));

                    tmp = eps_Inf + sum(real(tmp));
                    obj.epsilon_xx(conditionx) = tmp;
                    obj.epsilon_yy(conditiony) = tmp;
                    obj.epsilon_zz(conditionz) = tmp;

                    obj.CCPR.gammaEx = 1+obj.CCPR.alphaEx;
                    obj.CCPR.gammaEy = 1+obj.CCPR.alphaEy;
                    obj.CCPR.gammaEz = 1+obj.CCPR.alphaEz;

                    obj.CCPR.Jx = zeros(size(obj.CCPR.alphaEx));
                    obj.CCPR.Jy = zeros(size(obj.CCPR.alphaEy));
                    obj.CCPR.Jz = zeros(size(obj.CCPR.alphaEz));

                else
                    % lossy dielectric
                    obj.epsilon_xx(conditionx) = otmp.RI^2;
                    obj.sigma_xx(conditionx)   = otmp.loss;

                    % lossy dielectric
                    obj.epsilon_yy(conditiony) = otmp.RI^2;
                    obj.sigma_yy(conditiony)   = otmp.loss;

                    % lossy dielectric
                    obj.epsilon_zz(conditionz) = otmp.RI^2;
                    obj.sigma_zz(conditionz)   = otmp.loss;
                end
        end
 	end
end

% set electric-field update coefficients
obj.coefEx   = (1-obj.sigma_xx)./(1+obj.sigma_xx);
obj.coefEx_y = obj.Courant(3)*imp0 ./obj.epsilon_xx./(1+obj.sigma_xx);
obj.coefEx_z = obj.Courant(2)*imp0 ./obj.epsilon_xx./(1+obj.sigma_xx);

obj.coefEy   = (1-obj.sigma_yy)./(1+obj.sigma_yy);
obj.coefEy_x = obj.Courant(3)*imp0 ./obj.epsilon_yy./(1+obj.sigma_yy);
obj.coefEy_z = obj.Courant(1)*imp0 ./obj.epsilon_yy./(1+obj.sigma_yy);
    
obj.coefEz   = (1-obj.sigma_zz)./(1+obj.sigma_zz);
obj.coefEz_y = obj.Courant(1)*imp0 ./obj.epsilon_zz./(1+obj.sigma_zz);
obj.coefEz_x = obj.Courant(2)*imp0 ./obj.epsilon_zz./(1+obj.sigma_zz);

% set magnetic-field update coefficients */
obj.coefHx_y = obj.Courant(3)/imp0 ./obj.mu_xx;
obj.coefHx_z = obj.Courant(2)/imp0 ./obj.mu_xx;

obj.coefHy_x = obj.Courant(3)/imp0 ./obj.mu_yy;
obj.coefHy_z = obj.Courant(1)/imp0 ./obj.mu_yy;

obj.coefHz_x = obj.Courant(2)/imp0 ./obj.mu_zz;
obj.coefHz_y = obj.Courant(1)/imp0 ./obj.mu_zz;

% initializations for PML
switch obj.PMLtype
    case 1
		% PML = struct('type','split PML', 'layers',10, 'm',3.5, 'coef',0.8);

		% Hx grid
		f = funcPML1(obj.yHx, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
		obj.PML.byHx = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(2)*f); % c*obj.dt/obj.dy

		f = funcPML1(obj.zHx, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m);
		obj.PML.bzHx = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(3)*f); % c*obj.dt/obj.dz

		% Hy grid
		f = funcPML1(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
		obj.PML.bxHy = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx

		f = funcPML1(obj.zHy, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m);
		obj.PML.bzHy = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(3)*f); % c*obj.dt/obj.dz

		% Hz grid
		f = funcPML1(obj.xHz, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
		obj.PML.bxHz = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx

		f = funcPML1(obj.yHz, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
		obj.PML.byHz = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(2)*f); % c*obj.dt/obj.dy

		% Ex grid (:,2:end-1,2:end-1)
		f = funcPML1(obj.yEx(:,2:end-1,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
		obj.PML.byEx = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(2)*f); % c*obj.dt/obj.dy

		f = funcPML1(obj.zEx(:,2:end-1,2:end-1), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m);
		obj.PML.bzEx = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(3)*f); % c*obj.dt/obj.dz

		% Ey grid (2:end-1,:,2:end-1)
		f = funcPML1(obj.xEy(2:end-1,:,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
		obj.PML.bxEy = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx

		f = funcPML1(obj.zEy(2:end-1,:,2:end-1), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m);
		obj.PML.bzEy = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(3)*f); % c*obj.dt/obj.dz

		% Ez grid (2:end-1,2:end-1,:)
		f = funcPML1(obj.xEz(2:end-1,2:end-1,:), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m);
		obj.PML.bxEz = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(1)*f); % c*obj.dt/obj.dx

		f = funcPML1(obj.yEz(2:end-1,2:end-1,:), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m);
		obj.PML.byEz = exp(-(obj.PML.m+1)*obj.PML.coef*obj.Courant(2)*f); % c*obj.dt/obj.dy
    case 2
%         test = 0;
        test = -1;
        gamma = 0*obj.PML.gamma/eps0*obj.dt;
        % PML = struct('type','UPML', 'layers',10, 'm',4, 'R_err',1e-16, 'ka_max',1);
		
		sigma_max_x = -log(obj.PML.R_err)*(obj.PML.m + 1)/(2*obj.PML.Wx) *c*obj.dt; %    imp0<-eta = sqrt(mu0/epsilon0*Material(1,1)/Material(1,2));
		sigma_max_y = -log(obj.PML.R_err)*(obj.PML.m + 1)/(2*obj.PML.Wy) *c*obj.dt;
		sigma_max_z = -log(obj.PML.R_err)*(obj.PML.m + 1)/(2*obj.PML.Wz) *c*obj.dt;

		% Ex grid (:,2:end-1,2:end-1)
		[fx,kx,f_kx,condPMLx] = funcPML2(obj.xEx(:,2:end-1,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
		[fy,ky,f_ky,condPMLy] = funcPML2(obj.yEx(:,2:end-1,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
		[fz,kz,f_kz,condPMLz] = funcPML2(obj.zEx(:,2:end-1,2:end-1), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max);
        obj.conditions.ExPML = condPMLx | condPMLy | condPMLz;
        switch test
            case -1
                obj.PML.cGx_1 = (ky - (sigma_max_y*fy+gamma*ky)/2)./(ky + (sigma_max_y*fy+gamma*ky)/2);
                obj.PML.cGx_2 =                      (1 + gamma/2)./(ky + (sigma_max_y*fy+gamma*ky)/2);
                obj.PML.cGx_3 =                      (1 - gamma/2)./(ky + (sigma_max_y*fy+gamma*ky)/2);
                obj.PML.cEx_1 = (kz - (sigma_max_z*fz+gamma*kz)/2)./(kz + (sigma_max_z*fz+gamma*kz)/2);
                obj.PML.cEx_2 = (kx + (sigma_max_x*fx+gamma*kx)/2)./(kz + (sigma_max_z*fz+gamma*kz)/2);
                obj.PML.cEx_3 = (kx - (sigma_max_x*fx+gamma*kx)/2)./(kz + (sigma_max_z*fz+gamma*kz)/2);
            case 0
                obj.PML.cFx_1 = (obj.epsilon_xx - obj.sigma_xx*obj.dt/2)./(obj.epsilon_xx + obj.sigma_xx*obj.dt/2);
                obj.PML.cFx_2 =                                     imp0./(obj.epsilon_xx + obj.sigma_xx*obj.dt/2);
%                 obj.PML.cFx_2 = (1 + obj.PML.cFx_1)*imp0/2;
                obj.PML.cGx_1 = (ky - sigma_max_y/2*fy)./(ky + sigma_max_y/2*fy);
                obj.PML.cGx_2 =                       1./(ky + sigma_max_y/2*fy);
                obj.PML.cEx_1 = (kz - sigma_max_z/2*fz)./(kz + sigma_max_z/2*fz);
                obj.PML.cEx_2 = (kx + sigma_max_x/2*fx)./(kz + sigma_max_z/2*fz);
                obj.PML.cEx_3 = (kx - sigma_max_x/2*fx)./(kz + sigma_max_z/2*fz);
            case 1
                obj.PML.cFx_1 = exp(-obj.sigma_xx./obj.epsilon_xx*obj.dt);
                obj.PML.cFx_2 = (1 + obj.PML.cFx_1)*imp0/2;
                obj.PML.cGx_1 = exp(-sigma_max_y.*f_ky);
                obj.PML.cGx_2 = (1 + obj.PML.cGx_1)/2;
                obj.PML.cEx_1 = exp(-sigma_max_z.*f_kz);
                obj.PML.cEx_2 = kx./kz.*(1 + exp(-sigma_max_z*f_kz))./(1 + exp(-sigma_max_x*f_kx));
                obj.PML.cEx_3 = obj.PML.cEx_2.*exp(-sigma_max_x*f_kx);
            case 2
                obj.PML.cFx_1 = exp(-obj.sigma_xx./obj.epsilon_xx*obj.dt);
                obj.PML.cFx_2 = (1 + obj.PML.cFx_1)*imp0/2;
                obj.PML.cFx_1 = obj.PML.cFx_1(obj.conditions.ExPML);
                obj.PML.cGx_1 = exp(-sigma_max_y.*f_ky(obj.conditions.ExPML));
                obj.PML.cGx_2 = (1 + obj.PML.cGx_1)/2;
                obj.PML.cEx_1 = exp(-sigma_max_z.*f_kz(obj.conditions.ExPML));
                obj.PML.cEx_2 = kx(obj.conditions.ExPML)./kz(obj.conditions.ExPML).*(1 + exp(-sigma_max_z*f_kz(obj.conditions.ExPML)))./(1 + exp(-sigma_max_x*f_kx(obj.conditions.ExPML)));
                obj.PML.cEx_3 = obj.PML.cEx_2.*exp(-sigma_max_x*f_kx(obj.conditions.ExPML));
        end
        
		% Ey grid (2:end-1,:,2:end-1)
		[fx,kx,f_kx,condPMLx] = funcPML2(obj.xEy(2:end-1,:,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
		[fy,ky,f_ky,condPMLy] = funcPML2(obj.yEy(2:end-1,:,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
		[fz,kz,f_kz,condPMLz] = funcPML2(obj.zEy(2:end-1,:,2:end-1), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max);
        obj.conditions.EyPML = condPMLx | condPMLy | condPMLz;
        switch test
            case -1
                obj.PML.cGy_1 = (kz - (sigma_max_z*fz+gamma*kz)/2)./(kz + (sigma_max_z*fz+gamma*kz)/2);
                obj.PML.cGy_2 =                      (1 + gamma/2)./(kz + (sigma_max_z*fz+gamma*kz)/2);
                obj.PML.cGy_3 =                      (1 - gamma/2)./(kz + (sigma_max_z*fz+gamma*kz)/2);
                obj.PML.cEy_1 = (kx - (sigma_max_x*fx+gamma*kx)/2)./(kx + (sigma_max_x*fx+gamma*kx)/2);
                obj.PML.cEy_2 = (ky + (sigma_max_y*fy+gamma*ky)/2)./(kx + (sigma_max_x*fx+gamma*kx)/2);
                obj.PML.cEy_3 = (ky - (sigma_max_y*fy+gamma*ky)/2)./(kx + (sigma_max_x*fx+gamma*kx)/2);
            case 0
                obj.PML.cFy_1 = (obj.epsilon_yy - obj.sigma_yy*obj.dt/2)./(obj.epsilon_yy + obj.sigma_yy*obj.dt/2);
                obj.PML.cFy_2 =                                     imp0./(obj.epsilon_yy + obj.sigma_yy*obj.dt/2);
%                 obj.PML.cFy_2 = (1 + obj.PML.cFy_1)*imp0/2;
                obj.PML.cGy_1 = (kz - sigma_max_z/2*fz)./(kz + sigma_max_z/2*fz);
                obj.PML.cGy_2 =                       1./(kz + sigma_max_z/2*fz);
                obj.PML.cEy_1 = (kx - sigma_max_x/2*fx)./(kx + sigma_max_x/2*fx);
                obj.PML.cEy_2 = (ky + sigma_max_y/2*fy)./(kx + sigma_max_x/2*fx);
                obj.PML.cEy_3 = (ky - sigma_max_y/2*fy)./(kx + sigma_max_x/2*fx);
            case 1
                obj.PML.cFy_1 = exp(-obj.sigma_yy./obj.epsilon_yy*obj.dt);
                obj.PML.cFy_2 = (1 + obj.PML.cFy_1)*imp0/2;
                obj.PML.cGy_1 = exp(-sigma_max_z*f_kz);
                obj.PML.cGy_2 = (1 + obj.PML.cGy_1)/2;
                obj.PML.cEy_1 = exp(-sigma_max_x*f_kx);
                obj.PML.cEy_2 = ky./kx.*(1 + exp(-sigma_max_x*f_kx))./(1 + exp(-sigma_max_y*f_ky));
                obj.PML.cEy_3 = obj.PML.cEy_2.*exp(-sigma_max_y*f_ky);
            case 2
                obj.PML.cFy_1 = exp(-obj.sigma_yy./obj.epsilon_yy*obj.dt);
                obj.PML.cFy_2 = (1 + obj.PML.cFy_1)*imp0/2;
                obj.PML.cFy_1 = obj.PML.cFy_1(obj.conditions.EyPML);
                obj.PML.cGy_1 = exp(-sigma_max_z*f_kz(obj.conditions.EyPML));
                obj.PML.cGy_2 = (1 + obj.PML.cGy_1)/2;
                obj.PML.cEy_1 = exp(-sigma_max_x*f_kx(obj.conditions.EyPML));
                obj.PML.cEy_2 = ky(obj.conditions.EyPML)./kx(obj.conditions.EyPML).*(1 + exp(-sigma_max_x*f_kx(obj.conditions.EyPML)))./(1 + exp(-sigma_max_y*f_ky(obj.conditions.EyPML)));
                obj.PML.cEy_3 = obj.PML.cEy_2.*exp(-sigma_max_y*f_ky(obj.conditions.EyPML));
        end
        
		% Ez grid (2:end-1,2:end-1,:)
		[fx,kx,f_kx,condPMLx] = funcPML2(obj.xEz(2:end-1,2:end-1,:), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
		[fy,ky,f_ky,condPMLy] = funcPML2(obj.yEz(2:end-1,2:end-1,:), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
		[fz,kz,f_kz,condPMLz] = funcPML2(obj.zEz(2:end-1,2:end-1,:), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max);
        obj.conditions.EzPML = condPMLx | condPMLy | condPMLz;
        switch test
            case -1
                obj.PML.cGz_1 = (kx - (sigma_max_x*fx+gamma*kx)/2)./(kx + (sigma_max_x*fx+gamma*kx)/2);
                obj.PML.cGz_2 =                      (1 + gamma/2)./(kx + (sigma_max_x*fx+gamma*kx)/2);
                obj.PML.cGz_3 =                      (1 - gamma/2)./(kx + (sigma_max_x*fx+gamma*kx)/2);
                obj.PML.cEz_1 = (ky - (sigma_max_y*fy+gamma*ky)/2)./(ky + (sigma_max_y*fy+gamma*ky)/2);
                obj.PML.cEz_2 = (kz + (sigma_max_z*fz+gamma*kz)/2)./(ky + (sigma_max_y*fy+gamma*ky)/2);
                obj.PML.cEz_3 = (kz - (sigma_max_z*fz+gamma*kz)/2)./(ky + (sigma_max_y*fy+gamma*ky)/2);
            case 0
                obj.PML.cFz_1 = (obj.epsilon_zz - obj.sigma_zz*obj.dt/2)./(obj.epsilon_zz + obj.sigma_zz*obj.dt/2);
                obj.PML.cFz_2 =                                     imp0./(obj.epsilon_zz + obj.sigma_zz*obj.dt/2);
%                 obj.PML.cFz_2 = (1 + obj.PML.cFz_1)*imp0/2;
                obj.PML.cGz_1 = (kx - sigma_max_x/2*fx)./(kx + sigma_max_x/2*fx);
                obj.PML.cGz_2 =                       1./(kx + sigma_max_x/2*fx);
                obj.PML.cEz_1 = (ky - sigma_max_y/2*fy)./(ky + sigma_max_y/2*fy);
                obj.PML.cEz_2 = (kz + sigma_max_z/2*fz)./(ky + sigma_max_y/2*fy);
                obj.PML.cEz_3 = (kz - sigma_max_z/2*fz)./(ky + sigma_max_y/2*fy);
            case 1
                obj.PML.cFz_1 = exp(-obj.sigma_zz./obj.epsilon_zz*obj.dt);
                obj.PML.cFz_2 = (1 + obj.PML.cFz_1)*imp0/2;
                obj.PML.cGz_1 = exp(-sigma_max_x*f_kx);
                obj.PML.cGz_2 = (1 + obj.PML.cGz_1)/2;
                obj.PML.cEz_1 = exp(-sigma_max_y*f_ky);
                obj.PML.cEz_2 = kz./ky.*(1 + exp(-sigma_max_y*f_ky))./(1 + exp(-sigma_max_z*f_kz));
                obj.PML.cEz_3 = obj.PML.cEz_2.*exp(-sigma_max_z*f_kz);
            case 2
                obj.PML.cFz_1 = exp(-obj.sigma_zz./obj.epsilon_zz*obj.dt);
                obj.PML.cFz_2 = (1 + obj.PML.cFz_1)*imp0/2;
                obj.PML.cFz_1 = obj.PML.cFz_1(obj.conditions.EzPML);
                obj.PML.cGz_1 = exp(-sigma_max_x*f_kx(obj.conditions.EzPML));
                obj.PML.cGz_2 = (1 + obj.PML.cGz_1)/2;
                obj.PML.cEz_1 = exp(-sigma_max_y*f_ky(obj.conditions.EzPML));
                obj.PML.cEz_2 = kz(obj.conditions.EzPML)./ky(obj.conditions.EzPML).*(1 + exp(-sigma_max_y*f_ky(obj.conditions.EzPML)))./(1 + exp(-sigma_max_z*f_kz(obj.conditions.EzPML)));
                obj.PML.cEz_3 = obj.PML.cEz_2.*exp(-sigma_max_z*f_kz(obj.conditions.EzPML));
        end
        
		% Hx grid
		[fx,kx,f_kx,condPMLx] = funcPML2(obj.xHx, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
		[fy,ky,f_ky,condPMLy] = funcPML2(obj.yHx, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
		[fz,kz,f_kz,condPMLz] = funcPML2(obj.zHx, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max);
        obj.conditions.HxPML = condPMLx | condPMLy | condPMLz;
        switch test
            case {-1,0}
                obj.PML.cMx_1 = (ky - sigma_max_y/2*fy)./(ky + sigma_max_y/2*fy);
                obj.PML.cMx_2 =       1/imp0./obj.mu_xx./(ky + sigma_max_y/2*fy);
%                 obj.PML.cMx_2 = (1 + obj.PML.cMx_1)/imp0/2./obj.mu_xx;
                obj.PML.cHx_1 = (kz - sigma_max_z/2*fz)./(kz + sigma_max_z/2*fz);
                obj.PML.cHx_2 = (kx + sigma_max_x/2*fx)./(kz + sigma_max_z/2*fz);
                obj.PML.cHx_3 = (kx - sigma_max_x/2*fx)./(kz + sigma_max_z/2*fz);
            case 1
                obj.PML.cMx_1 = exp(-sigma_max_y*f_ky);
                obj.PML.cMx_2 = (1 + obj.PML.cMx_1)/imp0/2./obj.mu_xx;
                obj.PML.cHx_1 = exp(-sigma_max_z*f_kz);
                obj.PML.cHx_2 = kx./kz.*(1 + exp(-sigma_max_z*f_kz))./(1 + exp(-sigma_max_x*f_kx));
                obj.PML.cHx_3 = obj.PML.cHx_2.*exp(-sigma_max_x*f_kx);
            case 2
                obj.PML.cMx_1 = exp(-sigma_max_y*f_ky);
                obj.PML.cMx_2 = (1 + obj.PML.cMx_1)/imp0/2./obj.mu_xx;
                obj.PML.cMx_1 = obj.PML.cMx_1(obj.conditions.HxPML);
                obj.PML.cHx_1 = exp(-sigma_max_z*f_kz(obj.conditions.HxPML));
                obj.PML.cHx_2 = kx(obj.conditions.HxPML)./kz(obj.conditions.HxPML).*(1 + exp(-sigma_max_z*f_kz(obj.conditions.HxPML)))./(1 + exp(-sigma_max_x*f_kx(obj.conditions.HxPML)));
                obj.PML.cHx_3 = obj.PML.cHx_2.*exp(-sigma_max_x*f_kx(obj.conditions.HxPML));
        end
        
		% Hy grid
		[fx,kx,f_kx,condPMLx] = funcPML2(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
		[fy,ky,f_ky,condPMLy] = funcPML2(obj.yHy, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
		[fz,kz,f_kz,condPMLz] = funcPML2(obj.zHy, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max);
        obj.conditions.HyPML = condPMLx | condPMLy | condPMLz;
        switch test
            case {-1,0}
                obj.PML.cMy_1 = (kz - sigma_max_z/2*fz)./(kz + sigma_max_z/2*fz);
                obj.PML.cMy_2 =       1/imp0./obj.mu_yy./(kz + sigma_max_z/2*fz);
%                 obj.PML.cMy_2 = (1 + obj.PML.cMy_1)/imp0/2./obj.mu_yy;
                obj.PML.cHy_1 = (kx - sigma_max_x/2*fx)./(kx + sigma_max_x/2*fx);
                obj.PML.cHy_2 = (ky + sigma_max_y/2*fy)./(kx + sigma_max_x/2*fx);
                obj.PML.cHy_3 = (ky - sigma_max_y/2*fy)./(kx + sigma_max_x/2*fx);
            case 1
                obj.PML.cMy_1 = exp(-sigma_max_z*f_kz);
                obj.PML.cMy_2 = (1 + obj.PML.cMy_1)/imp0/2./obj.mu_yy;
                obj.PML.cHy_1 = exp(-sigma_max_x*f_kx);
                obj.PML.cHy_2 = ky./kx.*(1 + exp(-sigma_max_x*f_kx))./(1 + exp(-sigma_max_y*f_ky));
                obj.PML.cHy_3 = obj.PML.cHy_2.*exp(-sigma_max_y*f_ky);
            case 2
                obj.PML.cMy_1 = exp(-sigma_max_z*f_kz);
                obj.PML.cMy_2 = (1 + obj.PML.cMy_1)/imp0/2./obj.mu_yy;
                obj.PML.cMy_1 = obj.PML.cMy_1(obj.conditions.HyPML);
                obj.PML.cHy_1 = exp(-sigma_max_x*f_kx(obj.conditions.HyPML));
                obj.PML.cHy_2 = ky(obj.conditions.HyPML)./kx(obj.conditions.HyPML).*(1 + exp(-sigma_max_x*f_kx(obj.conditions.HyPML)))./(1 + exp(-sigma_max_y*f_ky(obj.conditions.HyPML)));
                obj.PML.cHy_3 = obj.PML.cHy_2.*exp(-sigma_max_y*f_ky(obj.conditions.HyPML));
        end
        
		% Hz grid
		[fx,kx,f_kx,condPMLx] = funcPML2(obj.xHz, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max);
		[fy,ky,f_ky,condPMLy] = funcPML2(obj.yHz, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max);
		[fz,kz,f_kz,condPMLz] = funcPML2(obj.zHz, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max);
        obj.conditions.HzPML = condPMLx | condPMLy | condPMLz;
        switch test
            case {-1,0}
                obj.PML.cMz_1 = (kx - sigma_max_x/2*fx)./(kx + sigma_max_x/2*fx);
                obj.PML.cMz_2 =       1/imp0./obj.mu_zz./(kx + sigma_max_x/2*fx);
%                 obj.PML.cMz_2 = (1 + obj.PML.cMz_1)/imp0/2./obj.mu_zz;
                obj.PML.cHz_1 = (ky - sigma_max_y/2*fy)./(ky + sigma_max_y/2*fy);
                obj.PML.cHz_2 = (kz + sigma_max_z/2*fz)./(ky + sigma_max_y/2*fy);
                obj.PML.cHz_3 = (kz - sigma_max_z/2*fz)./(ky + sigma_max_y/2*fy);
            case 1
                obj.PML.cMz_1 = exp(-sigma_max_x*f_kx);
                obj.PML.cMz_2 = (1 + obj.PML.cMz_1)/imp0/2./obj.mu_zz;
                obj.PML.cHz_1 = exp(-sigma_max_y*f_ky);
                obj.PML.cHz_2 = kz./ky.*(1 + exp(-sigma_max_y*f_ky))./(1 + exp(-sigma_max_z*f_kz));
                obj.PML.cHz_3 = obj.PML.cHz_2.*exp(-sigma_max_z*f_kz);
            case 2
                obj.PML.cMz_1 = exp(-sigma_max_x*f_kx);
                obj.PML.cMz_2 = (1 + obj.PML.cMz_1)/imp0/2./obj.mu_zz;
                obj.PML.cMz_1 = obj.PML.cMz_1(obj.conditions.HzPML);
                obj.PML.cHz_1 = exp(-sigma_max_y*f_ky(obj.conditions.HzPML));
                obj.PML.cHz_2 = kz(obj.conditions.HzPML)./ky(obj.conditions.HzPML).*(1 + exp(-sigma_max_y*f_ky(obj.conditions.HzPML)))./(1 + exp(-sigma_max_z*f_kz(obj.conditions.HzPML)));
                obj.PML.cHz_3 = obj.PML.cHz_2.*exp(-sigma_max_z*f_kz(obj.conditions.HzPML));
        end
        
		%{
		%% General FDTD coefficients 
		obj.PML.K_a = (2*eps0*obj.epsilon - obj.sigma*obj.dt)./(2*eps0*obj.epsilon + obj.sigma*obj.dt);
		obj.PML.K_b =                                2*obj.dt./(2*eps0*obj.epsilon + obj.sigma*obj.dt);
		obj.PML.K_c = obj.mu;
   
		%% PML coefficients along x-axis
		sigma_max = -(obj.PML.m + 1)*log(obj.PML.R_err)/(2*eta*obj.PML.Wx*obj.dx);

		% Ex grid (:,2:end-1,2:end-1)
		[f,k] = funcPML2(obj.xEz(:,2:end-1,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Ex_c =   2*eps0*k + sigma_max*f*obj.dt;
        obj.PML.k_Ex_d = -(2*eps0*k - sigma_max*f*obj.dt);

		% Ey grid (2:end-1,:,2:end-1)
		[f,k] = funcPML2(obj.xEy(2:end-1,:,2:end-1), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Ey_a = (2*eps0*k - sigma_max*f*obj.dt)./(2*eps0*k + sigma_max*f*obj.dt);
        obj.PML.k_Ey_b =                               1./(2*eps0*k + sigma_max*f*obj.dt);

		% Ez grid (2:end-1,2:end-1,:)
		[f,kx_Ez] = funcPML2(obj.xEz(2:end-1,2:end-1,:), obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Gz_a = (2*eps0*k - sigma_max*f*obj.dt)./(2*eps0*k + sigma_max*f*obj.dt);
        obj.PML.k_Gz_b =                          2*eps0./(2*eps0*k + sigma_max*f*obj.dt);
        
		% Hx grid
		[f,k] = funcPML2(obj.xHx, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max)
		obj.PML.k_Hx_c =  (2*eps0*k + sigma_max*f*obj.dt)/mu0;
        obj.PML.k_Hx_d = -(2*eps0*k - sigma_max*f*obj.dt)/mu0;

		% Hy grid
		[f,k] = funcPML2(obj.xHy, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Hy_a = (2*eps0*k - sigma_max*f*obj.dt)./(2*eps0*k + sigma_max*f*obj.dt);
        obj.PML.k_Hy_b =                               1./(2*eps0*k + sigma_max*f*obj.dt);

		% Hz grid
		[f,k] = funcPML2(obj.xHz, obj.PML.x1,obj.PML.x2,obj.PML.Wx, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Bz_a = (2*eps0*k - sigma_max*f*obj.dt)./(2*eps0*k + sigma_max*f*obj.dt);
        obj.PML.k_Bz_b =                   2*eps0*obj.dt./(2*eps0*k + sigma_max*f*obj.dt);

		%% PML coefficients along y-axis
		sigma_max = -(obj.PML.m + 1)*log(obj.PML.R_err)/(2*eta*obj.PML.Wy*obj.dy);

		% Ex grid (:,2:end-1,2:end-1)
		[f,k] = funcPML2(obj.yEx(:,2:end-1,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max)
		obj.PML.k_Gx_a = (2*eps0*k - sigma_max*f*obj.dt)./(2*eps0*k + sigma_max*f*obj.dt);
        obj.PML.k_Gx_b =                          2*eps0./(2*eps0*k + sigma_max*f*obj.dt);

		% Ey grid (2:end-1,:,2:end-1)
		[f,k] = funcPML2(obj.yEy(2:end-1,:,2:end-1), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max)
		obj.PML.k_Ey_c =   2*eps0*k + sigma_max*f*obj.dt;
		obj.PML.k_Ey_d = -(2*eps0*k - sigma_max*f*obj.dt);

		% Ez grid (2:end-1,2:end-1,:)
		[f,k] = funcPML2(obj.yEz(2:end-1,2:end-1,:), obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Ez_a = (2*eps*k - sigma_max*f*obj.dt)./(2*eps*k + sigma_max*f*obj.dt);
        obj.PML.k_Ez_b =                              1./(2*eps*k + sigma_max*f*obj.dt);

		% Hx grid
		[f,k] = funcPML2(obj.yHx, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Bx_a = (2*eps*k - sigma_max*f*obj.dt)./(2*eps*k + sigma_max*f*obj.dt);
        obj.PML.k_Bx_b =                  2*eps0*obj.dt./(2*eps*k + sigma_max*f*obj.dt);

		% Hy grid
		[f,k] = funcPML2(obj.yHy, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Hy_c =  (2*eps*k + sigma_max*f*obj.dt)./mu0;
        obj.PML.k_Hy_d = -(2*eps*k - sigma_max*f*obj.dt)./mu0;

		% Hz grid
		[f,k] = funcPML2(obj.yHz, obj.PML.y1,obj.PML.y2,obj.PML.Wy, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Hz_a = (2*eps*k - sigma_max*f*obj.dt)./(2*eps*k + sigma_max*f*obj.dt);
        obj.PML.k_Hz_b =                              1./(2*eps*k + sigma_max*f*obj.dt);

		%% PML coefficients along z-axis 
		sigma_max = -(obj.PML.m + 1)*log(obj.PML.R_err)/(2*eta*obj.PML.Wz*obj.dz);

		% Ex grid (:,2:end-1,2:end-1)
		[f,k] = funcPML2(obj.zEx(:,2:end-1,2:end-1), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Ex_a = (2*eps*k - sigma_max*f*obj.dt)./(2*eps*k + sigma_max*f*obj.dt);
        obj.PML.k_Ex_b =                              1./(2*eps*k + sigma_max*f*obj.dt);

		% Ey grid (2:end-1,:,2:end-1)
		[f,k] = funcPML2(obj.zEy(2:end-1,:,2:end-1), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max)
		obj.PML.k_Gy_a = (2*eps0*k - sigma_max*f*obj.dt)./(2*eps0*k + sigma_max*f*obj.dt);
        obj.PML.k_Gy_b =                          2*eps0./(2*eps0*k + sigma_max*f*obj.dt);

		% Ez grid (2:end-1,2:end-1,:)
		[f,k] = funcPML2(obj.zEz(2:end-1,2:end-1,:), obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max)
		obj.PML.k_Ez_c =   2*eps0*k + sigma_max*f*obj.dt;
		obj.PML.k_Ez_d = -(2*eps0*k - sigma_max*f*obj.dt);
		
		% Hx grid
		[f,k] = funcPML2(obj.zHx, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Hx_a = (2*eps*k - sigma_max*f*obj.dt)./(2*eps*k + sigma_max*f*obj.dt);
        obj.PML.k_Hx_b =                              1./(2*eps*k + sigma_max*f*obj.dt);

		% Hy grid
		[f,k] = funcPML2(obj.zHy, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_By_a = (2*eps*k - sigma_max*f*obj.dt)./(2*eps*k + sigma_max*f*obj.dt);
        obj.PML.k_By_b =                  2*eps0*obj.dt./(2*eps*k + sigma_max*f*obj.dt);

		% Hz grid
		[f,k] = funcPML2(obj.zHz, obj.PML.z1,obj.PML.z2,obj.PML.Wz, obj.PML.m,obj.PML.ka_max)
        obj.PML.k_Hz_c =  (2*eps*k + sigma_max*f*obj.dt)./mu0;
        obj.PML.k_Hz_d = -(2*eps*k - sigma_max*f*obj.dt)./mu0;
		%}
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

function [f,k,f_k,condition] = funcPML2(x,x1,x2,Wx,m,ka_max)
    f = zeros(size(x));
    condition1 = x1>x;
    f(condition1) = ((x1-x(condition1))/Wx).^m;
    condition2 = x2<x;
    f(condition2) = ((x(condition2)-x2)/Wx).^m;
    k = 1 + (ka_max - 1)*f;

    % output
    condition = condition1 | condition2;
%     f = f(condition);
%     k = k(condition);
    f_k = f./k;
end
