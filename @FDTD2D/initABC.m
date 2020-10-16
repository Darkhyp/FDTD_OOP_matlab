%% initialize Absorbing boundary conditioms
function obj = initABC(obj,type)
%#codegen

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
			obj.ABC.coef(:,2) = coefABC2(obj.Courant(2),obj.ABC.type); % ((start:end),y) <- along y-axis
			obj.ABC.coef(:,1) = coefABC2(obj.Courant(1),obj.ABC.type); % ((start:end),x) <- along x-axis
		else
			obj.ABC.oldEy_x	= zeros(2, obj.Ny); % (m(start:end),n(all)) <- along x-axis
			obj.ABC.oldEx_y	= zeros(obj.Nx, 2); % (m(all),n(start:end)) <- along y-axis
			% calculate coefficient on left-end and right-end of grid
			%{
			obj.ABC.coef(:,2) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % ((start:end),y) <- along y-axis
			obj.ABC.coef(:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % ((start:end),x) <- along x-axis
			%}
			obj.ABC.coef(:,2) = coefABC2(obj.Courant(2),obj.ABC.type); % ((start:end),y) <- along y-axis
			obj.ABC.coef(:,1) = coefABC2(obj.Courant(1),obj.ABC.type); % ((start:end),x) <- along x-axis
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
			obj.ABC.coef(:,:,2) = coefABC2(obj.Courant(2),obj.ABC.type); % (3,start:end,y) <- along y-axis
			obj.ABC.coef(:,:,1) = coefABC2(obj.Courant(1),obj.ABC.type); % (3,start:end,x) <- along x-axis
		else
			obj.ABC.oldEy_x	= zeros(3,obj.Ny,2,2); % (m(1,2,3),n(all),prev:next,start:end) <- along x-axis
			obj.ABC.oldEx_y	= zeros(obj.Nx,3,2,2); % (m(all),n(1,2,3),prev:next,start:end) <- along x-axis
			% calculate coefficient on left-end and right-end of grid
			%{
			obj.ABC.coef(:,:,2) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (3,start:end) <- along y-axis
			obj.ABC.coef(:,:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (3,start:end) <- along x-axis
			%}
			obj.ABC.coef(:,:,2) = coefABC2(obj.Courant(2),obj.ABC.type); % (3,start:end) <- along y-axis
			obj.ABC.coef(:,:,1) = coefABC2(obj.Courant(1),obj.ABC.type); % (3,start:end) <- along x-axis
		end
end

end
