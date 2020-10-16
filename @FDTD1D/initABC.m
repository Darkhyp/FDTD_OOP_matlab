%% initialize Absorbing boundary conditioms
function obj = initABC(obj,type)
%#codegen

switch type
    case '1abc' % first order ABC
        obj.ABC.type = 1;
        obj.ABC.oldEz_x	= [0, 0];
    otherwise % second order ABC
        obj.ABC.type = 2;
        obj.ABC.oldEz_x = zeros(3,2,2);
end                   
%                     obj.ABC.coef = (obj.Courant-1)/(obj.Courant+1);
% calculate coefficient on left-end and right-end of grid
%                     obj.ABC.coef = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end)
obj.ABC.coef = coefABC2(obj.Courant(1),obj.ABC.type); % (3,start:end)

end
