%% apply Absorbing boundary conditioms
function obj = applyABC(obj)
%#codegen

if isempty(obj.ABC) % check if initABC() has been called
    fprintf('\nUse initABC() method before applyABC()!...')
%                 exit(-1);
end

switch obj.ABC.type
    case 0    % simple ABC for left side of grid
        obj.Ez(1) = obj.Ez(2);
    case 1
        % ABC for left and right sides of grid
        obj.Ez([1,obj.Nx]) = obj.ABC.oldEz_x + obj.ABC.coef.*(obj.Ez([2,obj.Nx+1-2]) - obj.Ez([1,obj.Nx]));
        % update stored fields
        obj.ABC.oldEz_x = obj.Ez([2,obj.Nx-1]);
    case 2 % Second-order ABC.
        % ABC for left side of grid
        obj.Ez(1) = ...
            obj.ABC.coef(1,1) * (obj.ABC.oldEz_x(1,2,1) +  obj.Ez(3)) ...
          + obj.ABC.coef(2,1) * (obj.ABC.oldEz_x(1,1,1) + obj.ABC.oldEz_x(3,1,1) ...
                               - obj.ABC.oldEz_x(2,2,1) -  obj.Ez(2)) ...
          + obj.ABC.coef(3,1) *  obj.ABC.oldEz_x(2,1,1) - obj.ABC.oldEz_x(3,2,1);
        % update stored fields
        obj.ABC.oldEz_x1(:,2) = obj.ABC.oldEz_x(:,1,1);
        obj.ABC.oldEz_x1(:,1) = obj.Ez(1:3);

        % ABC for right side of grid
        m = obj.Nx;
        obj.Ez(m) = ...
            obj.ABC.coef(1,2) * (obj.ABC.oldEz_x(1,2,2) +  obj.Ez(m-2)) ...
          + obj.ABC.coef(2,2) * (obj.ABC.oldEz_x(1,1,2) + obj.ABC.oldEz_x(3,1,2) ...
                               - obj.ABC.oldEz_x(2,2,2) -  obj.Ez(m-1)) ...
          + obj.ABC.coef(3,2) *  obj.ABC.oldEz_x(2,1,2) - obj.ABC.oldEz_x(3,2,2);
        % update stored fields
        obj.ABC.oldEz_x(:,2,2) = obj.ABC.oldEz_x(:,1,2);
        obj.ABC.oldEz_x(:,1,2) = obj.Ez(m-(0:2));
end

end
