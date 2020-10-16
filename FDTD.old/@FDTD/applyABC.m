%% apply Absorbing boundary conditioms
function obj = applyABC(obj)
    if isempty(obj.ABC) % check if initABC() has been called
        fprintf('\nUse initABC() method before applyABC()!...')
%                 exit(-1);
    end
    switch obj.Dimensionality
        case 1
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
        case 2
            switch obj.ABC.type
                case 1
                    % first-order ABC.
                    if obj.isTE
                        % ABC at x1 (start) and x2 (end)
                        m_1 = 1;
                        m_2 = obj.Nx;
                        obj.Ez(m_1,:)	= obj.ABC.oldEz_x(1,:) + obj.ABC.coef(1)*(obj.Ez(m_1+1,:) - obj.Ez(m_1,:));
                        obj.Ez(m_2,:)	= obj.ABC.oldEz_x(2,:) + obj.ABC.coef(2)*(obj.Ez(m_2-1,:) - obj.Ez(m_2,:));
                        % update stored fields
                        obj.ABC.oldEz_x = obj.Ez([2,obj.Nx-1],:);

                        % ABC at y1 (start) and y2 (end)
                        n_1 = 1;
                        n_2 = obj.Ny;
                        obj.Ez(:,n_1)	= obj.ABC.oldEz_y(:,1) + obj.ABC.coef(1)*(obj.Ez(:,n_1+1) - obj.Ez(:,n_1));
                        obj.Ez(:,n_2)	= obj.ABC.oldEz_y(:,2) + obj.ABC.coef(2)*(obj.Ez(:,n_2-1) - obj.Ez(:,n_2));
                        % update stored fields
                        obj.ABC.oldEz_y = obj.Ez(:,[2,obj.Ny-1]);
                    else
                        % ABC at x1 (start) and x2 (end)
                        m_1 = 1;
                        m_2 = obj.Nx;
                        obj.Ey(m_1,:)	= obj.ABC.oldEy_x(1,:) + obj.ABC.coef(1)*(obj.Ey(m_1+1,:) - obj.Ey(m_1,:));
                        obj.Ey(m_2,:)	= obj.ABC.oldEy_x(2,:) + obj.ABC.coef(2)*(obj.Ey(m_2-1,:) - obj.Ey(m_2,:));
                        % update stored fields
                        obj.ABC.oldEy_x = obj.Ey([2,obj.Nx-1],:);

                        % ABC at y1 (start) and y2 (end)
                        n_1 = 1;
                        n_2 = obj.Ny;
                        obj.Ex(:,n_1)	= obj.ABC.oldEx_y(:,1) + obj.ABC.coef(1)*(obj.Ex(:,n_1+1) - obj.Ex(:,n_1));
                        obj.Ex(:,n_2)	= obj.ABC.oldEx_y(:,2) + obj.ABC.coef(2)*(obj.Ex(:,n_2-1) - obj.Ex(:,n_2));
                        % update stored fields
                        obj.ABC.oldEx_y = obj.Ex(:,[2,obj.Ny-1]);
                    end
                case 2
                    % Second-order ABC.
                    if obj.isTE
                        % ABC at x1 side of grid
                        obj.Ez(1,:) = ...
                            obj.ABC.coef(1,1)*(obj.ABC.oldEz_x(1,:,2,1) +  obj.Ez(3,:)) ...
                          + obj.ABC.coef(2,1)*(obj.ABC.oldEz_x(1,:,1,1) + obj.ABC.oldEz_x(3,:,1,1) ...
                                             - obj.ABC.oldEz_x(2,:,2,1) -  obj.Ez(2,:)) ...
                          + obj.ABC.coef(3,1)* obj.ABC.oldEz_x(2,:,1,1) - obj.ABC.oldEz_x(3,:,2,1);
                        % memorize old fields
                        obj.ABC.oldEz_x(:,:,2,1) = obj.ABC.oldEz_x(:,:,1,1);
                        obj.ABC.oldEz_x(:,:,1,1) = obj.Ez(1:3,:);

                        % ABC at x2 side of grid
                        m = obj.Nx;
                        obj.Ez(m,:) = ...
                            obj.ABC.coef(1,2)*(obj.ABC.oldEz_x(1,:,2,2) +  obj.Ez(m-2,:)) ...
                          + obj.ABC.coef(2,2)*(obj.ABC.oldEz_x(1,:,1,2) + obj.ABC.oldEz_x(3,:,1,2) ...
                                             - obj.ABC.oldEz_x(2,:,2,2) -  obj.Ez(m-1,:)) ...
                          + obj.ABC.coef(3,2)* obj.ABC.oldEz_x(2,:,1,2) - obj.ABC.oldEz_x(3,:,2,2);
                        % memorize old fields
                        obj.ABC.oldEz_x(:,:,2,2) = obj.ABC.oldEz_x(:,:,1,2);
                        obj.ABC.oldEz_x(:,:,1,2) = obj.Ez(m-(0:2),:);

                        % ABC at y1 of grid
                        obj.Ez(:,1) = ...
                            obj.ABC.coef(1,1)*(obj.ABC.oldEz_y(:,1,2,1) +  obj.Ez(:,3)) ...
                          + obj.ABC.coef(2,1)*(obj.ABC.oldEz_y(:,1,1,1) + obj.ABC.oldEz_y(:,3,1,1) ...
                                             - obj.ABC.oldEz_y(:,2,2,1) -  obj.Ez(:,2)) ...
                          + obj.ABC.coef(3,1)* obj.ABC.oldEz_y(:,2,1,1) - obj.ABC.oldEz_y(:,3,2,1);
                        % memorize old fields
                        obj.ABC.oldEz_y(:,:,2,1) = obj.ABC.oldEz_y(:,:,1,1);
                        obj.ABC.oldEz_y(:,:,1,1) = obj.Ez(:,1:3);

                        % ABC at y2 of grid
                        n = obj.Ny;
                        obj.Ez(:,n) = ...
                            obj.ABC.coef(1,2)*(obj.ABC.oldEz_y(:,1,2,2) +  obj.Ez(:,n-2)) ...
                          + obj.ABC.coef(2,2)*(obj.ABC.oldEz_y(:,1,1,2) + obj.ABC.oldEz_y(:,3,1,2) ...
                                             - obj.ABC.oldEz_y(:,2,2,2) -  obj.Ez(:,n-1)) ...
                          + obj.ABC.coef(3,2)* obj.ABC.oldEz_y(:,2,1,2) - obj.ABC.oldEz_y(:,3,2,2);
                        % memorize old fields
                        obj.ABC.oldEz_y(:,:,2,2) = obj.ABC.oldEz_y(:,:,1,2);
                        obj.ABC.oldEz_y(:,:,1,2) = obj.Ez(:,n-(0:2));
                    else
                        % ABC at x1 of grid
                        obj.Ey(1,:) = ...
                            obj.ABC.coef(1,1)*(obj.ABC.oldEy_x(1,:,2,1) +  obj.Ey(3,:)) ...
                          + obj.ABC.coef(2,1)*(obj.ABC.oldEy_x(1,:,1,1) + obj.ABC.oldEy_x(3,:,1,1) ...
                                             - obj.ABC.oldEy_x(2,:,2,1) -  obj.Ey(2,:)) ...
                          + obj.ABC.coef(3,1)* obj.ABC.oldEy_x(2,:,1,1) - obj.ABC.oldEy_x(3,:,2,1);
                        % memorize old fields
                        obj.ABC.oldEy_x(:,:,2,1) = obj.ABC.oldEy_x(:,:,1,1);
                        obj.ABC.oldEy_x(:,:,1,1) = obj.Ey(1:3,:);

                        % ABC at x2 of grid
                        m = obj.Nx;
                        obj.Ey(m,:) = ...
                            obj.ABC.coef(1,2)*(obj.ABC.oldEy_x(1,:,2,2) +  obj.Ey(m-2,:)) ...
                          + obj.ABC.coef(2,2)*(obj.ABC.oldEy_x(1,:,1,2) + obj.ABC.oldEy_x(3,:,1,2) ...
                                             - obj.ABC.oldEy_x(2,:,2,2) -  obj.Ey(m-1,:)) ...
                          + obj.ABC.coef(3,2)* obj.ABC.oldEy_x(2,:,1,2) - obj.ABC.oldEy_x(3,:,2,2);
                        % memorize old fields
                        obj.ABC.oldEy_x(:,:,2,2) = obj.ABC.oldEy_x(:,:,1,2);
                        obj.ABC.oldEy_x(:,:,1,2) = obj.Ey(m-(0:2),:);

                        % ABC at y1 of grid
                        obj.Ex(:,1) = ...
                            obj.ABC.coef(1,1)*(obj.ABC.oldEx_y(:,1,2,1) +  obj.Ex(:,3)) ...
                          + obj.ABC.coef(2,1)*(obj.ABC.oldEx_y(:,1,1,1) + obj.ABC.oldEx_y(:,3,1,1) ...
                                             - obj.ABC.oldEx_y(:,2,2,1) -  obj.Ex(:,2)) ...
                          + obj.ABC.coef(3,1)* obj.ABC.oldEx_y(:,2,1,1) - obj.ABC.oldEx_y(:,3,2,1);
                        % memorize old fields
                        obj.ABC.oldEx_y(:,:,2,1) = obj.ABC.oldEx_y(:,:,1,1);
                        obj.ABC.oldEx_y(:,:,1,1) = obj.Ex(:,1:3);

                        % ABC at y2 of grid
                        n = obj.Ny;
                        obj.Ex(:,n) = ...
                            obj.ABC.coef(1,2)*(obj.ABC.oldEx_y(:,1,2,2) +  obj.Ex(:,n-2)) ...
                          + obj.ABC.coef(2,2)*(obj.ABC.oldEx_y(:,1,1,2) + obj.ABC.oldEx_y(:,3,1,2) ...
                                             - obj.ABC.oldEx_y(:,2,2,2) -  obj.Ex(:,n-1)) ...
                          + obj.ABC.coef(3,2)* obj.ABC.oldEx_y(:,2,1,2) - obj.ABC.oldEx_y(:,3,2,2);
                        % memorize old fields
                        obj.ABC.oldEx_y(:,:,2,2) = obj.ABC.oldEx_y(:,:,1,2);
                        obj.ABC.oldEx_y(:,:,1,2) = obj.Ex(:,n-(0:2));
                    end
            end
        case 3
            switch obj.ABC.type
                case 1  % first-order ABC.
                    % ABC at x1 (start) and x2 (end)
                    m = 1;
                    obj.Ey(m,:,:) = obj.ABC.oldEy_x(1,:,:) + obj.ABC.coef(1,1)*(obj.Ey(m+1,:,:) - obj.Ey(m,:,:));
                    obj.Ez(m,:,:) = obj.ABC.oldEz_x(1,:,:) + obj.ABC.coef(1,2)*(obj.Ez(m+1,:,:) - obj.Ez(m,:,:));
                    m = obj.Nx;
                    obj.Ey(m,:,:) = obj.ABC.oldEy_x(2,:,:) + obj.ABC.coef(2,1)*(obj.Ey(m-1,:,:) - obj.Ey(m,:,:));
                    obj.Ez(m,:,:) = obj.ABC.oldEz_x(2,:,:) + obj.ABC.coef(2,2)*(obj.Ez(m-1,:,:) - obj.Ez(m,:,:));
                    obj.ABC.oldEy_x(:,:,:) = obj.Ey([2,obj.Nx-1],:,:);
                    obj.ABC.oldEz_x(:,:,:) = obj.Ez([2,obj.Nx-1],:,:);

                    % ABC at y1 (start) and y2 (end)
                    n = 1;
                    obj.Ez(:,n,:) = obj.ABC.oldEz_y(:,1,:) + obj.ABC.coef(1,3)*(obj.Ez(:,n+1,:) - obj.Ez(:,n,:));
                    obj.Ex(:,n,:) = obj.ABC.oldEx_y(:,1,:) + obj.ABC.coef(1,4)*(obj.Ex(:,n+1,:) - obj.Ex(:,n,:));
                    n = obj.Ny;
                    obj.Ez(:,n,:) = obj.ABC.oldEz_y(:,2,:) + obj.ABC.coef(2,3)*(obj.Ez(:,n-1,:) - obj.Ez(:,n,:));
                    obj.Ex(:,n,:) = obj.ABC.oldEx_y(:,2,:) + obj.ABC.coef(2,4)*(obj.Ex(:,n-1,:) - obj.Ex(:,n,:));
                    obj.ABC.oldEz_y(:,:,:) = obj.Ez(:,[2, obj.Ny-1],:);
                    obj.ABC.oldEx_y(:,:,:) = obj.Ex(:,[2, obj.Ny-1],:);

                    % ABC at z1 (start) and z2 (end)
                    p = 1;
                    obj.Ex(:,:,p) = obj.ABC.oldEx_z(:,:,1) + obj.ABC.coef(1,5)*(obj.Ex(:,:,p+1) - obj.Ex(:,:,p));
                    obj.Ey(:,:,p) = obj.ABC.oldEy_z(:,:,1) + obj.ABC.coef(1,6)*(obj.Ey(:,:,p+1) - obj.Ey(:,:,p));
                    p = obj.Nz;
                    obj.Ex(:,:,p) = obj.ABC.oldEx_z(:,:,2) + obj.ABC.coef(2,5)*(obj.Ex(:,:,p-1) - obj.Ex(:,:,p));
                    obj.Ey(:,:,p) = obj.ABC.oldEy_z(:,:,2) + obj.ABC.coef(2,6)*(obj.Ey(:,:,p-1) - obj.Ey(:,:,p));
                    obj.ABC.oldEx_z(:,:,:) = obj.Ex(:,:,[2, obj.Nz-1]);
                    obj.ABC.oldEy_z(:,:,:) = obj.Ey(:,:,[2, obj.Nz-1]);
                case 2 % Second-order ABC.
                    % ABC at x1 side of grid
                    obj.Ey(1,:,:) = ...
                        obj.ABC.coef(1,1,1)*(obj.ABC.oldEy_x(1,:,:,2,1) +  obj.Ey(3,:,:)) ...
                      + obj.ABC.coef(2,1,1)*(obj.ABC.oldEy_x(1,:,:,1,1) + obj.ABC.oldEy_x(3,:,:,1,1) ...
                                           - obj.ABC.oldEy_x(2,:,:,2,1) -  obj.Ey(2,:,:)) ...
                      + obj.ABC.coef(3,1,1)* obj.ABC.oldEy_x(2,:,:,1,1) - obj.ABC.oldEy_x(3,:,:,2,1);
                    obj.Ez(1,:,:) = ...
                        obj.ABC.coef(1,1,2)*(obj.ABC.oldEz_x(1,:,:,2,1) +  obj.Ez(3,:,:)) ...
                      + obj.ABC.coef(2,1,2)*(obj.ABC.oldEz_x(1,:,:,1,1) + obj.ABC.oldEz_x(3,:,:,1,1) ...
                                           - obj.ABC.oldEz_x(2,:,:,2,1) -  obj.Ez(2,:,:)) ...
                      + obj.ABC.coef(3,1,2)* obj.ABC.oldEz_x(2,:,:,1,1) - obj.ABC.oldEz_x(3,:,:,2,1);
                    % memorize old fields
                    obj.ABC.oldEy_x(:,:,:,2,1) = obj.ABC.oldEy_x(:,:,:,1,1);
                    obj.ABC.oldEy_x(:,:,:,1,1) = obj.Ey(1:3,:,:);
                    obj.ABC.oldEz_x(:,:,:,2,1) = obj.ABC.oldEz_x(:,:,:,1,1);
                    obj.ABC.oldEz_x(:,:,:,1,1) = obj.Ez(1:3,:,:);

                    % ABC at x2 side of grid
                    m = obj.Nx;
                    obj.Ey(m,:,:) = ...
                        obj.ABC.coef(1,2,1)*(obj.ABC.oldEy_x(1,:,:,2,2) +  obj.Ey(m-2,:,:)) ...
                      + obj.ABC.coef(2,2,1)*(obj.ABC.oldEy_x(1,:,:,1,2) + obj.ABC.oldEy_x(3,:,:,1,2) ...
                                           - obj.ABC.oldEy_x(2,:,:,2,2) -  obj.Ey(m-1,:,:)) ...
                      + obj.ABC.coef(3,2,1)* obj.ABC.oldEy_x(2,:,:,1,2) - obj.ABC.oldEy_x(3,:,:,2,2);
                    obj.Ez(m,:,:) = ...
                        obj.ABC.coef(1,2,2)*(obj.ABC.oldEz_x(1,:,:,2,2) +  obj.Ez(m-2,:,:)) ...
                      + obj.ABC.coef(2,2,2)*(obj.ABC.oldEz_x(1,:,:,1,2) + obj.ABC.oldEz_x(3,:,:,1,2) ...
                                           - obj.ABC.oldEz_x(2,:,:,2,2) -  obj.Ez(m-1,:,:)) ...
                      + obj.ABC.coef(3,2,2)* obj.ABC.oldEz_x(2,:,:,1,2) - obj.ABC.oldEz_x(3,:,:,2,2);
                    % memorize old fields
                    obj.ABC.oldEy_x(:,:,:,2,2) = obj.ABC.oldEy_x(:,:,:,1,2);
                    obj.ABC.oldEy_x(:,:,:,1,2) = obj.Ey(m-(0:2),:,:);
                    obj.ABC.oldEz_x(:,:,:,2,2) = obj.ABC.oldEz_x(:,:,:,1,2);
                    obj.ABC.oldEz_x(:,:,:,1,2) = obj.Ez(m-(0:2),:,:);

                    % ABC at y1 of grid
                    obj.Ex(:,1,:) = ...
                        obj.ABC.coef(1,1,3)*(obj.ABC.oldEx_y(:,1,:,2,1) +  obj.Ex(:,3,:)) ...
                      + obj.ABC.coef(2,1,3)*(obj.ABC.oldEx_y(:,1,:,1,1) + obj.ABC.oldEx_y(:,3,:,1,1) ...
                                           - obj.ABC.oldEx_y(:,2,:,2,1) -  obj.Ex(:,2,:)) ...
                      + obj.ABC.coef(3,1,3)* obj.ABC.oldEx_y(:,2,:,1,1) - obj.ABC.oldEx_y(:,3,:,2,1);
                    obj.Ez(:,1,:) = ...
                        obj.ABC.coef(1,1,4)*(obj.ABC.oldEz_y(:,1,:,2,1) +  obj.Ez(:,3,:)) ...
                      + obj.ABC.coef(2,1,4)*(obj.ABC.oldEz_y(:,1,:,1,1) + obj.ABC.oldEz_y(:,3,:,1,1) ...
                                           - obj.ABC.oldEz_y(:,2,:,2,1) -  obj.Ez(:,2,:)) ...
                      + obj.ABC.coef(3,1,4)* obj.ABC.oldEz_y(:,2,:,1,1) - obj.ABC.oldEz_y(:,3,:,2,1);
                    % memorize old fields
                    obj.ABC.oldEx_y(:,:,:,2,1) = obj.ABC.oldEx_y(:,:,:,1,1);
                    obj.ABC.oldEx_y(:,:,:,1,1) = obj.Ex(:,1:3,:);
                    obj.ABC.oldEz_y(:,:,:,2,1) = obj.ABC.oldEz_y(:,:,:,1,1);
                    obj.ABC.oldEz_y(:,:,:,1,1) = obj.Ez(:,1:3,:);

                    % ABC at y2 of grid
                    n = obj.Ny;
                    obj.Ex(:,n,:) = ...
                        obj.ABC.coef(1,2,3)*(obj.ABC.oldEx_y(:,1,:,2,2) +  obj.Ex(:,n-2,:)) ...
                      + obj.ABC.coef(2,2,3)*(obj.ABC.oldEx_y(:,1,:,1,2) + obj.ABC.oldEx_y(:,3,:,1,2) ...
                                           - obj.ABC.oldEx_y(:,2,:,2,2) -  obj.Ex(:,n-1,:)) ...
                      + obj.ABC.coef(3,2,3)* obj.ABC.oldEx_y(:,2,:,1,2) - obj.ABC.oldEx_y(:,3,:,2,2);
                    obj.Ez(:,n,:) = ...
                        obj.ABC.coef(1,2,4)*(obj.ABC.oldEz_y(:,1,:,2,2) +  obj.Ez(:,n-2,:)) ...
                      + obj.ABC.coef(2,2,4)*(obj.ABC.oldEz_y(:,1,:,1,2) + obj.ABC.oldEz_y(:,3,:,1,2) ...
                                           - obj.ABC.oldEz_y(:,2,:,2,2) -  obj.Ez(:,n-1,:)) ...
                      + obj.ABC.coef(3,2,4)* obj.ABC.oldEz_y(:,2,:,1,2) - obj.ABC.oldEz_y(:,3,:,2,2);
                    % memorize old fields
                    obj.ABC.oldEx_y(:,:,:,2,2) = obj.ABC.oldEx_y(:,:,:,1,2);
                    obj.ABC.oldEx_y(:,:,:,1,2) = obj.Ex(:,n-(0:2),:);
                    obj.ABC.oldEz_y(:,:,:,2,2) = obj.ABC.oldEz_y(:,:,:,1,2);
                    obj.ABC.oldEz_y(:,:,:,1,2) = obj.Ez(:,n-(0:2),:);

                    % ABC at z1 of grid
                    obj.Ex(:,:,1) = ...
                        obj.ABC.coef(1,1,5)*(obj.ABC.oldEx_z(:,:,1,2,1) +  obj.Ex(:,:,3)) ...
                      + obj.ABC.coef(2,1,5)*(obj.ABC.oldEx_z(:,:,1,1,1) + obj.ABC.oldEx_z(:,:,3,1,1) ...
                                           - obj.ABC.oldEx_z(:,:,2,2,1) -  obj.Ex(:,:,2)) ...
                      + obj.ABC.coef(3,1,5)* obj.ABC.oldEx_z(:,:,2,1,1) - obj.ABC.oldEx_z(:,:,3,2,1);
                    obj.Ey(:,:,1) = ...
                        obj.ABC.coef(1,1,6)*(obj.ABC.oldEy_z(:,:,1,2,1) +  obj.Ey(:,:,3)) ...
                      + obj.ABC.coef(2,1,6)*(obj.ABC.oldEy_z(:,:,1,1,1) + obj.ABC.oldEy_z(:,:,3,1,1) ...
                                           - obj.ABC.oldEy_z(:,:,2,2,1) -  obj.Ey(:,:,2)) ...
                      + obj.ABC.coef(3,1,6)* obj.ABC.oldEy_z(:,:,2,1,1) - obj.ABC.oldEy_z(:,:,3,2,1);
                    % memorize old fields
                    obj.ABC.oldEx_z(:,:,:,2,1) = obj.ABC.oldEx_z(:,:,:,1,1);
                    obj.ABC.oldEx_z(:,:,:,1,1) = obj.Ex(:,:,1:3);
                    obj.ABC.oldEy_z(:,:,:,2,1) = obj.ABC.oldEy_z(:,:,:,1,1);
                    obj.ABC.oldEy_z(:,:,:,1,1) = obj.Ey(:,:,1:3);

                    % ABC at z2 of grid
                    p = obj.Nz;
                    obj.Ex(:,:,p) = ...
                        obj.ABC.coef(1,2,5)*(obj.ABC.oldEx_z(:,:,1,2,2) +  obj.Ex(:,:,p-2)) ...
                      + obj.ABC.coef(2,2,5)*(obj.ABC.oldEx_z(:,:,1,1,2) + obj.ABC.oldEx_z(:,:,3,1,2) ...
                                           - obj.ABC.oldEx_z(:,:,2,2,2) -  obj.Ex(:,:,p-1)) ...
                      + obj.ABC.coef(3,2,5)* obj.ABC.oldEx_z(:,:,2,1,2) - obj.ABC.oldEx_z(:,:,3,2,2);
                    obj.Ey(:,:,p) = ...
                        obj.ABC.coef(1,2,6)*(obj.ABC.oldEy_z(:,:,1,2,2) +  obj.Ey(:,:,p-2)) ...
                      + obj.ABC.coef(2,2,6)*(obj.ABC.oldEy_z(:,:,1,1,2) + obj.ABC.oldEy_z(:,:,3,1,2) ...
                                           - obj.ABC.oldEy_z(:,:,2,2,2) -  obj.Ey(:,:,p-1)) ...
                      + obj.ABC.coef(3,2,6)* obj.ABC.oldEy_z(:,:,2,1,2) - obj.ABC.oldEy_z(:,:,3,2,2);
                    % memorize old fields
                    obj.ABC.oldEx_z(:,:,:,2,2) = obj.ABC.oldEx_z(:,:,:,1,2);
                    obj.ABC.oldEx_z(:,:,:,1,2) = obj.Ex(:,:,p-(0:2));
                    obj.ABC.oldEy_z(:,:,:,2,2) = obj.ABC.oldEy_z(:,:,:,1,2);
                    obj.ABC.oldEy_z(:,:,:,1,2) = obj.Ey(:,:,p-(0:2));
            end
    end
end
