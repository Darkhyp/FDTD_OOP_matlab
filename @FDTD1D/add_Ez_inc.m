%% calculate source function at given time and location
function obj = add_Ez_inc(obj,Position,time,sourceStr)
%             if (sourceStr.ppw <= 0)
%                 fprintf(stderr,'E_z(inc): ezIncInit() must be called before ezInc.\n Points per wavelength must be positive.\n');
%                 exit(-1);
%             end

source = SourceFunction(time,sourceStr,obj.Courant(1));
if sourceStr.isAdd
    obj.Ez(Position(1)) = obj.Ez(Position(1)) + source;
else
    obj.Ez(Position(1)) = source;
end
fprintf('%d, %g\n',time,obj.Ez(Position(1)))

end
