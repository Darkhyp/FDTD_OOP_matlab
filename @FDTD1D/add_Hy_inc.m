function obj = add_Hy_inc(obj,Position,time,sourceStr)
%#codegen

global imp0

source = SourceFunction(time,sourceStr,obj.Courant(1)) / imp0;

if sourceStr.isAdd
    obj.Hy(Position(1)) = obj.Hy(Position(1)) - source;
else
    obj.Hy(Position(1)) = -source;
end
fprintf('%d, %g\n',time,obj.Hy(Position(1)))

end
