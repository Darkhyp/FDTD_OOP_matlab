function source = SourceFunction(time,sourceStr,CourantNumber)
    switch sourceStr.Type
    case 'gauss_m'
	% modulated Gauss
        source = exp(-(CourantNumber/sourceStr.ppw*(time-sourceStr.timeDelay)).^2).* ...
		        sin(sourceStr.omega*time); % omega = (omega*dt)

    case 'gauss'
        source = exp(-(CourantNumber/sourceStr.ppw*(time-sourceStr.timeDelay)).^2);
    case 'Ricker'
        % The Ricker Wavelet
        arg = (pi*((CourantNumber*time - sourceStr.timeDelay)/sourceStr.ppw - 1)).^2;
        source = (1 - 2*arg).*exp(-arg);%.*sin(sourceStr.omega*time); % omega = (omega*dt)
    end
%     source = exp(-((time - delay - location / cdtds) / width)^2);
%     source = sin(2*pi/ppw * (cdtds * time - location));
%  arg = (pi * ((cdtds * time - location) / ppw - 1.0))^2;
%     source = (1 - 2*arg)*exp(-arg);
end