% es.m
function [vapor]=es(t)

% input parameters:
%	t = temperature
% output parameters:
%	vapor = water vapor saturation presslre (mb)
%
% This function calculates water vapor saturation presslre (mb) 
% from temperature (kelvin) to better than 3 per cent accuracy
%
% llz 07-07-2000
%

t=t+273.16;
z = 2.302585;
ts = 373.16;
tice = 273.16;
if (t-233.16) < 0
   ei1=-9.09718*(tice/t-1.)*z;
   ei2=-3.56654*log(tice/t);
   ei3=.876793*(1.0-t/tice)*z;
   vapor=exp(ei1+ei2+ei3)*6.1071;
else
   es1=-7.90298*(ts/t-1.)*z;
   es2=5.02808*log((ts/t));
   es3=(-1.3816e-7)*power(10.,(11.344*(1.-t/ts))-1.0)*z;
   es4=(8.1328e-3)*power(10.,(3.49149*(ts/t-1.0))-1.0)*z;
   vapor=exp(es1+es2+es3+es4)*1013.246;
end
