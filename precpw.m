% precpw.m
function [pcpw]=precpw(td,p,n)
%
% input parameters:
%	td=dew point (celsius)
%	p =pressure (millibars)
%	n =number of levels
% output parameters:
%	pcpw = precipitable water precpw (cm)
%
% This function computes total precipitable water precpw (cm) in a 
% vertical column of air based upon sounding data at n levels.
% Calculations are done in cgs units.
%
% llz 07-07-2000
%

% initialize value of precipitable water
pw=0.;
nl=n-1;

% calculte the mixing ratio at the lowest level
wbot = wb(td(1),p(1));

for i=1:nl
   wtop = wb(td(i+1),p(i+1));
   % calculate the layer-mean mixing ratio (g/kg)
   ww = .5*(wtop+wbot);
   % make the mixing ratio dimensionless
   wl=.001*ww;
   % calculate the specific humidity.
   ql=wl/(wl+1);
   % the factor of 1000. below converts from millibars to dynes/cm**2.
   dp = 1000.*(p(i)-p(i+1));
   pw = pw+(ql/980.616)*dp;
   wbot=wtop;
end
pcpw = pw;

