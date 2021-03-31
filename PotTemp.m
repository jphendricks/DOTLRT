	function [pottemp] = pt(t,p)

%this function returns potential temperature (celsius) given
%temperature t (celsius) and pressure p (mb) by solving the poisson
%equation.

tk      = t+273.15;
pottemp = tk*((1000./p)^0.286);
