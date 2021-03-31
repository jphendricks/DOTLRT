% wb.m
function [mixratio]=wb(td,p)
% input parameters:
%	td = temperature or dew point;
%	p = pressure;
% output parameters:
%	mixratio = mixing ratio 
% This function returns the mixing ratio (grams of water vapor per
% kilogram of dry air) given the dew point(celsius) 
% and pressure (millibars). If the temperature is input 
% instead of the dew point, then saturation mixing ratio (same units) 
% is returned.
%
% llz 07-07-2000
%
  [x]=es(td);
  mixratio = 621.97*x/(p-x);
