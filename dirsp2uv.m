% dirsp2uv.m
function [u,v]=dirsp2uv(dir,spd)
% input parameters:
%	dir = array of direction
%	spd = array of speed
%
% output parameters:
%	u = array of u components.
%	v = array of v components.
% convert arrays of directions and speeds dimensioned to jcount
% elements in each array to arrays of u and v coordinates
%
% llz 07/10/00.
%

     theta=dir;
     theta(dir>=0 & dir<=90)=(-1*(theta(dir>=0 & theta<=90)+90));
     theta(dir<0 | dir>90)=270-theta(dir<0 | theta>90);
     u=cos(theta*pi/180).*spd;
     v=sin(theta*pi/180).*spd;
