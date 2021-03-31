function [press,temp,dpt,u,v,lat,lon,ht,raobinfo]=reformatNASA(fid)

% input parameters:
%	fid = pointer to data file.
% output parameters:
%	press
%	temperature
%	dewp = dew points
%	u    = u-component of witn in kt
%	v    = v-component of witn in kt
%	lat  = latitude
%	lon  = longitude
%	alt  = elevation
% This function reads a NASA ASCII RAOB data file, and gets press, height,
% temperature, dew points, direction and speeds.
% Then it calculates u and v from direction and speed, calculates
% lat and lon of the drift.
% bbs 12/06/04

%Conversion from m to feet
C1 = 3.2808399;
%Conversion from kt to m/s
C2 = 0.5144;

for i=1:3
    line = fgetl(fid);
end
line = fgetl(fid);
site = line(1:7);

switch site
case 'POINT M'
    lat0 = 34.06; lon0 = -119.06; site = 'Pt. Mugu';
case 'SAN NIC'
    lat0 = 33.26; lon0 = -119.49; site = 'San Nicolas Island';
case 'EDWARDS'
   lat0 = 34.9; lon0 = -117.92;   site = 'Edwards Air Force Base';
case 'VANDENB'   
   lat0 = 34.75;lon0 = -120.57;   site = 'Vandenberg Air Force Base'; 
end

line = fgetl(fid);
time  = line(1:4);
if site(1) ~= 'E'
    date  = line(8:18);
else
    date  = line(8:16);
end
for i=1:4
    line = fgetl(fid);
end

ic = 0;
for i = 1:200
    line = fgetl(fid);
    ic = ic + 1;    
    if site(1) == 'E'
        if findstr(line(1:11),'TERMINATION')
            break
        end
    else
        if findstr(line(1:11),'           ')
            break
        end
    end
    dat(ic,1) = str2num(line(1,1:6));
    dat(ic,6) = str2num(line(1,7:10));
    dat(ic,7) = str2num(line(1,11:16));
    dat(ic,3) = str2num(line(1,22:27));
    dat(ic,4) = str2num(line(1,28:33));
    dat(ic,2) = str2num(line(1,34:41));
    dat(ic,5) = str2num(line(1,42:45));
    q(ic)     = str2num(line(1,46:51));
    rho(ic)   = str2num(line(1,52:59));
end
nl = length(dat)-1;
for i = 1:100
    line = fgetl(fid);
    if length(line) == 0
        line
    elseif findstr(line(1:11),'SIGNIFICANT') > 0
        line = fgetl(fid);
        line = fgetl(fid);
        line = fgetl(fid);
        break
    end
end
ic = nl;
for i = 1:200
    line = fgetl(fid);
    if findstr(line(1:11),'           ')
        break
    elseif findstr(line(1:11),'TERMINATION') > 0
       break 
    end
    ic = ic + 1;
    dat(ic,1) = str2num(line(1,1:6));
    dat(ic,6) = str2num(line(1,7:10));
    dat(ic,7) = str2num(line(1,11:14));
    dat(ic,3) = str2num(line(1,15:20));
    dat(ic,4) = str2num(line(1,21:26));
    dat(ic,2) = str2num(line(1,27:33));
    dat(ic,5) = str2num(line(1,38:41));
end
nl    = length(dat);
dat   = sortrows(dat,1);
ht    = dat(:,1)'/C1;
press = dat(:,2)';
temp  = dat(:,3)';
dpt   = dat(:,4)';
rh    = dat(:,5)';
Dir   = dat(:,6)';
Sp    = dat(:,7)'/C2;

dht = diff(ht); 

fclose(fid);
clear message i line dum

[u,v]   = dirsp2uv(Dir,Sp);
u(Dir==999)=0; v(Dir==999)=0;
[vapor] = precpw(dpt,press,nl);
times = [];  rtime = str2num(time); badflag = -999;
[lat,lon,times] = drift(rtime,lat0,lon0,nl,ht,u,v,badflag);
raobinfo = [site ':  ' date ', ' time  ' UTC   PWV=' num2str(vapor) 'cm'];

