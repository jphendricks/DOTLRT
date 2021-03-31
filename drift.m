% drift.m
function [lat,lon,times]=drift(time0,lat0,lon0,n,z,u,v,misflag)
% input parameters:
%	time0=starting time.
%	lat0,lon0 = initial latitude/longitude.
%	n = number of levels in raob.
%	misflag = missing value for z,u,v arrays.
%	z = press, heights, in meters.
%	u/v = u/v winds, in m/s unit.
% output parameters:
%	lat, lon = computed lat/lon at each level of report.
%	times = computed time at each level.
%
% 07-07-2000 llz
%
% This function computes the latitude and longitude of a rawinsonde by
% assuming a constant ascent rate of 4.58 m/s and by using the reported winds,
% in addition, it computes a time for each level with reported winds by 
% adding to the actual release time or by using an assumed release time 
% of 50 mins before the synoptic time.

size=500;
ascrat = 4.58;

if lat0==misflag | lon0==misflag
   return
end

ihr=time0/100;
imin=time0-(ihr*100);

% compute the longitude increment based on the latitude of the station.
dxlon = 111200.*cos(lat0*pi/180.);
lat(1)=lat0;
lon(1)=lon0;
times(1)=time0;

for i=2:n
   iflag=0;
   if u(i)==misflag | v(i)==misflag | z(i)==misflag
      iflag=1;
   end
   if iflag==1
      lat(i)=lat(i-1);
      lon(i)=lon(i-1);
      times(i)=times(i-1);
   else	% find a level below the starting level with good data
      for j=(i-1):-1:1
         if j==0
            lat(i)=lat0;
            lon(i)=lon0;
            times(i)=time0;
         else
            jflag=0;
            if u(j)==misflag | v(j)==misflag | z(j)==misflag
               jflag=1;
            end
            if jflag~=1		% compute a new lat/lon/time for level i.
               dt=(z(i)-z(j))/ascrat;
               imin=imin+(dt/60. + .5);
               if imin>=60
                  ihr=ihr+1;
                  if ihr == 24
                     ihr==0;
                     imin=imin-60;
                  end
               end
               times(i)=ihr*100+imin;
               lat(i)=lat(j)+.5*(v(i)+v(j))*dt/111200.;
               lon(i)=lon(j)+.5*(u(i)+u(j))*dt/dxlon;
               break
            end
         end
      end
   end
end
