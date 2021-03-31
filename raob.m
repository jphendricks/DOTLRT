function raob()

maxlength = 200;

minpress = 100;
maxpress = 1030;

warning off all

datfile = 'txt';
datfile = ['*.' datfile]; 
D       = dir(datfile);
num_snd = length(D);

%FOR EACH SOUNDING
for j = 1:num_snd
    press = [];
    temp  = [];
    dewp  = [];
    u     = [];
    v     = [];
    lat   = [];
    lon   = [];
    lon   = [];
    alt   = [];

    %GET THE DATA to PLOT
    datfile = D(j).name;
    len     = length(datfile);
        % Loop to read NASA sounding files 
        [fid,message] = fopen(datfile,'rt');
        if fid == -1;
            msgbox('Can not open the data file!','Error Window','Error');
            return;
        end
        %Reformat the data from the coded format
        [press,temp,dewp,u,v,lat,lon,alt,RaobInfo]=reformatNASA(fid);
        RaobInfo

        tk = temp + 273.16;
        %get the mixing ratio (g/kg), potential temperature (K), relative humidity,
        %vapor pressure (mb), and vapor density (g/m^3)
        for i=1:length(press)
            wl(i)  = wb(dewp(i),press(i));
            pt(i)  = PotTemp(temp(i),press(i));
            rh(i)  = es(dewp(i))/es(temp(i));
            [e(i),rho(i)] = vapor(tk(i),rh(i),0);
        end
   % END of GETTING THE DATA
   
   % create a common height grid and interpolate profiles onto it
   z = alt(1):150:alt(length(alt));
   atemp = interp1(alt,temp,z);
   adewp = interp1(alt,dewp,z);
   arho  = interp1(alt,rho,z);
   apres = interp1(alt,press,z);
      
   %fill the arrays to be used in averaging with a NaN
   if j == 1
      tt  = ones(num_snd,maxlength)*nan;
      dd  = ones(num_snd,maxlength)*nan;
      rr  = ones(num_snd,maxlength)*nan;
      pp  = ones(num_snd,maxlength)*nan;
   end
   
   %get rid of bad data, i.e. outside the max and min pressure
   u(isnan(lat)      | press>maxpress | press<minpress)=[];
   v(isnan(lat)      | press>maxpress | press<minpress)=[];
   temp(isnan(lat)   | press>maxpress | press<minpress)=[];
   dewp(isnan(lat)   | press>maxpress | press<minpress)=[];
   alt(isnan(lat)    | press>maxpress | press<minpress)=[];
   lon(isnan(lat)    | press>maxpress | press<minpress)=[];
   lat2=lat;
   lat(isnan(lat)    | press>maxpress | press<minpress)=[];
   press(isnan(lat2) | press>maxpress | press<minpress)=[];
   clear lat2

   dump        = ones(1,maxlength-length(u))*nan;
   mu(j,:)     = [u(1,:) dump];
   mv(j,:)     = [v(1,:) dump];
   mdewp(j,:)  = [dewp(1,:) dump];
   mtemp(j,:)  = [temp(1,:) dump];
   mlat(j,:)   = [lat(1,:) dump];
   mlon(j,:)   = [lon(1,:) dump];
   malt(j,:)   = [alt(1,:) dump];
   mpress(j,:) = [press(1,:) dump];

   dump        = ones(1,maxlength-length(atemp))*nan;
   tt(j,:)     = [atemp  dump];
   dd(j,:)     = [adewp  dump];
   rr(j,:)     = [arho   dump];
   pp(j,:)     = [apres  dump];
   zz(j,:)     = [z      dump];

   clear press temp dewp wl u v rho lat lon atemp adewp arho apres
end     % END of For J Loop

%MeanAtm(tt,dd,rr,pp,zz,j,datfile);
%Average p,t,td,rho, to create mean atmosphere
for ii = 1:length(pp)
    zm(ii)  = mean(zz(:,ii));
    pm(ii)  = mean(pp(:,ii));
    tcm(ii) = mean(tt(:,ii));
    dm(ii)  = mean(dd(:,ii));
    rm(ii)  = mean(rr(:,ii));
end
pm(isnan(zm)) = []; tcm(isnan(zm)) = []; dm(isnan(zm)) = []; rm(isnan(zm)) = []; zm(isnan(zm)) = [];
zpl = (zm-zm(1))/1000;

tm = tcm + 273.16;

