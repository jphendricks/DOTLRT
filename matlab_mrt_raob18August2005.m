% matlab_mrt_raob.m

% get raob data
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

%Average p,t,td,rho, to create mean atmosphere
for ii = 1:length(pp)
    zm(ii)  = mean(zz(:,ii));
    pm(ii)  = mean(pp(:,ii));
    tcm(ii) = mean(tt(:,ii));
    dm(ii)  = mean(dd(:,ii));
    rm(ii)  = mean(rr(:,ii));
end
pm(isnan(zm)) = []; tcm(isnan(zm)) = []; dm(isnan(zm)) = []; rm(isnan(zm)) = []; zm(isnan(zm)) = [];
zpl = (zm-zm(1));

tm = tcm + 273.15; % convert from C to K
nz = length(tm);

% DOTLRTv1.0   INPUT
% observation parameters
inp_height = 6.5;  % observation height (km AGL)
inp_theta = 25;    % nadir angle (degrees)

% instrument parameters
num_sb_freqs = 1;  % number of frequency poins per sideband
instr_spec(1) = 166; % chan_lo_freq
instr_spec(2) = 0;   % chan_if1_freq
instr_spec(3) = 0;  % chan_if2_freq
instr_spec(4) = 6000; % chan_bandwidth
instr_spec(5) = 0.07; % chan_dtrms

% surface parameters
num_streams = 8; % total number of stream angles
surface_reflectivity = 0.05 * ones(num_streams/2);

% atmosphere parameters
for j = 1 : nz
  atm_inp(j,1) = zpl(j) / 1000; % level altitude (km AGL)
  atm_inp(j,2) = pm(j); % layer pressure (mb)
  atm_inp(j,3) = tm(j); % layer temperature (K)
  atm_inp(j,4) = rm(j); % layer water vapor density (g/m^3)
  atm_inp(j,5) = 0; % cloud liquid density (g/m^3)
  atm_inp(j,6) = 0; % rain density (g/m^3)
  atm_inp(j,7) = 0; % ice density (g/m^3)
  atm_inp(j,8) = 0; % snow density (g/m^3)
  atm_inp(j,9) = 0; % graupel density (g/m^3)
end

% DOTLRTv1.0    OUTPUT

% Tb_inp(jpol) = brightness temperature for observation level at observation angle
%    jpol = 1 (horizontal), 2 (vertical)
% tau(jud) = opacity for observation level at observation angle
%     jud  = 1 (upwelling), 2 (downwelling)
% Tb[jz,jstream,jpol] = brightness temperature profile at observation angle
%    jz = 1 ... nz+1    at levels   atm_inp_altitude(jz-1)  (first level at surface)
%    jud  = 1 (upwelling), 2 (downwelling)
%    jpol = 1 (horizontal), 2 (vertical)
% dTb_dT(jz,jstream,jpol), dTb_dp(1:nz,jstream,jpol), dTb_dq(1:nz,jstream,jpol), dTb_dw(1:nz,jstream,jhyd,jpol)
%     = Jacobian for observation level at observation angle
%    jz = 1, nz   at levels     atm_inp_altitude(jz)
%    jhyd  = 1 ... rnum_hydro_phases
%     jstream  = 1 ... num_streams
%     jpol = 1 (horizontal), 2 (vertical)
% stream_angles( 1 ... num_streams )

[Tb_inp,tau,Tb,dTb_dT,dTb_dp,dTb_dq,dTb_dw,stream_angles] = mex_mrt( inp_height, inp_theta, num_sb_freqs, instr_spec, ...
                                                                     num_streams, surface_reflectivity, nz, atm_inp);