function drawskewt()

global drawmapflg hdfile datfile minpress maxpress num_snd num_prf title
global mbarbclr mflength mfintv mbarb_size sameframe gtlakesflg

barbclr = ['b';'r';'g';'m';'c';'y';'k';'w'];
maxlength = 200;
epoch = datenum(1970,1,2,0,0,0);

warning off all

SkewtDir = cd;
if datfile == 'txt' 
   datfile = ['*.' datfile]; 
   D       = dir(datfile);
   num_snd = length(D);
end


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
    if strcmp('tmp',datfile(len-2:len))
        % Loop to read FSL's ASCII format sounding data from the FSL Web site
        if j == 1
            [fid,message] = fopen(datfile,'rt');
            if fid == -1;
                msgbox('Can not open the data file!','Error Window','Error');
                return;
            end
        end
        %Reformat the data from the coded format
        [press,temp,dewp,u,v,lat,lon,alt,RaobInfo]=reformatFSL(fid);
        RaobInfo

        tk = temp + 273.16;
        %get the mixing ratio (g/kg), potential temperature (K), relative humidity,
        %vapor pressure (mb), and vapor density (g/m^3)
        for i=1:length(press)
            wl(i)  = wb(dewp(i),press(i));
            pt(i) = PotTemp(temp(i),press(i));
            rh(i)  = es(dewp(i))/es(temp(i));
            [e(i),rho(i)] = vapor(tk(i),rh(i),0);
        end

        num   = findstr(RaobInfo,':');
        info  = RaobInfo(1,(num(4)-9):length(RaobInfo));
        dmy   = RaobInfo(1,(num(1)+2):(num(2)-5));
        hrstr = RaobInfo(1,(num(2)+1):(num(3)-1));
        hr    = str2num(hrstr);

        if j==1
            begt  = (datenum(dmy)-epoch)*86400+hr*3600;
        end
        secnd(j) = (datenum(dmy)-epoch)*86400+hr*3600;
        RaobInfo = [];
        RaobInfo = strcat(dmy,hrstr,info);
        u        = u*1.94;
        v        = v*1.94;
        lenth = length(alt);
        if ~feof(fid)
            alt   = alt*1000.;
        else
            break
        end
    elseif strcmp('txt',datfile(len-2:len))
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
    end			% END of GETTING THE DATA
   
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
   
   %fill the arrays to be used in a map plot with a NaN
   if j == 1
      mu     = ones(num_snd,maxlength)*nan;
      mv     = ones(num_snd,maxlength)*nan;
      mtemp  = ones(num_snd,maxlength)*nan;
      mdewp  = ones(num_snd,maxlength)*nan;
      mlat   = ones(num_snd,maxlength)*nan;
      mlon   = ones(num_snd,maxlength)*nan;
      malt   = ones(num_snd,maxlength)*nan;
      mpress = ones(num_snd,maxlength)*nan;
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
   
   % Draw Graphs   
   orient landscape;
   if ~ sameframe
      % not same frame make a new figure for each sounding
      skewt = figure(j+1);
      set(skewt, 'Inverthardcopy','on', ...
         'Position',[120 150 880 500]);
      clf
   else
      % same frame do a new figure and the background only once
      if j == 1
         skewt = figure(2);
         set(skewt,'Position',[120 150 880 500]);
         set(skewt,'DoubleBuffer','off','BackingStore',...
                   'off','Renderer','zbuffer')
         clf;
         drawbg(minpress,maxpress,skewt)
         whitebg(skewt,[.89 .91 .93])
         drawnow
      end
   end
   if sameframe
      % same frame draw profiles
      drawgraph(minpress,maxpress,press,temp,dewp,u,v,j,RaobInfo)
      drawnow
   else
      % not same frame draw background and profiles
      drawbg(minpress,maxpress,skewt)
      whitebg(skewt,[.89 .91 .93])
      drawgraph(minpress,maxpress,press,temp,dewp,u,v,1,RaobInfo)
   end

   %Do the soundings on the map
   set(skewt, 'Inverthardcopy','off')
   dump        = ones(1,maxlength-length(u))*nan;
   mu(j,:)     = [u(1,:) dump];
   mv(j,:)     = [v(1,:) dump];
   mdewp(j,:)  = [dewp(1,:) dump];
   mtemp(j,:)  = [temp(1,:) dump];
   mlat(j,:)   = [lat(1,:) dump];
   mlon(j,:)   = [lon(1,:) dump];
   malt(j,:)   = [alt(1,:) dump];
   mpress(j,:) = [press(1,:) dump];
   if drawmapflg == 1
      if num_snd == 1 | ~ sameframe
         %determine lat/lon limits of the map
         minlat = min(lat);
         maxlat = max(lat);
         minlon = min(lon);
         maxlon = max(lon);
         latlim = [minlat-.75  maxlat+.75];
         lonlim = [minlon-.75  maxlon+.75];
         zlimit = [min(alt)-.1 max(alt)+.1];
         [mapax,map,maplegend] = plterrain(skewt,latlim,lonlim,zlimit);
         drawnow; cd(SkewtDir);
         % draw lines, stems and barbs
         barbhdl = drawbarbonmap(u,v,lat,lon,alt,press,mbarbclr,mapax,map,maplegend);
      end
      suptitle(title);
      if sameframe ~= 1
          orient(skewt,'portrait');
          print(gcf,'-djpeg','-r300','-painters','-adobecset',datfile(1:len-4));
      end
   end

   dump        = ones(1,maxlength-length(atemp))*nan;
   tt(j,:)     = [atemp  dump];
   dd(j,:)     = [adewp  dump];
   rr(j,:)     = [arho   dump];
   pp(j,:)     = [apres  dump];
   zz(j,:)     = [z      dump];

   clear press temp dewp wl u v rho lat lon atemp adewp arho apres
end     % END of For J Loop


if drawmapflg & sameframe
    %draw map with all the soundings on it
    minlat = min(mlat(:));
    maxlat = max(mlat(:));
    minlon = min(mlon(:));
    maxlon = max(mlon(:));
    latlim = [minlat-.75   maxlat+.75];
    lonlim = [minlon-.75   maxlon+.75];
    zlimit = [min(min(malt))-.5 max(max(malt)+.5)];
    [mapax,map,maplegend] = plterrain(skewt,latlim,lonlim,zlimit);
    drawnow;  cd(SkewtDir);
    for k = 1:j
        barbhdl=drawbarbonmap(mu(k,:),mv(k,:),mlat(k,:),mlon(k,:),...
            malt(k,:),mpress(k,:),barbclr(k),...
            mapax,map,maplegend);
    end
    suptitle(title);
    orient(skewt,'portrait');
    filename = ['AllSond_' datfile(15:len-4)]
    print(gcf,'-djpeg','-r300','-painters','-adobecset',filename);
end


MeanAtm(tt,dd,rr,pp,zz,j,datfile);

