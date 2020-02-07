rm(list = ls(all.names = TRUE))
library(data.table)

##PARAM
#Nice: 43.655217, 7.211717
#Spa: 50.457312, 5.953582
#Paris: 48.825447, 2.344051
#Tours: 47.372513, 0.711754

lat = 49.78126; #latitude (decimal degrees)
lon = 15.90820; #longitude (decimal degrees)
elevTarget = 20; # (degrees from horizontal, +above, -below)
aziTarget = 180; # (degrees from North, east to west)
zone = 0; # GMT + zone ? keep 0 because downloaded file is UTC
heatingTemp = 15; # °C
mirrorArea = 1; # m² #0.2827 m²
removeSummer = FALSE; #remove period outside of startHeating -> stopHeating
startHeating = '15/10'; #DD/MM #only for north hemisphere
stopHeating = '15/04'; #DD/MM
takeBuildingInertia = FALSE; #to take into account thermal inertia of building
gvalue = 0.75;  #transmissivity of the whole solar spectrum through glazing
# https://energie.wallonie.be/fr/07-06-facteur-solaire-g.html?IDC_PEB=9491&IDD=113659&IDC=9094
## END PARAM

# download data of a typical climatologic year, whose year is written as 2011
df <- read.table(paste0('http://re.jrc.ec.europa.eu/api/tmy?lat=',lat,'&lon=',lon), sep=",", skip=16, nrows=8760, header=TRUE)

lat = lat*pi/180;
df$datenum = ISOdatetime( 2011, substr(df$time.UTC., 5, 6), substr(df$time.UTC., 7, 8), substr(df$time.UTC., 10, 11),00,00, tz = "GMT")
df$dayOfYear = as.integer( (df$datenum - df$datenum[1])/86400 + 1 )

startHeating = ISOdate(2011,substr(startHeating, 4, 5),substr(startHeating, 1, 2), hour = 0, min = 0, sec = 0, tz = "GMT" )
startHeating = as.integer( (startHeating - df$datenum[1]) + 1 )
stopHeating = ISOdate(2011,substr(stopHeating, 4, 5),substr(stopHeating, 1, 2), hour = 0, min = 0, sec = 0, tz = "GMT" )
stopHeating = as.integer( (stopHeating - df$datenum[1]) + 1 )

##supposed temp in building
df$tempBuilding = NA
for (i in (1+48):(8760) ) {
  df$tempBuilding[i] = mean( df$T2m[(i-24):i] )*0.67 + mean( df$T2m[(i-48):(i-25)] )*0.33
}

# dates to keep
if (removeSummer == TRUE){
  df = df[ df$dayOfYear < stopHeating | df$dayOfYear > startHeating , ]
}

df$decimalHour = as.integer(substr(df$time.UTC., 10, 11))/24
df$decli = 23.45*pi/180*sin(2*pi*(df$dayOfYear+284)/365.25); #solar declination in rad

## sun position
df$B = 2 * pi / 364 * (df$dayOfYear-81);
df$MST = df$decimalHour + (lon - zone * 15) / 361; # Mean solar time
df$EOT = (9.87 * sin(2*df$B) - 7.53 * cos(df$B) - 1.5 * sin(df$B)) / 1440; # Equation of time
df$t = df$MST + df$EOT ;
df$t[df$t < 0] = 1-abs(df$t[df$t < 0]);

df$elevSun = asin(sin(df$decli)*sin(lat)-cos(df$decli)*cos(lat)*cos(2*pi*df$t)); # solaire elevation in rad

#save only when sun is above horizon
df$elevSun[ df$elevSun <= 0 ] = NA;

# https://en.wikipedia.org/wiki/Solar_azimuth_angle, aziSun in degrees
df$aziSun = 180/pi*acos((sin(df$decli)-sin( df$elevSun )*sin(lat))/(cos( df$elevSun )*cos(lat))); #+morning; -afternoon

df$aziSun[df$t > 0.5 && df$t < 1] = 360 - df$aziSun ;
# equivalent to
# df$aziSun[df$t > 0.5 & df$t < 1 & !is.na(df$aziSun)] = 360 - df$aziSun[df$t > 0.5 & df$t < 1 & !is.na(df$aziSun)] ;

df$aziSun[df$t == 0.5] = 180;
df$aziSun[df$t == 0.0] = 0;
df$aziSun[df$t == 1.0] = 360;

## transfer to rad
df$compl2elevSun = (pi/2 - df$elevSun);
df$aziSun4comput = df$aziSun * pi/180;
compl2elevTarget = (90 - elevTarget) * pi/180;
aziTarget = aziTarget* pi/180;

# imagine the sun is due east 45° elevation, and the target is due west 0° elevation.
# the mirror azimuth will not be the mean azimuth of sun and target (180°)
# that's why we need these 4 lines from cerebralmeltdown.com, file Jim_s_Arduino_SunTracker_V99_2/Functions.pde, function FindHeliostatAngle()
df$yz = acos( cos(df$compl2elevSun) * cos(compl2elevTarget) + sin(df$compl2elevSun) * sin(compl2elevTarget) * cos(df$aziSun4comput-aziTarget) );
df$ya = asin( sin(df$compl2elevSun) * ( sin(df$aziSun4comput-aziTarget) / sin(df$yz) ) );
df$compl2elevMirror = acos( cos(compl2elevTarget) * cos(df$yz/2) + sin(compl2elevTarget) * sin(df$yz/2) * cos(df$ya) );
df$aziMirror = asin( sin(df$yz/2) * ( sin(df$ya)/sin(df$compl2elevMirror) ) ) + aziTarget;

# efficiency = cos(incidence)
# https://fr.wikipedia.org/wiki/Panneau_solaire , paragraph "Angle d'incidence du soleil"
# reminding that sin(compl_alpa) = cos(alpha)
df$efficiency = cos(df$compl2elevMirror)*cos(df$compl2elevSun) + sin(df$compl2elevMirror)*sin(df$compl2elevSun)*cos(df$aziSun4comput-df$aziMirror);

# tempBuilding and efficiency have NA
if (takeBuildingInertia == TRUE) {
  id = df$tempBuilding < heatingTemp & !is.na(df$tempBuilding) & !is.na(df$efficiency);
  print(sum(id))
  df$result[id] = df$Gb.n.[id] * df$efficiency[id]
} else {
  id = df$T2m < heatingTemp & !is.na(df$T2m) & !is.na(df$efficiency)
  df$result[id] = df$Gb.n.[id] * df$efficiency[id]
}

kwh = sum(df$result, na.rm = TRUE)*mirrorArea / 1000 * gvalue;

res = matrix(data = NA, nrow = 2, ncol = 4);
rownames(res) <- c('price','co2_kg');
colnames(res) <- c('gaz','fuel','elecFR','elecBE');

### https://elyotherm.fr/comparatif-cout-energies-kwh
res['price','gaz'] = kwh*0.086;
res['price','fuel'] = kwh*0.101;
res['price','elecFR'] = kwh*0.1714; #gouv.fr 2018
res['price','elecBE'] = kwh*0.2659; #gouv.fr 2018
res['co2_kg','gaz'] = kwh*0.229;
res['co2_kg','fuel'] = kwh*0.3;
res['co2_kg','elecFR'] = kwh*0.05599; #electricityMap json
res['co2_kg','elecBE'] = kwh*0.23252; #electricityMap json

# https://www.ademe.fr/sites/default/files/assets/documents/avis-ademe-modes-chauffage-individuels.pdf
# https://www.quelleenergie.fr/magazine/nouvelles-energie/quels-types-de-chauffage-energie-2017/


print(paste( 'energy received from heliostat = ', kwh, 'kwh'))
print(res)


