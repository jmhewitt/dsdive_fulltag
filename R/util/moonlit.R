moonlit = function(date, lat, lon) {
  # Uses the suncalc package to determine if there is bright moon illumination 
  # at a given location by determining if a) the moon is above the horizon, and 
  # b) the moon is at least 50% illuminated.
  # 
  # Parameters:
  #  date - vector of DateTime objects (can be a Date or character in UTC, 
  #    or a POSIXct) at which daylight will be evaluated
  #  lat - latitude of observation
  #  lon - longitude of observation
  
  moon_up = getMoonPosition(
    date = date, lat = lat, lon = lon, keep = 'altitude'
  )[,'altitude'] > 0
  
  moon_bright = getMoonIllumination(
    date = date, keep = 'fraction'
  )[,'fraction'] > .5
 
  moon_up & moon_bright  
}
