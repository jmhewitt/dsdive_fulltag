daytime = function(date, lat, lon) {
  # Uses the suncalc package to determine if it is daytime at a given location 
  # wrt. civil sunrise and sunset times.  So, daytime is any time where the 
  # sun's altitude is greater than -6 degrees.
  # 
  # Parameters:
  #  date - vector of DateTime objects (can be a Date or character in UTC, 
  #    or a POSIXct) at which daylight will be evaluated
  #  lat - latitude of observation
  #  lon - longitude of observation
  
  getSunlightPosition(
    date = date, lat = lat, lon = lon, keep = 'altitude'
  )[,'altitude'] > - 6 * pi / 180
}
