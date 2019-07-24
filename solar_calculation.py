class AngleValueType(object):
  def __init__(self, radians=None, degrees=None, ):
    """Constructs a AngularValueType instance.

	  The constructor should be called with either radians=x or degrees=y and the 
    constructor will calculate the complementary. If none is given or both are 
    given but do not match, a ValueError will be thrown.

    Args:
      radians: the radian value. 
      degrees: the degree value. 
		"""
    if not radians and not degrees:
      raise ValueError('Radians and Degress both not given')

    if radians and not degrees:
      self.radians = radians
      self.degrees = math.degrees(radians) 
    elif degrees and not radians:
      self.degrees = degrees
      self.radians = math.radians(degrees)
    else:
      if abs(math.degrees(radians) - degrees) > 0.01:
        raise ValueError("Radians and Degrees do not agree") 
      self.radians = radians
      self.degrees = degrees
  
  @property
  def sin(self):
    """Returns the sin value of the angle.
    """
    return math.sin(self.radians)

  @property
  def cos(self):
    """Returns the cos value of the angle.
    """
    return math.cos(self.radians)

class SolarCalculation(object):
  
  def __init__(self, local_time, latitude, longitude, is_daylight_saving_time,
               apparent_extraterrestrial_solar_intensity=None, 
               atmospheric_extinction_coefficient=None):
    """Initializes a SolarCalculation instance.
    
    Args:
      local_time: A datetime instance, the local time.
      latitude: A float value, the latitude in degrees
      longitude: A float value, the longitude in degrees.
      is_daylight_saving_time: A bool value, whether use dailight saving time.
      apparent_extraterrestrial_solar_intensity: A float number, the apparent 
          extraterrestrial solar intensity in Btu/(hr*ft^2). Its necessary to 
          calculate solar irradiance flux.
      atmospheric_extinction_coefficient: A float number, the atmospheric 
          extinction coefficient. Its necessary to calculate solar irradiance 
          flux.
    """
    self.local_time = local_time
    self.latitude = AngleValueType(degrees=latitude)
    self.longitude = AngleValueType(degrees=longitude) 
    self.is_daylight_saving_time = is_daylight_saving_time
    self.apparent_extraterrestrial_solar_intensity = (
        apparent_extraterrestrial_solar_intensity)
    self.atmospheric_extinction_coefficient = atmospheric_extinction_coefficient

  @property
  def day_of_year(self):
    day_of_year = self.local_time.timetuple().tm_yday
    return day_of_year

  @property
  def equation_of_time(self):
    """Reuturns the equation of time in minutes.

    The difference between the true solar time and the mean solar time changes 
    continuously day-to-day with an annual cycle. This quantity is known as the
    equation of time. The equation of time, ET in minutes is approximated by:

               $ET = 9.87\sin(2D)-7.53\cos(D)-1.5\sin(D)$

    where $D = 360^{o}\frac{(N-81)}{365}$ and N is the day of year. nterest
    LSTM = local longitude of standard t:ime mer

    Returns:
      The equation of time in minutes.
    """
    degrees = (self.day_of_year - 81.0) * (360.0/365.0)
    radians = math.radians(degrees)
    equation_of_time_minutes = (9.87 * math.sin(2*radians) - 
                                7.53 * math.cos(radians) - 
                                1.5 * math.sin(radians))
    return equation_of_time_minutes

  @property
  def declination_angle(self):
    """Returns the declination angel.

    The declination is the angular distance of the sun north or south of the 
    earth's equator. The declination angle, δ, for the Northern Hemisphere 
    (reverse the declination angle's sign for the Southern Hemisphere) is:
          $\delta = 23.45^{o}\sin[\frac{N+284}{365}\times360^{o}]$
    where N is the day number of the year, with January 1st. equal to 1.
    """
    declination_angle = 23.45 * AngleValueType(
        degrees=(self.day_of_year + 284.) / 365. * 360.).sin
    return AngleValueType(degrees=declination_angle)

  @property
  def apparent_solar_time(self):
    """Returns the apparent solar time.

    The apparent solar time, AST (or local solar time) in the western longitudes 
    is calculated from

               $AST = LST + (4min/deg)(LSTM - Long) + ET$

    where:
    LST = local standard time or clock time for that time zone (may need to 
          adjust for daylight savings time, DST, that is LST = DST – 1 hr)
    Long = local longitude at the position of interest
    LSTM = local longitude of standard time meridian
              $LSTM = 15^{o} \times (\frac{Long}{15^{o}})$
    
    Returns:i
      The apparent solar time as a datetime object.
    """
    if self.is_daylight_saving_time:
      lst = self.local_time- timedelta(hours=1)
    else:
      lst = self.local_time
    lstm = 15 * round(self.longitude.degrees / 15)
    ast = lst + timedelta(
        minutes=(4*(lstm-self.longitude.degrees)+self.equation_of_time))
    return ast
  
  @property
  def hour_angle(self):
    """Returns the hour angle.

    The hour angle, H, is the azimuth angle of the sun's rays caused by the 
    earth's rotation, and H can be computed from:

    $H = \fact{(No. of minutes past midnight, AST)-720mins}{4min/deg}$

    The hour angle defined here is negative in the morning and positive in the 
    afternoon.

    Returns:
      The hour angle as a AngleValueType.
    """
    num_mins_past_midnight = (self.apparent_solar_time.hour * 60. + 
                              self.apparent_solar_time.minute)
    hour_angle = (num_mins_past_midnight - 720.) / 4. 
    return AngleValueType(degrees=hour_angle) 

  @property
  def altitude_angle(self):
    """Returns the altitude angle.

    The solar altitude angle (β1) is the apparent angular height of the sun in the
    sky if you are facing it.
                sin(β1)=cos(L)cos(δ)cos(H)+sin(L)sin(δ)
    where
      L = latitude(positive in either hemisphere)
      δ = declination angle (negative for Southern Hemisphere)[-23.45 deg to +23.45 def]
      H = hour angle
    
    Returns
      The altitude angle as an AngleValueType.
    """
    # TODO(wenxuang): handel negative case.
    altitude_angle = math.asin(self.latitude.cos * self.declination_angle.cos * 
                               self.hour_angle.cos + self.latitude.sin * 
                               self.declination_angle.sin)
    return AngleValueType(radians=altitude_angle)
    
  @property
  def azimuth_angle(self):
    """Returns the solar azimuth angle.

    The solar azimuth, α1, is the angle away from south (north in the Southern 
    Hemisphere). 
               cos(α1) = (sin(β1)sin(L)-sin(δ))/(cos(β1)cos(L))
    Returns:
      The solar azimuth as an AngleValueType.
    """
    sign_hour_angle = 1 if self.hour_angle.degrees > 0 else -1
    azimuth_angle = math.acos(
        (self.altitude_angle.sin * self.latitude.sin - self.declination_angle.sin) 
        / (self.altitude_angle.cos * self.latitude.cos)) * sign_hour_angle
    return AngleValueType(radians=azimuth_angle)

  def collector_angle(self, surface_azimuth_angle, surface_tilt_angle):
    """Returns the collector angle of a surce at this location and time.

    The collector angle (θ) between the sun and normal to the surface is:
         cos(θ) = sin(β1)cos(β2) + cos(β1)sin(β2)cos(α1-α2)
    where α2 is the azimuth angle normal to the collector surface, and β2 is the
    tilt angle from the ground.

    Args:
      surface_azimuth_angle: An AngleValueType instance, the azimuth angle 
          normal to the collector surface.
      surface_tilt_angle: An AnglevalueType instance, the tilt angle from the 
          ground.
    
    Returns: 
      An AngleValueType instance, the collector angle between the sun and normal 
      to the surface.
    
    Raises:
      ValueError: If the surface_azimuth_angle and surface tilt_angle are not 
          AngleValueType instances.
    """
    if (type(surface_azimuth_angle) is not AngleValueType or 
        type(surface_tilt_angle) is not AngleValueType):
      raise ValueError("The angles must are AngleValueType") 
    
    collector_angle = math.acos(
        (self.altitude_angle.sin * surface_tilt_angle.cos) + 
        (self.altitude_angle.cos * surface_tilt_angle.sin) * math.cos(
            math.radians(self.azimuth_angle.degrees - 
                         surface_azimuth_angle.degrees)))
    
    return AngleValueType(radians=collector_angle)

  def direct_normal_irradiance(self, altitude):
    """Calculates and returns the direct normal irradiance to the ground.

    the direct normal irradiance to the ground is:
    
        $I_{DN} = A \exp(-\frac{p}{p_{0}}\frac{B}{\sin(\beta_{1})})$
    
    where A is the apparent extraterrestrial solar intensity, B is the 
    atmospheric extinction coefficient (mainly due to changes in atmospheric 
    moisture), and p/p0 is the pressure at the location of interest relative to 
    a standard atmosphere, given by p/p0 = exp(-0.0000361z) where z is the 
    altitude in feet above sea level.

    Args:
      altitude: A float number, the altitude in feet above the sea level.
    
    Returns:
      The direct normal irradiance as a float number in Btu/(hr*ft^2).
    
    Raises:
      ValueError: If the apparent extraterrestrial solar intensity or the 
          atmospheric extinction coefficien has not been initialied.
    """
    a = self.apparent_extraterrestrial_solar_intensity
    b = self.atmospheric_extinction_coefficient
    if not a or not b:
      raise ValueError("The apparent extraterrestrial solar intensity or the atmospheric extinction coefficien has not been initialied.")
    direct_normal_irradiance = a * math.exp(
        -math.exp(-0.0000361 * altitude) * b / self.altitude_angle.sin)
    return direct_normal_irradiance

  def direct_radiation_onto_collector(self, surface_azimuth_angle, 
                                      surface_tilt_angle, altitude):
    """Returns the direct radiation flux onto a collector.

    The direct radiation flux onto the collector is:
                      $I_{D} = I_{DN}\cos(\theta)$
    Args:
      surface_azimuth_angle: An AngleValueType instance, the azimuth angle 
          normal to the collector surface.
      surface_tilt_angle: An AnglevalueType instance, the tilt angle from the 
          ground.
      altitude: A float number, the altitude in feet above the sea level.    
    """
    direct_radiation_onto_collector = (
        self.collector_angle(surface_azimuth_angle, surface_tilt_angle).cos * 
        self.direct_normal_irradiance(altitude))
    return direct_radiation_onto_collector

  def diffuse_scattered_radiation(self, surface_azimuth_angle, 
                                 surface_tilt_angle, altitude):
    """Returns the diffuse-scattered radiation flux onto a collector.

    The direct radiation flux onto the collector is:
            $I_{DS} = C \times I_{DN}[\frac{1+\cos(\beta_{2})}{2}]$
    where C is the ratio of diffuse radiation on a horizontal surface to the 
    direct normal irradiation. 

    Args:
      surface_azimuth_angle: An AngleValueType instance, the azimuth angle 
          normal to the collector surface.
      surface_tilt_angle: An AnglevalueType instance, the tilt angle from the 
          ground.
      altitude: A float number, the altitude in feet above the sea level.    
    """
    dn = self.direct_normal_irradiance(altitude)
    c = 0.136  # constant  
    ds = c * dn * (1 + surface_tilt_angle.cos) / 2
    return ds
