class AngleValueType(object):
  def __init__(self, radians=None, degrees=None):
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
    return math.sin(self.degrees)

  @property
  def cos(self):
    """Returns the cos value of the angle.
    """
    return math.cos(self.degrees)

class SolarCalculation(object):
  
  def __init__(self, local_time, latitude, longitude, is_daylight_saving_time):
    self.local_time = local_time
    self.latitude = latitude
    self.longitude = longitude
    self.is_daylight_saving_time = is_daylight_saving_time

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
    lstm = 15 * round(self.longitude / 15)
    ast = lst + timedelta(
        minutes=(4*(lstm-self.longitude)+self.equation_of_time))
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
    return AngularValueType(degrees=hour_angle) 

  @property
  def altitude_angle(self):
    """Returns the altitude angle.

    The solar altitude angle (β) is the apparent angular height of the sun in the
    sky if you are facing it.
                sin(β)=cos(L)cos(δ)cos(H)+sin(L)sin(δ)
    where
      L = latitude(positive in either hemisphere)
      δ = declination angle (negative for Southern Hemisphere)[-23.45 deg to +23.45 def]
      H = hour angle
    
    Returns
      The altitude angle as an AngleValueType.
    """
