Tutorial
========

Basic usage with an example dataset. Note, this tutorial uses the
matplotlib graphing module which is *not* a dependency of ``AeroEvap``,
be sure to install it to your environment before running this tutorial
if you want the plots to display correctly.

    >>> import pandas as pd
    >>> import numpy as np
    >>> from aeroevap import Aero
    >>> from IPython.display import IFrame
    >>> import matplotlib.pyplot as plt
    >>> %matplotlib inline

Example data
------------

This example uses buoy data from a location near Savannah, GA (NOAA station ID
is 41008). The buoy is maintained by the National Data Buoy Center (NDBC), more
buoy information is shown in the embedded page below. The meterological data
used in this example is hosted by NOAA and downloaded directly and formatted
for a month of data.

.. raw:: html
    
    <iframe
        width="700"
        height="500"
        src="https://www.ndbc.noaa.gov/station_page.php?station=41008"
        frameborder="0"
        allowfullscreen
    ></iframe>


The line below downloads the time series of current year buoy standard
meterological data directly from the NDBC.

Input units:

==== ==== === ==== === === === ==== ==== ==== ==== === ====
WDIR WSPD GST WVHT DPD APD MWD PRES ATMP WTMP DEWP VIS TIDE
==== ==== === ==== === === === ==== ==== ==== ==== === ====
degT m/s  m/s m    sec sec deg hPa  degC degC degC nmi ft
==== ==== === ==== === === === ==== ==== ==== ==== === ====

    
    >>> # get standard meterological data from National Data Buoy Center
    >>> met_df = pd.read_csv(
    >>>     'https://www.ndbc.noaa.gov/data/l_stdmet/41008.txt', 
    >>>     delim_whitespace=True, skiprows=[1], na_values=[999.0]
    >>> )

Make a datetime index and clean up the dataframe.


    >>> met_df.index = pd.to_datetime(
    >>>     dict(
    >>>         year=met_df['#YY'], 
    >>>         month=met_df.MM, 
    >>>         day=met_df.DD, 
    >>>         hour=met_df.hh,
    >>>         minute=met_df.mm
    >>>     )
    >>> )
    >>> met_df.index.name = 'date'
    >>> met_df.drop(['#YY','MM','DD','hh','mm'], axis=1, inplace=True)
    >>> met_df.head()




.. raw:: html

     <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }

        .dataframe tbody tr th {
            vertical-align: top;
        }

        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>WDIR</th>
          <th>WSPD</th>
          <th>GST</th>
          <th>WVHT</th>
          <th>DPD</th>
          <th>APD</th>
          <th>MWD</th>
          <th>PRES</th>
          <th>ATMP</th>
          <th>WTMP</th>
          <th>DEWP</th>
          <th>VIS</th>
          <th>TIDE</th>
        </tr>
        <tr>
          <th>date</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>2020-01-31 23:50:00</td>
          <td>359</td>
          <td>10.1</td>
          <td>11.8</td>
          <td>1.80</td>
          <td>8.33</td>
          <td>5.01</td>
          <td>101.0</td>
          <td>1016.8</td>
          <td>11.0</td>
          <td>13.9</td>
          <td>10.4</td>
          <td>99.0</td>
          <td>99.0</td>
        </tr>
        <tr>
          <td>2020-02-01 00:50:00</td>
          <td>349</td>
          <td>8.1</td>
          <td>9.6</td>
          <td>1.68</td>
          <td>5.56</td>
          <td>4.99</td>
          <td>69.0</td>
          <td>1017.6</td>
          <td>10.6</td>
          <td>13.9</td>
          <td>10.0</td>
          <td>99.0</td>
          <td>99.0</td>
        </tr>
        <tr>
          <td>2020-02-01 01:50:00</td>
          <td>13</td>
          <td>8.5</td>
          <td>10.1</td>
          <td>1.61</td>
          <td>7.69</td>
          <td>5.10</td>
          <td>105.0</td>
          <td>1016.4</td>
          <td>10.1</td>
          <td>14.0</td>
          <td>9.5</td>
          <td>99.0</td>
          <td>99.0</td>
        </tr>
        <tr>
          <td>2020-02-01 02:50:00</td>
          <td>24</td>
          <td>7.8</td>
          <td>9.1</td>
          <td>1.68</td>
          <td>7.14</td>
          <td>5.17</td>
          <td>103.0</td>
          <td>1015.8</td>
          <td>10.1</td>
          <td>14.1</td>
          <td>9.5</td>
          <td>99.0</td>
          <td>99.0</td>
        </tr>
        <tr>
          <td>2020-02-01 03:50:00</td>
          <td>35</td>
          <td>7.1</td>
          <td>9.0</td>
          <td>1.59</td>
          <td>6.67</td>
          <td>5.13</td>
          <td>103.0</td>
          <td>1015.6</td>
          <td>10.6</td>
          <td>14.2</td>
          <td>10.0</td>
          <td>99.0</td>
          <td>99.0</td>
        </tr>
      </tbody>
    </table>
    </div>
    <br>



Because the input dataset does not include relative humitidy we can
estimate it using an approximation to the Clausius–Clapeyron relation
using air and dewpoint temperatures. Relative humitidy is needed in the
aerodynamic mass-transfer evaporation calculations.


    >>> # vapor pressure and saturation vapor pressure using Clausius–Clapeyron relation
    >>> met_df['e'] = 0.611 * np.exp( 5423 * ((1/273) - (1/(met_df.DEWP+273.15))) )
    >>> met_df['es'] = 0.611 * np.exp( 5423 * ((1/273) - (1/(met_df.ATMP+273.15))) )


    >>> # calculate relative humitidy
    >>> met_df['RH'] = 100 * (met_df.e/met_df.es)
    >>> plt.figure(figsize=(8,4))
    >>> met_df.RH.plot()
    >>> plt.ylabel('estimated relative humitidy')


.. figure:: _static/RH.png

In this case we do *not* need to convert air pressure to millibars
because 1 hPa = 1 mbar.

Create an :obj:`.Aero` object
-----------------------------

The :obj:`.Aero` object allows for loading a :obj:`pandas.DataFrame` containing
meterological data required for calculating aerodynamic mass-transfer
open water evaporation in parrallel. The object can be initialized from
a :obj:`pandas.DataFrame` or the :obj:`pandas.DataFrame` can be assigned
later, e.g.


    >>> Aero_empty = Aero()
    >>> Aero_with_df = Aero(met_df)


    >>> Aero_empty.df is None
        True



    >>> # the df property can be assigned after initialization:
    >>> Aero_empty.df = met_df


    >>> # the data has been added
    >>> Aero_empty.df.head()


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }

        .dataframe tbody tr th {
            vertical-align: top;
        }

        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>WDIR</th>
          <th>WSPD</th>
          <th>GST</th>
          <th>WVHT</th>
          <th>DPD</th>
          <th>APD</th>
          <th>MWD</th>
          <th>PRES</th>
          <th>ATMP</th>
          <th>WTMP</th>
          <th>DEWP</th>
          <th>VIS</th>
          <th>TIDE</th>
          <th>e</th>
          <th>es</th>
          <th>RH</th>
        </tr>
        <tr>
          <th>date</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>2020-01-31 23:50:00</td>
          <td>359</td>
          <td>10.1</td>
          <td>11.8</td>
          <td>1.80</td>
          <td>8.33</td>
          <td>5.01</td>
          <td>101.0</td>
          <td>1016.8</td>
          <td>11.0</td>
          <td>13.9</td>
          <td>10.4</td>
          <td>99.0</td>
          <td>99.0</td>
          <td>1.279457</td>
          <td>1.332185</td>
          <td>96.042019</td>
        </tr>
        <tr>
          <td>2020-02-01 00:50:00</td>
          <td>349</td>
          <td>8.1</td>
          <td>9.6</td>
          <td>1.68</td>
          <td>5.56</td>
          <td>4.99</td>
          <td>69.0</td>
          <td>1017.6</td>
          <td>10.6</td>
          <td>13.9</td>
          <td>10.0</td>
          <td>99.0</td>
          <td>99.0</td>
          <td>1.245352</td>
          <td>1.296822</td>
          <td>96.031065</td>
        </tr>
        <tr>
          <td>2020-02-01 01:50:00</td>
          <td>13</td>
          <td>8.5</td>
          <td>10.1</td>
          <td>1.61</td>
          <td>7.69</td>
          <td>5.10</td>
          <td>105.0</td>
          <td>1016.4</td>
          <td>10.1</td>
          <td>14.0</td>
          <td>9.5</td>
          <td>99.0</td>
          <td>99.0</td>
          <td>1.203866</td>
          <td>1.253801</td>
          <td>96.017309</td>
        </tr>
        <tr>
          <td>2020-02-01 02:50:00</td>
          <td>24</td>
          <td>7.8</td>
          <td>9.1</td>
          <td>1.68</td>
          <td>7.14</td>
          <td>5.17</td>
          <td>103.0</td>
          <td>1015.8</td>
          <td>10.1</td>
          <td>14.1</td>
          <td>9.5</td>
          <td>99.0</td>
          <td>99.0</td>
          <td>1.203866</td>
          <td>1.253801</td>
          <td>96.017309</td>
        </tr>
        <tr>
          <td>2020-02-01 03:50:00</td>
          <td>35</td>
          <td>7.1</td>
          <td>9.0</td>
          <td>1.59</td>
          <td>6.67</td>
          <td>5.13</td>
          <td>103.0</td>
          <td>1015.6</td>
          <td>10.6</td>
          <td>14.2</td>
          <td>10.0</td>
          <td>99.0</td>
          <td>99.0</td>
          <td>1.245352</td>
          <td>1.296822</td>
          <td>96.031065</td>
        </tr>
      </tbody>
    </table>
    </div>

.. raw:: html

   <br>

You may only assign a :obj:`pandas.DataFrame` to :attr:`.Aero.df`,

    >>> # this will not work, df needs to be a dataframe
    >>> Aero_empty.df = 'high five'

::

    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    <ipython-input-13-5de371e56275> in <module>
          1 # this will not work, df needs to be a dataframe
    ----> 2 Aero_empty.df = 'high five'
    

    ~/AeroEvap/aeroevap/aero.py in df(self, df)
        122     def df(self, df):
        123         if not isinstance(df, pd.DataFrame):
    --> 124             raise TypeError("Must assign a pandas.DataFrame object")
        125         self._df = df
        126 


    TypeError: Must assign a pandas.DataFrame object


.. Tip:: 
   The ``df`` is a property of the :obj:`.Aero` class which means it can be
   assigned or reassigned if, for example, you wanted to run the evaporation
   calculations on a modified version of input meterological time series
   without creating a new :obj:`.Aero` instance.

Input variables and units
-------------------------

The meterological variables needed for running the aerodynamic
mass-transfer estimation of evaporation are the following:

================= ===== ======
variable          units naming
================= ===== ======
wind speed        m/s   WS
air pressure      mbar  P
air temperature   C     T_air
skin temperature  C     T_skin
relative humidity 0-100 RH
================= ===== ======

where the “naming” column refers to the internal names expected by the
:meth:`.Aero.run` method, i.e. the column headers in the dataframe should
either be named accordingly or a dictionary that maps your column names
to those internal names can be passed (see examples below).

To run the evaporation calculation you will also need the anemometer
height in meters and the temporal sampling frequency of the data in
seconds.

Run calculation on time series
------------------------------

As mentioned, this dataset has unique naming conventions, therefore we need to tell ``AeroEvap`` which variables are which with a dictionary,

    >>> # make a naming dict to match up columns with Aero variable names
    >>> names = {
    >>>     'WSPD' : 'WS',
    >>>     'ATMP' : 'T_air',
    >>>     'WTMP' : 'T_skin',
    >>>     'PRES' : 'P'
    >>> }

Alternatively you could rename wind speed, air and surface temperature, and air pressure columns to the apprpriate names specified in the table above in :ref:`Input variables and units`.

Now we are ready to run the aerodynamic mass-transer evaporation on the full
time series in our dataframe. Lastly, the sensor height of the anemometer and
temporal sampling frequency of the data needs to be supplied, in this case the
height is 4 meters and the data frequency is 10 minutes or 600 seconds.

This example assumes there are 8 physical or logical processors
available for parallelization, if not specified the :meth:`.Aero.run` routine
will attempt to use half of the available processors.

    >>> np.seterr('ignore')
    >>> # create a new Aero object and calculate evaporation on all rows
    >>> A = Aero(met_df)
    >>> A.run(sensor_height=4, timestep=600, variable_names=names)

After the calculations are complete three new time series will be added to the
:attr:`.Aero.df` dataframe: ‘E’, ‘Ce’, and ‘VPD’ which are open-water evaporation (mm/timestep), bulk transfer coefficient, and vapor pressure deficit
(kPa).

    >>> A.df[['E', 'Ce', 'VPD']].head()

.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>E</th>
          <th>Ce</th>
          <th>VPD</th>
        </tr>
        <tr>
          <th>date</th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>2020-01-31 23:50:00</td>
          <td>0.023573</td>
          <td>0.001552</td>
          <td>0.327503</td>
        </tr>
        <tr>
          <td>2020-02-01 00:50:00</td>
          <td>0.020519</td>
          <td>0.001527</td>
          <td>0.360776</td>
        </tr>
        <tr>
          <td>2020-02-01 01:50:00</td>
          <td>0.024901</td>
          <td>0.001545</td>
          <td>0.411624</td>
        </tr>
        <tr>
          <td>2020-02-01 02:50:00</td>
          <td>0.023312</td>
          <td>0.001538</td>
          <td>0.422028</td>
        </tr>
        <tr>
          <td>2020-02-01 03:50:00</td>
          <td>0.019437</td>
          <td>0.001519</td>
          <td>0.391987</td>
        </tr>
      </tbody>
    </table>
    </div>


View the calculated evaporation,

    >>> plt.figure(figsize=(8,4))
    >>> A.df.E.plot()
    >>> plt.ylabel('evaporation mm/10 min')


.. figure:: _static/evap_10min.png


The calculated open-water evaporation is shown below after creating a
daily sum.

    >>> plt.figure(figsize=(8,4))
    >>> A.df.E.resample('D').sum().plot()
    >>> plt.ylabel('evaporation mm/day')


.. figure:: _static/evap_daily.png

And the wind speed relation versus the calculated evaporation.


    >>> plt.figure(figsize=(8,4))
    >>> plt.scatter(A.df.WSPD.resample('D').mean(), A.df.E.resample('D').sum())
    >>> plt.ylabel('evaporation mm/day')
    >>> plt.xlabel('mean daily wind speed m/s')


.. figure:: _static/wind_vs_evap.png


Single calculation
------------------

The :obj:`.Aero` class also provides a method :meth:`.Aero.single_calc` that can
be used on a single set of meterological data to calculate the
instantaneous open-water evaporation. It requires the same inputs as
:meth:`.Aero.run` however the inputs are scalars as opposed to time series.
For example using the first timestamp of our example buoy data we can
calculate E, Ce, and VPD:

    >>> datetime = '2019-08-01 00:00:00'
    >>> wind = 3.3
    >>> pressure = 1021.2
    >>> T_air = 18.1
    >>> T_skin = 18.4
    >>> RH = 80.26
    >>> sensor_height = 4
    >>> timestep = 600
    >>> E, Ce, VPD = Aero.single_calc(
    >>>     datetime,
    >>>     wind,
    >>>     pressure,
    >>>     T_air,
    >>>     T_skin,
    >>>     RH,
    >>>     sensor_height,
    >>>     timestep
    >>> )

    >>> E, Ce, VPD
        (0.008724959939647368, 0.001310850807452679, 0.44947250457458576)

Theory behind calculations
--------------------------

This is a work in progress, for now please refer to `references hosted on GitHub <https://github.com/WSWUP/AeroEvap/tree/master/references>`_ about the methodologies used.
