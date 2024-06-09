README
======

Data package to allow replication of the results in the Bourne & Oates article: 
"Extreme threshold failures within a heterogeneous elastic thin-sheet and the 
spatial-temporal development of induced seismicity within the Groningen gas field"


ReservoirModel_24bcm_pressure.csv
=================================
Comma separated values of instantaneous mean reservoir pressure [bar] 
within regular reservoir grid blocks that span the reservoir thickness
and an annual time sampling from 1st January 1958 to 1st January 2022.
Grid center coordinates [m] are provided in the Netherlands National
Triangulation Coordinate System (Rijksdriehoek). 

The period from 1958 to 2017 is history-matched against reservoir gas 
production and reservoir pressure measurements. The period from 2017 to 
2022 represent a forecast based on 24 billion cubic meters (bcm) per year
of total gas production from the Groningen field.

ReservoirModel_24bcm_compaction_linear.csv
==========================================
Comma separated values of instantaneous mean reservoir compaction [m] 
within regular reservoir grid blocks in the same grid format as reservoir 
pressures. Compaction is history matched to surface geodesy observations
between 1960 and 2017. Outside this range reservoir compaction is forecast
using a linear relationship between reservoir pressure and compaction.

ReservoirModel_thickness.csv
============================
Comma separated values of reservoir thickness [m] within regular reservoir 
grid blocks in the same grid format as reservoir pressures. These are based
on direct measurements via well penetrations and indirect measurements from
interpretation of a reflection seismic image of the reservoir.

FaultModel_geometries.csv
=========================
Comma separated values of geometric fault properties for all faults interpretated
from the reflection seismic image at the top Top Rotliegend horizon. 
Grid coordinates [m] are provided in the Netherlands National
Triangulation Coordinate System (Rijksdriehoek). The geometric attributes include
elevation [m], reservoir thickness [m] measured as net (1) and gross (2), 
dip azimuth [degrees] and dip [degrees].

Horizon_Top_Rotliegend_50x50.txt
================================
Space separated values of depth below datum [m] of the Top Rotliegend horizon 
within a regular 50 m by 50 m grid obtained from interpretation of a reflection 
seismic image of the reservoir.

EventCatalogue_20161231.csv
===========================
Extract from the KNMI earthquake catalogue for events within the Groningen field area.
Comma separated values characterize every report event with the following attributes:
Place [string], Date [DD/MM/YYYY], Time [hh:mm:ss], Magnitude [local magnitude], Depth [km], 
Easting [m], Northing [m].






