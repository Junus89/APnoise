ADScase3 ...................................... Case name
1 ............................................. Solver ID [0: Barry and Magliozzi ;1: Hanson]
#================ PROPELLER AERO DATA================
2.5 ........................................... Propeller blade diameter (m)
5 ............................................. Number of blades
1900.0 ........................................ RPM
0.2 ........................................... Free stream Mach number
0.0 ........................................... pitch angle of propeller shaft axis relative to flight direction (deg)
10 ............................................ Acoustic harmonic numbers
BladeGeoData.txt .............................. Blade geometry data file
BladeAeroData.txt ............................. Blade aerodynamics data file
FD_analytical_blade_in_NB5RPM1900.txt ......... Input file for aerodynamic sources
#================= ATMOSPHERE SECTION ===============
2000.0 ........................................ Flight altitude (m)
275.15 ........................................ Temperature (K)
79500.0 ....................................... Pressure (Pa)
1.007 ......................................... Density (Kg/m^3)
332.65 ........................................ Sound speed at this atltitu (m/s)
#================= RECEIVER SECTION =================
1 ............................................. Logic for observer mics [0: automatically generated; 1: read from file]
ObsrvrGeomZ.txt ................................ Receiver geometry file
51 ............................................ Polar angle num
101 ........................................... Azimuthal angle num
#============ DEBUG & Dev SECTION ===================
DummyFile...................................... Added at line 12
