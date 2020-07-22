# Swimmer-Spatial-Tracker-with-Phase-Contrast-Imaging
The tracking algorithm is base on <Taute, K., Gude, S., Tans, S. et al. High-throughput 3D tracking of bacteria on a standard phase contrast microscope. Nat Commun 6, 8776 (2015)>

## Protocol
First, Using <SwimmerSpatialTracker_Crop.m> to crop suiteble cells as reference cells.\
Next, Using <SwimmerSpatialTracker_LibBuilder.m> to align reference cells with maximal central intensity and using rotating operation to compensate optical path tilt.\
After build library, you can start to track your cell by <SwimmerSpatialTracker.m>
## Example 1: E. coli strain RP437
## Example 2: V. alginolyticus strain VIO5
## Example 3: V. fischeri strain ATCC7744
