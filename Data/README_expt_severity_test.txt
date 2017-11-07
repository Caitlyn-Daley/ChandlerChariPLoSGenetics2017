
Some information on data collection for GS_Wings_RawArea_Oct2015.csv
Data Provenance (how data was collected):
via dissections, area imaging
Author: 
Performed - Ian Dworkin, Gayatri Sivaratnam, Samiksha Kaul 
Dissections - Gayatri Sivaratnam, Samiksha Kaul
Imaging - I. Dworkin, G. Sivaratnam, Samiksha Kaul 

Dates 
samples Collected: July 21, 2015
Dissection: August 17 - September25, 2015
Imging: September 30, 2015 (ID, GS, SK)
(Date samples were collected, dates for imaging)

A bit about the experiment:
16 DGRP lines (males) were crossed to 3 mutations (sd 1, sd E3, sd 58d) and a wildtype (SAM) (females). 
The effects of intermediate mutation on different lines, whether the severity of this effect varied across the different lines is being investigated. Only hemizygous (for sd alleles) males were phenotyped.

Microscope, camera used. Leica m125 with a Leica DFC400 camera

Magnification : 50X (0.5X?)
aperture about halfway across (59%)
Brightness: 68%
Gamma: 0.80

Notes:
naming: SAM/sd1/sde3/sd58d_line#_M_1/2/3/4/5/6.TIF
Imaging from slides
1	4
2	5
3	6

Septembr 30 - GS, SK, ID
sdE3X313 #5 ripped - no image acquired
sd1x313 #1 unknown if ripped - image acquired
sd1x313 #3 scratched coverslip - image acquired
sd1x313 #4 unknown thread under coverslip - image acquired
sde3x341 #1 and #2 , #4 folded, #6 is torn - no image acquired
sd58dx341 #2, #3 and #6 air bubble very close to wing edge 
sd1x341 #5 has small wing tear in the middle
sd58dx379 #6 - bubble surrounding - image acquired
sde3x399 #1 - potential rip, confirm with ID
sde3x399 #2 - potential rip/fold, ask ID - ripped, don't use
sde3x399 #5 - potential rip, ask ID 
left off at sde3x399

October 1 - SK
Started at sd58dx399
sd58dx399 #6 - Possible foreign particles on the wing - image acquired  
sd1x399 #6 - foreign particle on the wing - image acquired 
SAMx365 #5 - particle blocking part of wing outline - image acquired 
sde3x365 #1 - unknown thread under coverslip - image acquired 
sde3x365 #2 - air bubble covering left side of the wing - image acquired 
sde3x365 #3 - check with ID for ripped wing - image acquired 
sd58dx365 #2 - check wtih ID for folded wing - no image acquired
sd58dx365 #6 - damaged/folded wing - no image acquired
sde3x362 #3 - ripped wing, check with ID - no image acquired
sde3x362 #4 - bottom wing edges folded - no image acquired 
sd58dx362 #3, #4 - folded or damaged?, check with ID 
left off at sde3x315 

October 2 - SK
Started at sd58dx315
sd58dx315 #2 - ripped wing - no image acquired 
sd58dx315 #4 - foldd wing - no image acquired 
sd1x315 #3 - wing border unclear - no image acquired 
sde3x380 #1 - ripped and folded wing - no image acquired 
sde3x380 #4 - wing border folded - no image acquired 
sde3x380 #6 - ripped wing - no image acquired 
sd1x375 #6 - foreign object on the wing - image acquired 
sd1x357 #6 - foreign object on wing - image acquired 
sd58d #3 - wing folded slightly, check with ID - image acquired 
*Two sd1x375 slides so only one imaged 
SAMx380 #3 - coverslip scratched 
left off at SAMx360

October 5 - GS
Started with sd58dx360
sd58dx360 #5 - bubble impeding edge - image acquired
sd1x360 #4 - small fold, ask ID - image acquired
sde3x357 #2 - potential tear in wing, refer to #3 (similar), ask ID - image acquired
sde3x357 #6 - wing too big for screen, potentially female? - no image acquired
sd58dx357 #1 - wing torn - no image acquired
sd58dx357 #5 - small fold, potential vein bubble - image aquired
sd1x357 #1 - alar lobe folded - image acquired

2 slides of sd1x357 - used slide made September 17, 2015
2 slides of sd1x375 - used slide made September 14, 2015


sd1x357 #4 - folded alar lobe - image acquired
samX358 #6 - tip was torn - image acquired
sde3x358 #6 - wing partially folded, ask ID - image acquired
sd1x358 #4 - wing torn - image acquired
finished with sd1x358, next imaging should begin with SAMX324

October 9 - GS
started with SAMX324

2 Slides of sd1x324 - used slide made Spetember 17, 2015

sd1x324 #4 - folded alar lobe - image acquired
sde3x304 #1 - female wing - no image acquired
sde3x324 #2 - wing tip torn - image acquired
sd1x304 #2 - slight rip - image acquired

dissection imaging complete - phase 1

October 19,2015
Slides to Re-do

sde3x341 - no more males to dissect
sde3x362 - x
sde3x399 - done
sd58dx365 - done
sde3x362 - x
sd58dx362 - x
sd58dx315 - done
sde3x380 - done
sde3x357 - done using 18/8/15

362 crosses were discarded due to experimental errors

	7	10	
13	8	11
14	9	12

extra dissections completed

area recording next



Note on extracting data
Used custom macro in ImageJ v1.50c  the code is below
// Macro to compute wing area for Drosophila images.

macro "WingSizeMacro2 [q]" {  
run("Find Edges");
run("Enhance Contrast...", "saturated=0.4");
run("Sharpen");
run("Sharpen");
run("Make Binary");
run("Close-");
run("Dilate");
run("Fill Holes");
run("Despeckle");
run("Remove Outliers...", "radius=50 threshold=50 which=Dark");
run("Analyze Particles...", "size=50-Infinity pixel display summarize");
// run("Close");
//
