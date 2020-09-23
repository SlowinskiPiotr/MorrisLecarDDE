# MorrisLecarDDE

Repository with code used to for numerical bifurcation analysis and simulations presented in: "Robust spike timing in a self-coupled excitable
cell with delay". Preprint availible at ...

The analysis was performed using the DDE-BIFTOOL v3.1.1 package (http://ddebiftool.sourceforge.net) in Matlab R2020a on macOS 10.14.6.

The code assumes familiarity with Matlab and the DDE-BIFTOOL package.

Any questions/ comments/ bugs please get in touch at p.m.slowinski@exeter.ac.uk

## Files and folders
* simulations.m - matlab script with a code to simulate the Morris–Lecar model with self coupling and delay (Figure 1 E and G in the paper)
* bifurcation_diagrams.m - matlab script with a code to perfrom numerical bifurcation analysis of the Morris–Lecar model with self coupling and delay (Figure 1 A-D and F in the paper)
* br_plot3.m - modified function br_plot.m from DDE-BIFTOOL v3.1.1 to allow 3D plots of the branches
* br_splot3.m - modified function br_splot.m from DDE-BIFTOOL v3.1.1 to allow 3D plots of the branches with stability information
* _System_Files_ - folder with matlab files that define the Morris–Lecar model with self coupling and delay (the files are required by the simulations.m and bifurcation_diagrams.m scripts)
* _Data_ - folder with data used in Figure 3 and Figures S1-S3
* _DynamicClamp_ - folder with code used by the dynamic clamp
