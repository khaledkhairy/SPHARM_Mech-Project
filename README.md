# SPHARM_Mech-Project
Spherical harmonics-assisted modeling of tissue mechanics.
![Alt text](https://github.com/khaledkhairy/SPHARM_Mech-Project/blob/master/clips/Screen%20Shot%202017-07-06%20at%209.25.56%20AM.png "SPHARM-MECH screenshot")

-----------------------------------------------------------------------------
Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  
* Neither the name HHMI-Janelia Research Campus nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL HHMI-JANELIA RESEARCH CAMPUS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



## Status: 
In production use at Janelia. This is a nascent set of tools that is undergoing large changes and code cleanup. We consider the library suitable for use by our collaborators as well as other research groups. Due to limited staffing, we do not guarantee support for outside groups.

## Summary

SPHARM-MECH is an application aimed at exposing functionality in the general shape_tools library, and making it easier for the user to try different configurations and parameters. SPHARM-MECH is taylored towards the specific problem of modeling tissue mechanics in biological systems. This application is a "sandbox" and not all features have proven useful for modeling tissue mechanics.

Building notes:
------------------------------------------------------------------------------
To compile SPHARM-MECH please use cmake.

Modify lines 3, 21, 22 and 59 in CMakeLists.txt.

For MacOSX Sierra, it is easier to just use the supplied binary

SPHARM-MECH was compiled using Qt 5.7 and VTK 7.1.1 (although other versions might also work)


The application has been built using the following C/C++ libraries on both MacOSX and Windows 7:

[1] Eigen 3.1.1: http://eigen.tuxfamily.org

[2] VTK 7.1.1: http://www.vtk.org/

[3] Qt 5.7: http://www.qt.io/

[4] shape tools: https://github.com/khaledkhairy/shape_tools


Building has not been tested extensively on any other configuration but should be straightforward.

Basic modeling example steps:
----------------------------------------------------
(Note before you start: At step [5] you will running an optimization that can take some time. During this optimzation, it is possible that the shape disappears from the display window. Press <R> in this case to re-center the shape.)
[1] Launch the application.

[2] Click on "Load undeformed surface" OR drag a surface (shp3) file onto the application window. The example fly embryo provided, with gene expression patterns is <path to repo>/MACOSXBinary/test/Fly_embryo_Berkeley_stage5_last_Lmax80_with_dorsal_L56.shp3. The shape loads and shows up in the upper left window. You can rotate and zoom in/out with the mouse. You can also click on "curvature" to see the shape colored according to local mean curvature.

[3] Select the gene expression scalar field dorsal_s from the drop-down menu under "Fold 1" and then enter 1.0 in the edit field next to set that the dorsal_s expression will be used 100% (there are no other scalar fields for this run). The configuration panel should look similar to this:
![Alt text](https://github.com/khaledkhairy/SPHARM_Mech-Project/blob/master/clips/example_configuration.jpg "example configuration")

[4] In the optimization panel make sure the SUBPLEX(NLOPT) tab is selected, set the number of iterations to 500,000 and press "start" to start the optimization.

[5] The shape updates during optimization and current energy is plotted in the lower left panel. The text feedback provides information about <total shape energy> <volume constraint energy> <self-intersection (should be 0)> <self-check of curvature calculation using total Gaussian curvature (should be close to 1.0 throughout)>. Quirk:Should the shape disappear during optimization, please click anywhere into the upper left window and press R.

[6] At the end of optimization, or when the user presses "stop" the current lowest energy shape will be shown. Here is an example when stopping processing after 3 hours of step [4].

![Alt text](https://github.com/khaledkhairy/SPHARM_Mech-Project/blob/master/clips/Screen_after_3h.png "SPHARM-MECH screenshot")

The result above is saved in: <path to repo>/MACOSXBinary/test/Result_after_3_h_Fly_embryo_Berkeley_stage5_last_Lmax80_with_dorsal_L56.shp3

[7] Export/save the shape from the "File" menu as the common obj format for viewing/importing into other applications or shp3 format which can be used with SHAPE, SPHARM-MECH or the Matlab tools accompanying this work.

If a morphogen gradient-based analysis is desired, please follow the simple instructions in the ![Alt text](https://github.com/khaledkhairy/shape_tools "shape_tools") repository for modification of the shell class. 

Flow field visualization 
------------------------------------------------------------------------------
It has not been implemented in the SPHARM-MECH GUI itself yet. 
If you have Matlab available, after generating a minimum energy shape, to view the tissue flow field please ensure that the "matlab_code" folder is in your path and in Matlab execute commands:

```json
fn1 = '<path to file>/Fly_embryo_Berkeley_stage5_last_Lmax80_with_dorsal_L56.shp3';
fn2 = '<path to file>/vfi_result.shp3';

s1 = shp_surface;
s1 = s1.read_shp_surface_ascii(fn1);
m1 = get_mesh(s1);
m1 = m1.translate_to_center_of_mass();

s2 = shp_surface;
s2 = s2.read_shp_surface_ascii(fn2);
m2 = get_mesh(s2);
m2 = m2.translate_to_center_of_mass();
plot_difference(m1,m2);
```
