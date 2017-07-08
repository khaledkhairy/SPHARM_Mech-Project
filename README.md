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
------------------------------------------------------------------------------

## Status: 
In production use at Janelia. This is a nascent set of tools that is undergoing large changes and code cleanup. We consider the library suitable for use by our collaborators as well as other research groups. Due to limited staffing, we do not guarantee support for outside groups.

## Summary

SPHARM-MECH is an application aimed at exposing functionality in the general shape_tools library, and making it easier for the user to try different configurations and parameters. SPHARM-MECH is taylored towards the specific problem of modeling tissue mechanics in biological systems. This application is a "sandbox" and not all features have proven useful for modeling tissue mechanics.

To compile SPHARM-MECH please use cmake.

Modify lines 3, 21, 22 and 59 in CMakeLists.txt.

----------------------------------------------------------------------------------
For MacOSX Sierra, it is easier to just use the supplied binary
----------------------------------------------------------------------------------
SPHARM-MECH was compiled using Qt 5.7 and VTK 7.1.1 (although other versions might also work)

Building notes:

The application has been built using the following C/C++ libraries on both MacOSX and Windows 7:

[1] Eigen 3.1.1: http://eigen.tuxfamily.org

[2] VTK 7.1.1: http://www.vtk.org/

[3] Qt 5.7: http://www.qt.io/

[4] shape tools: https://github.com/khaledkhairy/shape_tools

Building has not been tested extensively on any other configuration but should be straightforward.


After generating a minimum energy shape, to view the tissue flow field please download put the Matlab folder in your path and execute commands:

'''json
fn1 = '/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/Fly_embryo_Berkeley_stage5_last_Lmax80_with_dorsal_L56.shp3';
fn2 = '/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/vfi_result.shp3';

s1 = shp_surface;
s1 = s1.read_shp_surface_ascii(fn1);
m1 = get_mesh(s1);
m1 = m1.translate_to_center_of_mass();

s2 = shp_surface;
s2 = s2.read_shp_surface_ascii(fn2);
m2 = get_mesh(s2);
m2 = m2.translate_to_center_of_mass();
plot_difference(m1,m2);
'''
