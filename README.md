# Systems analysis toolbox

## Overview of the repostitory
<div align="justify">
Simulating systems response is crucial to understand them and to develop suitable controllers depending on application requirements. For instance in 5-axis CNCs, the control parameters used for each one of their axes could have a considerable impact on the overall manufacturing resolution. In certain robotics applications, controllers need to operate with high precision. 
<br />
Several high performance platforms (e.g. Matlab) offer a control systems toolbox which can be used reliably to model and control systems. Nevertheless, these toolboxes tend to be built on proprietary functions, which can't be customized depending on user needs. This repository is intended to provide an open source platform with similar functionalities for representing systems responses and understanding what controllers are suitable for them. Additionally, thanks to its versatility and structure, this repository can be used for teaching purposes, helping students understand the magice happening behing control toolboxes.
<br /> 
This toolbox was used on for the development of the <a href="https://www.sciencedirect.com/science/article/pii/S2212827117302305">scientific article</a>, which shows an approach to reduce inaccuracies on five-axis CNCs.
<br /> 
<br /> 
 <p align="center">
   <img src="/Visualizations/Time_response_cursor.png" width="700" />
</p>
<br />

## Understanding repository

The repository was developed in C#, using the following software version:
```
- Microsoft Visual Studio Community 2019 (Version 16.11.29)
- Target Framework: .NET Framework 4.6.1
- OxyPlot.Wpf (Version 2.1.2) --> Can be installed using NuGet package manager
- OxyPlot (Version 2.1.2) --> Can be installed using NuGet package manager
- OxyPlot.Wpf.Shared (Version 2.1.2) --> Can be installed using NuGet package manager
```

The project has been compiled using the aforementioned libraries into a self contained project. Nevertheless, if some compilation errors arise, please check the aforementioned versions.
<br />
The files on the repository are detailed as follows:

```
- Microsoft Visual Studio Community 2019 (Version 16.11.29)
- Target Framework: .NET Framework 4.6.1
- OxyPlot.Wpf (Version 2.1.2) --> Can be installed using NuGet package manager
- OxyPlot (Version 2.1.2) --> Can be installed using NuGet package manager
- OxyPlot.Wpf.Shared (Version 2.1.2) --> Can be installed using NuGet package manager
```





## Experimental setup for repository

Although this repository is easily extendable, it is worth noting that during the development of this code a given fixed setup was used. The setup consisted of:

```
- VICON motion tracking capture with 16 IR-cameras
- 12 purposely selected 3D-printed (PLA) objects
- Artec Eva Lite 3D - Scanner
- 26 4-mm marker hand setup (see image below)
```

<p align="center">
   <img src="/Visualizations/Frequency_spectrum.PNG" width="400" />
</p>
 
<!---
your comment goes here
and here
![This is an image](/Visualizations/vis_1.png)
![This is an image](/Visualizations/grasp_wine_glass.gif) ![This is an image](/Visualizations/grasp_cup.gif)
<img src = "/Visualizations/grasp_wine_glass.gif" width="400"> <img src = "/Visualizations/grasp_cup.gif" width="570">

-->
</div>

Nevertheless, the contributions of our work provide tools to extend our results to different setups, provided that changes are made on the pertinent locations. If you are not sure how to do so, please contact the author (see details below).

## Contributions

The contributions of this repository can be summarized as follows:

```
- A library for post processing and visualization of motion capture data --> HMCL_lib
- A set of purposely selected objects, which guarantee grasping postures variability --> Info_objects
- Algorithms to track known 3D objects with random marker placements --> Track_object_lib
- Algorithms for mesh manipulation --> Mesh_Manip_lib
- Pseudonymized data from trials (saved as Study objects, see HMCL_lib) --> Trials
- Algorithms to calculate contact surfaces information betwen a hand mesh and an object mesh --> see main
- Python algorithms to handle the MANO model, including marker position optimization --> fit_MANO (Developed in 
  cooperation with Omid Taheri (https://is.mpg.de/person/otaheri).
```

## Examples of hand-object contact level human manipulation

### Visualization of human hand-object interaction

Plots of the contact interaction of the subject in a specified frame. Contact surfaces are visible in read. 

<p align="center">
   <img src="/Visualizations/Sine_response.PNG" width="650" />
</p>

### Grasping of a cylinder

Plots of object with contact surfaces for grasping analysis. Contact surfaces are visible in red.

<p align="center">
  <img src="/Visualizations/grasp_cylinder.png" width="550" />  
</p>

### Grasping of a wine glass and a cup

Digital representation of the manipulation of a wine glass and a cup.

<p align="center">
  <img src="/Visualizations/grasp_wine_glass.gif" width="342" />
  <img src="/Visualizations/grasp_cup.gif" width="575" />   
</p>

### Contact surfaces with different hand parameters

A generic and a customized hand models were used to calculate the contact surfaces during the human manipulation recording. It can be seen that the hand model parameters play a role on the accuracy of the contact surfaces.

<p align="center">
  <img src="/Visualizations/Contact_surface_diff.PNG" width="600" />  
</p>

## License

Developed by Diego Hidalgo C. (2021). This repository is intended for research purposes only. If you wish to use any parts of the provided code for commercial purposes, please contact the author at hidalgocdiego@gmail.com.
