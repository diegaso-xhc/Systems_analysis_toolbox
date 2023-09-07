# Systems analysis toolbox

## Overview of the repostitory
<div align="justify">
Simulating systems response is crucial to understand them and to develop suitable controllers depending on application requirements. For instance in 5-axis CNCs, the control parameters used for each one of their axes could have a considerable impact on the overall manufacturing resolution. In certain robotics applications, controllers need to operate with high precision. 
<br />
<br />
Several high performance platforms (e.g. Matlab) offer a control systems toolbox which can be used reliably to model and control systems. Nevertheless, these toolboxes tend to be built on proprietary functions, which can't be customized depending on user needs. This repository is intended to provide an open source platform with similar functionalities for representing systems responses and understanding what controllers are suitable for them. Additionally, thanks to its versatility and structure, this repository can be used for teaching purposes, helping students understand the magice happening behing control toolboxes.
<br /> 
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

The project has been compiled using the aforementioned libraries into a self contained project. Nevertheless, if some compilation errors arise, please check the aforementioned versions. The most relevant files on the repository are detailed as follows:

```
- backupFunctions.cs --> Operates on parameters within system modeling, such as
raw phase and angle vectors extracted from a system's response.
- DifferentialEquations.cs --> Handles system model representations (e.g. from
transfer functions to state space representations, or to differential equations). It also contains Runge–Kutta methods
for finding approximate solutions of nonlinear equations. Additionally, it contains functions for transfer functions
algorithmics (e.g. series, or feedback loops). Finally, it contains print functions specifically tailored to each one
of the system model representations (e.g. state space).
- fastFourierTransform.cs --> Transforms time space vectors into the frequency
(fast Fourier transform) and viceversa (inverse fast Fourier transform).
- Matrix.cs --> Handles matrices in an efficient manner.
- MainWindow.xaml --> Contains the code required by WPF to launch the GUI for the user.
- MainWindow.xaml.cs --> Handles user requests and returns required outputs.
```
<br />
UPDATE (09.2023): Due to a recent change of the chart visualization library, it is worth noticing at the moment only the time response and frequency spectrum can be displayed graphically. The authors are currently updating the Bode and Nichols diagrams. Nevertheless, besides visualization, all functions are implemented and fully functional. The user simply needs to print the outputs or create methods to use them.
<br />
<br />

## Contributions

The contributions of this repository can be summarized as follows:

```
- Algorithms for handling differential equations, transfer functions, and state space representations
- Algorithms to manipulate matrix in an intuitive and efficient manner
- Algorithms to operate on system response parameters (e.g. angle, phase, etc.)
- An intuitive GUI (with similar nomenclature as MATLAB for transfer functions) to check systems time and frequency response
- An open source code, which can be used for teaching of fundamentals of control systems.
```

## Examples of GUI usage

### Visualization of human hand-object interaction

Plots of the contact interaction of the subject in a specified frame. Contact surfaces are visible in read. 

<p align="center">
   <img src="/Visualizations/Sine_response.PNG" width="650" />
</p>

### Grasping of a cylinder

<p align="center">
   <img src="/Visualizations/Frequency_spectrum.PNG" width="400" />
</p>

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
