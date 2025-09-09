# Modeling and simulating the dyanamics of a Skyhook

This project was part of Spacecraft Dynamics and Control course. The repo contains the report and all the codes used for the project.

https://github.com/user-attachments/assets/9591d3aa-4e19-46ff-b501-e89ff1fb0b1b


## Contents
- **Dumbell_sat_sim.py**

    1. This file contains all the function definitions for the simulations. When the file is run it plots the state variable time evolution.
  
    2. By enabling the save flags inside the file it also saves all the relavent data in .npz file and the plots as .pdf file.

- **Skyhook_vizualization.py**

     This file reads the time series data from the saved .npz file and shows a `3D vizualization` of the satellite in orbit.
     `VPython` library was used to make the above vizualization. The vizualization will not be upto scale.

- **Project_report.pdf**

    This file contains all the necessary information including the derivation for the simulation.
