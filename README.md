# Assignment - 9 
# Particle Methods for Fluid Flow Simulation

Run the bash script run.sh to generate report.pdf which contains the results by either using the command bash run.sh or ./run.sh
If the command results in any ownership errors then use chmod ug+x run.sh to give yourselves executing permission.

The python code will take very long time to run and even after using pypy the results were not improved so keep using python.

More error might be generated if required LATEX packages are not present. In that case the pictures of the pressure, density, energy and
velocity can be directly seen in the directory. report.pdf is also included from which the results generated can be seen without running the code

Time is printed in the python script for the user to see the progress of the program, remove this is not required in sod.py line 116. The total time of the program is 0.2 s and time step is 0.0001
