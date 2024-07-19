# BIGf90 0.3.1

* Users now can have executable files, input files, and output files in different locations
* Warnings were added in case of run_gibbs and run_postgibbs functions that requires inputs to be together with renum results
* Code changed to always use global path even if it is referring to working directory
* Verbose argument added to give option to print or not information on screen while running
* Default values added
* Output directories need to exist. The function will stop if they don´t
* Automatic installation checks included
* NEWS file included