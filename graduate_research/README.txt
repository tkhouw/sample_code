This directory contains a sample of the code I've written for my graduate research in statistical mechanics.
For this project, I simulate an agent-based economic model called the Growth, Exchange, and Distribution (GED) model.

Simulations are run by main.cpp using structs defined in exchange.h and exchange.cpp.
Running multiple simulations in parallel with different parameters is made easy by calling ./bash.sh in terminal.
bash.sh passes parameters into read.in which is read by main.cpp.

After simulations are complete, Python scripts in analysis/ handle the data analysis.  I've included only a few example scripts.
