Running chronosta.py

Ensure that Python 3, is installed along with the Pyhton and R packages indicated in the installation instructions.

Unix and macOS

1. In a terminal, create a directory (e.g., 'simulated_timetree_collection') and copy the constituent timetrees in Newick format with the extension .nwk, along with the 'chronosta.py' script located inside the 'code' directory:

    mkdir /Users/barba/chrono-sta/examples/simulated_timetree_collection
    cp *.nwk /Users/barba/chrono-sta/examples/simulated_timetree_collection
    cp chronosta.py /Users/barba/chrono-sta/examples/simulated_timetree_collection

2. Type 'cd' followed by the path to the created directory to change the working directory to the folder that contains the chronosta.py script:	

    cd /Users/barba/chrono-sta/examples/simulated_timetree_collection

3. To execute the Chrono-STA analysis, type:

    python3 chronosta.py

Windows

1. In a command prompt, create a directory (e.g., 'simulated_timetree_collection') and copy the constituent timetrees in Newick format with the extension .nwk, along with the 'chronosta.py' script located inside the 'code' directory:

    mkdir C:\Users\barba\chrono-sta\examples\simulated_timetree_collection
    copy *.nwk C:\Users\barba\chrono-sta\examples\simulated_timetree_collection
    copy chronosta.py C:\Users\chrono-sta\barba\examples\simulated_timetree_collection

2. Type 'cd' followed by the path to the created directory to change the working directory to the folder that contains the chronosta.py script:	

    cd C:\Users\barba\chrono-sta\examples\simulated_timetree_collection

3. To execute the Chrono-STA analysis, type:

    python3 chronosta.py

For questions or comments please contact:

jbarba@amnh.org
