Running chronosta.py

Ensure that Python 3, is installed along with the Pyhton and R packages indicated in the installation instructions.

Unix and macOS

  1.- In a terminal, create a directory (e.g., chronosta_example) and copy the partial timetrees in Newick format with the extention .nwk, along with the chronosta.py script:
    
    mkdir chronosta_example
    cp *.nwk /Users/barba/chronosta_example
    cp chronosta.py /Users/barba/chronosta_example
  
  2.- Type cd followed by the path to the directory to change the working directory to the folder that contains the chronosta.py script:
    
    cd /Users/barba/chronosta_example
  
  3.- To execute the Chrono-STA analysis type:
    python3 chronosta.py

Windows

  1.- In a command prompt, create a directory (e.g., chronosta_example) and copy the partial timetrees in Newick format with the extention .nwk, along with the chronosta.py script:

    mkdir chronosta_example
    copy *.nwk C:\Users\barba\chronosta_example
    copy chronosta.py C:\Users\barba\chronosta_example
  
  2.- Type cd followed by the path to the directory to change the working directory to the folder that contains the chronosta.py script:

    cd C:\Users\barba\chronosta_example
  
  3.- To execute the Chrono-STA analysis type:

    python3 chronosta.py

For questions or comments please contact:

jbarba@amnh.org
