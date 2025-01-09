Running chronosta.py

Ensure that Python 3, is installed along with the Pyhton and R packages indicated in the installation instructions.

**Unix and macOS**

1.	In a terminal, create a directory (e.g., `chronosta_analysis`) and copy all the constituent timetrees in the `simulated_timetree_collection directory` (in Newick format with the `.nwk` extension), as well as the `chronosta.py` script located in the `code` directory (suppose you have extracted the archive into `~/chrono-sta`)

    mkdir -p ~/chronosta_analysis
    cp ~/chrono-sta/examples/simulated_timetree_collection/*.nwk ~/chronosta_analysis
    cp ~/chrono-sta/code/chronosta.py ~/chronosta_analysis

2.	Type `cd` followed by the path to the created directory to change the working directory to the folder that contains the `chronosta.py` script and the constituent timetrees:	

    cd ~/chronosta_analysis

3.	To execute the Chrono-STA analysis, type:

    chronosta.py

If `-bash: ./chronosta.py: Permission denied` is displayed, make the `chronosta.py` script executable by typing `chmod +x chronosta.py`.

**Windows**

1.	In a command prompt, create a directory (e.g., `chronosta_analysis`) and copy all the constituent timetrees in the `simulated_timetree_collection` directory (in Newick format with the `.nwk` extension), as well as the `chronosta.py` script located in the `code` directory (suppose you have extracted the archive into `C:\Programs\chrono-sta\`):

    mkdir C:\chronosta_analysis
    copy C:\Programs\chrono-sta\examples\simulated_timetree_collection\*.nwk C:\chronosta_analysis
    copy C:\Programs\chrono-sta\code\chronosta.py C:\chronosta_analysis


2.	Type `cd` followed by the path to the created directory to change the working directory to the folder that contains the `chronosta.py` script and the constituent timetrees:	

    cd C:\chronosta_analysis


3.	To execute the Chrono-STA analysis, type:

    python chronosta.py

For questions or comments please contact:

jbarba@amnh.org
