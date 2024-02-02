# User Guide
In this section, you will learn how to operate AutoEpiCollect and access the outputted files to view your final 
epitope results.
## Opening AutoEpiCollectGUI
First, navigate to the AutoEpiCollect directory using these command::

Windows:
```cmd
cd %HOMEPATH%
cd AutoEpiCollect
```
Mac and Linux:
```bash
cd ~
cd AutoEpiCollect
```
To open the GUI program, enter the following commands into your Terminal or Anaconda Prompt::
```bash
conda activate aec_venv
python AutoEpiCollectGUI.py
```
Please be prepared to wait up to 2 minutes for the GUI to open on its first run. You should see a screen 
similar to this one pop up.
![Labeled GUI for AutoEpiCollect](GUI-home-annotated.png)
## Inputting Required Data and Configuring Run Settings
The numbered boxes in the home screen above show the order of inputs:

1. These two fields are for specifying the gene target. The white box on the left is where you can type in the 
   name of your gene target. AutoEpiCollect will retrieve a FASTA-formatted file of this gene's peptide sequence 
   from UniProt. The button on the right is if you want to input your own FASTA-formatted gene sequence. Simply 
   click this button and choose your desired FASTA-formatted file. 
2. This field specifies how you want your epitope data to be collected. There are two options from the drop-down 
   menu, "Collect New Data" and "Update Existing Data". Choose the "Collect New Data" option if you want to start 
   the epitope selection process from scratch. Choose the "Update Existing Data" option if you have collected some 
   epitope data, but want to add to this data with available epitope characteristics not chosen in the past. If you 
   have selected the "Update Existing Data" option, please click the "Input Existing Class I (or Class II) Data" 
   button and choose your incomplete data file. Remember that this file must be the **unranked and unfiltered 
   epitope Excel file** for this option to work. A summary of all the initial and final filenames will be discussed 
   down below.
3. The third field is a drop-down menu asking which MHC Class you would like to collect epitope data from. You can 
   choose to collect Class I, Class II, or Class I and II data.
4. Field four is for inputting cancer types and point mutations of interest. Input a name or abbreviation of 
   each cancer type you are studying, as well as the relevant missense mutations associated with each cancer. Please 
   format the cancers and mutations exactly how you see in the example, with a new line separating each cancer.
5. Field five is where you specify the run configurations for the epitope selection process. Each epitope 
   characteristic is grouped into categories based on its properties. Binding affinity is obtained by default and 
   cannot be switched off. Click on the epitope characteristics that you want AutoEpiCollect to obtain when 
   selecting the top epitopes. The four-variable scoring function requires the immune characteristics in order to 
   work. Because of this, checking the four-variable scoring function box will automatically select the immune 
   characteristics. Filtering is based off of the epitope stability characteristics, toxicity, and IFN-gamma release.
   If filtering is checked, AutoEpiCollect will filter epitope data according to the filtering characteristics that 
   are checked. For example, if you want to filter epitope data and only toxicity and half-life are checked, then 
   filtering will take place based on these two characteristics. If instability is checked as well, then filtering 
   would be based off of all three characteristics. The last checkbox is for population coverage analysis. This setting 
   allows AutoEpiCollect to calculate the coverage of a vaccine composed of the top epitopes predicted by 
   AutoEpiCollect. This setting is optional, but is important when determining the utility of your finalized vaccine. 
6. The final button on the home screen labeled six is the submit button. Pressing this button will take you to the 
   output window shown below where you can view the progress of AutoEpiCollect's epitope selection method with live updates.

## Output Window and Viewing Final Epitope Results
![Output Window](output-window.png)
After pressing the submit button, you will be taken to this output window which shows live updates of 
AutoEpiCollect's progress. During this time, you will be allowed to hit the back button and view what options you 
have selected, but you cannot restart the epitope selection process until the program is complete or until you close 
out the program while it is running. Clicking on the "Click for Output" button when the program is completed will 
open an Excel file of all the pre-ranked and pre-filtered epitope data collected, organized by point mutation. View 
the table below to learn about what possible output files are made during AutoEpiCollect's epitope selection process.
You can find all the files if you navigate to the created AutoEpiCollect directory.

| File*                                            | Description             |
|--------------------------------------------------|-------------------------|
| {gene_target}.fasta                              | This FASTA-formatted fi |
| all_variables_{mhci or mhcii}.xlsx               | mac-arm64               |
| all_epitopes_by_cancer_{mhci or mhcii}.xlsx      | win64                   |
| all_variables_ranked_{mhci or mhcii}.xlsx        | linux64                 |
| ranked_epitopes_by_cancer_{mhci or mhcii}.xlsx   | linux64                 |
| top_epitopes_{mhci or mhcii}.xlsx                | linux64                 |
| top_epitopes_by_cancer_{mhci or mhcii}.xlsx      | linux64                 |
| filtered_epitopes_{cancer}_{mhci or mhcii}.txt   | linux64                 |
| optimized_epitopes_{mhci or mhcii}.xlsx          | linux64                 |
| optimized_epitopes_{cancer}_{mhci or mhcii}.txt  | linux64                 |
| population_coverage_results_{mhci or mhcii}.xlsx | linux64                 |


