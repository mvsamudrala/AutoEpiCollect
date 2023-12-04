# Welcome to AutoEpiCollect!
This is the AutoEpiCollect documentation site. Here, you will learn how to download, install, and use both the 
console and GUI versions of AutoEpiCollect in order to collect immunogenic epitopes for a vaccine 
targeting any cancer(s) of your choice. 

AutoEpiCollect is a software that uses web-scraping, tools from the Immune Epitope Database and Analysis Resource 
(IEDB), and machine learning to predict potentially immunogenic neoantigens for any oncogene. 

AutoEpiCollect GUI homepage:
![GUI for AutoEpiCollect](GUI-home.png)

AutoEpiCollect workflow chart:
![Logic workflow for AutoEpiCollect program](workflow-chart.png)

Please click the tabs above to go through the installation process and use documentation for AutoEpiCollect.

# Installation
## Miniconda
AutoEpiCollect is best run using Python 3.8.x. In order to have the required dependencies for AutoEpiCollect to run 
smoothly, we strongly recommend that you install Miniconda, the lighter version of Anaconda Python distribution. 
Miniconda is used to create multiple virtual Python environments, which are run on separate projects that 
require different dependencies. Following these steps will allow you to create a virtual environment with the 
exact dependencies needed. It is possible to run AutoEpiCollect by installing each dependency one by one, but we 
recommend following the steps outlined in this guide to set up AutoEpiCollect with minimal effort. 

First, check if a version of Miniconda or Anaconda is installed in your system path and is compatible with 
AutoEpiCollect. If you are on macOS or Linux, open Terminal. If you are on Windows, open Anaconda Prompt, if 
installed (if not installed, this means you must download Miniconda). Run the following command::

```bash
conda search python
```
If a list of Python versions >=3.8.0 appear, you may skip the rest of the Miniconda installation and 
proceed to the [git installation process](#git).


### Linux
To install Miniconda on Linux, run the following commands in Terminal::
```bash
cd ~
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
### Windows
To install Miniconda on Windows, run the following commands in Command Prompt::
```cmd
cd %HOMEPATH%
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe
start /wait "" miniconda.exe /S
del miniconda.exe
```
### macOS
To install Miniconda on macOS, run the following commands in Terminal::
```bash
cd ~
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
**Note: These commands install the latest Miniconda compatible with 64-bit Linux and Windows systems, as well 
as Apple Silicon M1/M2. If you have an older architecture or the commands above are not working properly, please go 
to the [Miniconda installation](https://docs.conda.io/projects/miniconda/en/latest/miniconda-other-installer-links.html)
website, download the Python 3.x installer compatible with your system, open the downloaded file, and follow the 
instructions to finish the Miniconda installation.**

## Git
Git is a version control system that tracks changes in computer files, mainly source code. AutoEpiCollect is hosted 
on GitHub, a platform that stores Git repositories. Installing Git will allow you to easily download the most 
up-to-date files needed to successfully run AutoEpiCollect from Terminal or Command Prompt. First, check if you 
already have Git installed on your system by running the following command in your Terminal or Command Prompt::
```bash
git --version
```
If the Git version number appears, then you can move on to the 
[AutoEpiCollect installation process](#install-autoepicollect). If not, then please navigate to the [Git downloads 
website](https://git-scm.com/downloads) and follow the instructions to complete the installation process.
## ChromeDriver

## Install AutoEpiCollect
**Note: For Windows users. Anaconda Prompt is a shell program similar to Command Prompt that gives users access to 
command line input. Anaconda Prompt comes with the Miniconda installation and is necessary for running Git commands 
and creating your virtual environments. For the next steps of installing AutoEpiCollect and running the program, 
please use the Anaconda Prompt shell found by searching "miniconda3" in your start menu and clicking the Anaconda 
Prompt program.**

To install AutoEpiCollect, you will need to copy all the necessary files from AutoEpiCollect's GitHub repository. 
Then, you must enter AutoEpiCollect's directory to access all the main files and programs. Enter the 
commands below in Terminal or Anaconda Prompt (on Windows) to gain access to AutoEpiCollect's directory::
```bash
git clone "https://github.com/mvsamudrala/AutoEpiCollect"
cd AutoEpiCollect
```
As stated in the Miniconda section of the installation process, the use of a conda virtual environment to install 
all the packages needed for AutoEpiCollect to run smoothly is highly recommended. Located in the AutoEpiCollect 
directory is a .yml file that contains the necessary dependencies for creating a conda virtual environment 
compatible with AutoEpiCollect. Follow the commands below to create and activate this virtual environment on 
your machine::
```bash
conda env create -n aec_venv
conda activate aec_venv
```
**Note: Creating the virtual environment might take some time depending on your machine. Please be prepared to wait 
up to 30 minutes after entering the first command for all the required dependencies to be installed.**

If you wish to switch back to your default environment, run this command to deactivate the existing virtual 
environment::
```bash
conda deactivate
```
## Testing AutoEpiCollect
To run the GUI program, enter the following command while the aec_venv environment is running::
```bash
python AutoEpiCollectGUI.py
```
