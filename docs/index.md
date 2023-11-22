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
If a list of Python versions >=3.8.0 are shown, you may skip the rest of the Miniconda installation and 
proceed to the git installation process.

### Linux
To install Miniconda on Linus, run the following commands in Terminal::
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
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
**Note: These commands install the latest Miniconda compatible with arm64 or Apple Silicon M1/M2. If 
you have an older architecture,**