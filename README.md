# HybridRocketEngine
This repository currently contains program files for initial research into hybrid engine models. The purpose of this research is to test and compare mathematical models.
Folder src contains Initial-Research and Model

# Setup:
## SSH GitHub Setup
If you already have an SSH key and have linked it with GitHub, you can skip this step. However, if it is an older RSA key then you will need to generate a new SSH key with the ED25519 encryption since MATLAB does not support RSA.
 - #### Generating Key: [Link](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)  
   MacOS Terminal: `ssh-keygen -t ed25519 -C "your_email@example.com"`  
   Windows PowerShell: `ssh-keygen -t ed25519 -C "your_email@example.com"`  
 - #### Viewing Key:  
   MacOS Terminal: `cat ~/.ssh/id_ed25519.pub`  
   Windows PowerShell: `Get-Content C:\Users\Your_Username\.ssh\id_ed25519.pub`  
 - #### Adding Key to GitHub: [Link](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)  
   Profile -> Settings -> SSH and GPG keys -> New SSH key  
 - #### Setting Up Git with MATLAB [Link](https://www.mathworks.com/help/matlab/matlab_prog/set-up-git-source-control.html)
## VS Code Setup
*Suggested Extensions: Python, MATLAB, Rainbow CSV*

#### VS Code Download: [link](https://code.visualstudio.com/download)  
#### Python3.11 Download [here](https://www.python.org/downloads/)

 - #### Clone Repository
   In GitHub Repository -> <> Code -> Copy SSH Link
   In VS Code -> Welcome Page -> Clone Git Repository... -> Paste SSH Link -> Select Directory

 - #### Setup Git Config for Local Machine
   In VS Code -> Terminal -> New Terminal  
   `git config --global user.name "Your Name"`   
   `git config --global user.email "you@example.com"`   

 - #### Create Virtual Environment  
   MacOS Terminal: `python3.11 -m venv env`  
   Windows PowerShell: `python -m venv env`  

 - #### Activate Virtual Environment  
   MacOS Terminal: `source env/bin/activate`  
   Windows PowerShell: `env\Scripts\activate.bat`  

 - #### Deactivate Virtual Environment  
   MacOS Terminal: `deactivate`  
   Windows PowerShell: `deactivate`  

 - #### Download the required packages
   MacOS Terminal: `pip install -r requirements.txt`  
   Windows PowerShell: `pip install -r requirements.txt`  
