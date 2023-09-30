# HybridRocketEngine
This repository currently contains program files for initial research into hybrid engine models. The purpose of this research is to test and compare mathematical models.
Folder src contains Initial-Research and Model

## Setup:
### SSH GitHub Setup
If you already have an SSH key and have linked it with GitHub, you can skip this step. However, if it is an older RSA key then you will need to generate a new SSH key with the ED25519 encryption since MATLAB does not support RSA.
 - Generating Key: [Link](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
 - Adding Key to GitHub: [Link](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
 - Setting Up Git with MATLAB [Link](https://www.mathworks.com/help/matlab/matlab_prog/set-up-git-source-control.html)
### VS Code Setup
*Suggested Extensions: Python, MATLAB, Rainbow CSV*

## Clone this repository and set up a virtual environment for Python version 3.11. If you don't have it then you can install [here](https://www.python.org/downloads/)

```python3.11 -m venv env``` (MacOS) or ```python -m venv env``` (Windows)

## Open a Python Virtual Environment

MacOS
```source env/bin/activate```

Windows
```env\Scripts\activate.bat```

To deactivate a virtual Environment
```deactivate```

## Download the required packages

```pip install -r requirements.txt```