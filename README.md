## The files relevant to the study are:
*fem_beam.m, which generates the Timoshenko beam FEM model,
*beam_mor_pmor, which is the actual application of the MOR and pMOR methods on the model
*soirka folder, which contains the implementation of the MOR subspace algorithm


# bm-mfem

bm-mfem is the Chair of Structural Mechanics' in-house FEM code, not currently necessary for this project. Please check the wiki for further information.

## First steps
To install bm-mfem on your computer you need some software:
*  MATLAB version > R2017a
*  Git (The versioninig software used by Gitlab). Get it from https://git-scm.com/ and do NOT change any options during installation.
*  GitHub Desktop (A GUI for Git). After installation, do not sign in to GitHub, as we are using Gitlab
*  GiD, if you want to use a pre- and post-processor

After all required software is installed, you can set up bm-mfem:

1.  Start a command shell (Windows+R -> cmd) and enter the following commands: 
    * `git config --global user.name "USERNAME"`, where USERNAME is the name chosen in Gitlab (the name in bold, when you click on your avatar in the upper right corner)
    * `git config --global user.email "EMAIL"`, where EMAIL is your TUM email adress shown in your Gitlab profile
2.  Start GitHub Desktop. Do not sign in, but clone a respsitory using the URL of this respository (https://gitlab.lrz.de/ga69tiz/bm-mfem.git)
3.  Start Matlab and move to the bm-mfem folder
4.  Execute `bm_fem.m`. This adds all required folders to your search path
5.  Type `runTests` to see if everything is working correctly

## Rules for developing
*  We use the feature branch workflow. Please see https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow.
*  If you have not yet been using git, create a test repository and familiarize yourself with it.
*  **Do not add binary files like `*.mat, *.fig`!**. Please only add plain text files (e.g. `*.m, *.txt`) to the repository.
*  This repository is not indended for private backups. Only commit things belonging to the `bm-mfem` project. To backup results or similar files, use your NAS drive or LRZ Sync + Share.
*  If you are some commits behind `master`, you can get them into your branch via `Branch -> Update from default branch` in GitHub Desktop
