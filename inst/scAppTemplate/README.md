# Michael Lattke Astrocye App

## Run App Locally

When prompted to update packages, try first the 'None' option. 

## Option A: Database Version
```
devtools::install_github(
    "FrancisCrickInstitute/hc",
    ref = "main", 
    auth_token = "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"
)
```
Once the installation is complete, run:

```
library(hc)

hc::run_app()
```



## Working on the code of this app on your local computer

### Get the present app from github

Go to github repo for this project https://github.com/FrancisCrickInstitute/hc and select the branch you'd like to work on. 

Get all files to your computer using the path given in the green code button like so:

```
git clone git@github.com:FrancisCrickInstitute/hc.git
```

Create a new branch for your improvement work

```
git checkout -b app_improvements
```

To render and run the app locally from the downloaded code, e.g. after you've made changes, do the following 
```
setwd("path/to/project/folder/hc")

devtools::check()

# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# Detach all loaded packages and clean your environment
golem::detach_all_attached()
# rm(list=ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload()

# Run the application
run_app()
```

To check for the name of your current github branch you can do
```
git branch
```

## Deploying this app on a shiny server
To deploy this app on a shiny server, download the content of the repo using the git clone approach described above, make sure cytoscape is installed and copy the PPI folder onto the shiny server. Once all relevant R-packages are installed, the app should run. 

Once you are happy with your improvements, you can push your development branch to the remote github repo by doing (inside the PPI folder)

```
git add -A

git commit -m "Description of the changes you've made"

git push origin [name_of_your_new_branch]
```
