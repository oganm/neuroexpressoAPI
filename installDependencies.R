dir.create('packages',showWarnings = FALSE)
libPath = normalizePath('packages')

.libPaths(libPath)

if(!'devtools' %in% installed.packages()[,1]){
    install.packages('devtools',lib = libPath,repos = 'https://cran.rstudio.com')
}

dependencies = readLines('depends.txt')

local = dependencies[!grepl('/',dependencies)]
remote = dependencies[grepl('/',dependencies)]

installed = installed.packages()

if(any(!local %in% installed)){
    install.packages(local[!local %in% installed],lib = libPath,repos = 'https://cran.rstudio.com')
}

for(package in remote){
    devtools::install_github(package,lib = libPath)
}