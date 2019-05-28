dir.create('packages',showWarnings = FALSE)
libPath = normalizePath('packages')

.libPaths(libPath)

if(!'devtools' %in% installed.packages()[,1]){
    install.packages('devtools',lib = libPath)
}

dependencies = readLines('depends.txt')

local = dependencies[!grepl('/',dependencies)]
remote = dependencies[grepl('/',dependencies)]

installed = installed.packages()

install.packages(dependencies[!dependencies %in% installed],lib = libPath)

for(package in remote){
    devtools::install_github(package,lib = libPath)
}