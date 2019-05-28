if(!'devtools' %in% installed.packages()[,1]){
    install.packages('devtools')
}

dependencies = readLines('depends.txt')

local = dependencies[!grepl('/',dependencies)]
remote = dependencies[grepl('/',dependencies)]

installed = installed.packages()

install.packages(dependencies[!dependencies %in% installed])

for(package in remote){
    devtools::install_github(package)
}