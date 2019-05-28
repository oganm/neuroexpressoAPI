set_file_wd = function(){
    command = commandArgs(trailingOnly = FALSE)
    
    file = gsub('--file=','',command[grepl('--file',command)])
    if(length(file) == 1){
        setwd(dirname(file))
    }
}
set_file_wd()

libPath = normalizePath('packages')
.libPaths(libPath)

setwd(here::here())



pid = readLines('lastpid')

system(paste0('kill -9 ',pid))