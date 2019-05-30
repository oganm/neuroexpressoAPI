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

pid = Sys.getpid()
cat(pid,file = '/var/run/neapi/neapi.pid',append = FALSE)



library(plumber)
plum = plumb('server.R')
plum$run(port=8020,swagger = TRUE, host="0.0.0.0")

