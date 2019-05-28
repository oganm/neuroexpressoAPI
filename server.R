library(ggplot2)
library(magrittr)
library(jsonlite)
dir.create('packages',showWarnings = FALSE)
libPath = normalizePath('packages')

.libPaths(libPath)

library(cowplot)
library(viridis)
library(ogbox)
library(dplyr)
library(reshape2)

exprs  = readRDS('Data/exprs.rds')
designs = readRDS('Data/designs.rds')
genes = readRDS('Data/genes.rds')

minValue = readRDS('Data/minValue.rds')
maxValue = readRDS('Data/maxValue.rds')
prop ='ShinyNames'
regionNames = 'Region'
ogbox::sourceGithub('oganm/brainGenesManuscript/R/cellColors.R')
coloring = cellColors()
coloring = c(coloring,
             ShreejoyGabaergic = 'pink',
             ShreejoyPurkinje = 'pink',
             "*Purkinje" = 'pink',
             ShreejoyPyramidal = 'pink',
             ShreejoyOligo = 'pink',
             ShreejoyAstrocyte= 'pink',
             ShreejoyThPosLC = 'pink',
             'Layer 4 Pyra' = 'blue',
             'Layer 2 3 Pyra' = 'blue',
             'Layer 6a Pyra' = 'blue',
             'Layer 6b Pyra' = 'blue',
             'Oligodendrocyte precursors' ='darkgreen',
             Endothelial = 'yellow')


createFrame = function(gene,
                       geneList,
                       expression,
                       design,
                       prop,
                       reference = 'Reference',
                       pmid = 'PMID',
                       coloring,
                       field = 'Gene.Symbol',
                       regionSelect,
                       color = T,
                       order = 'Cell type',
                       treeChoice,
                       treeSelected){
    
    treeSelectedNames  = sapply(treeSelected,function(x){x[1]})
    names(treeSelected) = treeSelectedNames
    mouseExpr = expression[,!is.na(regionSelect),]
    mouseDes = design[!is.na(regionSelect),]
    mouseGene = geneList
    
    
    
    if (len(treeSelected)==0){
        # treeSelected = design[hierarchyNames[[treeChoice]][len(hierarchyNames[[treeChoice]])]] %>% unique %>% trimNAs
        treeSelected = hierarchies[[treeChoice]] %>% unlist %>% names %>% gsub("^.*[.](?![ ,])",'',.,perl = T)
    }
    
    # if (order == 'A-Z'){
    #     treeSelectedNames = sort(treeSelectedNames)
    # }
    
    tree = hierarchyNames[[treeChoice]]
    
    # to create groups to display have the fields relevant to the selected tree and find indexes of the choices in it
    selectFrom = mouseDes %>% select_(.dots=tree)
    # groups = lapply(treeSelectedNames, function(x){
    #     selectFrom %>% apply(1,function(y){x %in% y}) %>% which
    # })
    groups = lapply(treeSelected, function(x){
        if (is.null(attr(x,'ancestry'))){
            selectFrom[,len(tree)] %in% x %>% which
        } else{
            out = selectFrom[,len(attr(x,'ancestry'))+1] %in% x[1] %>% which
            # limit selection to its ancestor in case there are leaves with the same name
            if (len(attr(x,'ancestry'))>0){
                limitToParent = lapply(1:len(attr(x,'ancestry')),function(i){
                    selectFrom[,tree[i]] %in% attr(x,'ancestry')[i] %>% which
                })
                out = intersectList(c(list(out),limitToParent))
            }
            return(out)
        }
        # selectFrom %>% apply(1,function(y){x %in% y}) %>% which
    })
    names(groups) = treeSelected # in case 
    
    
    while(groups %>% names %>% duplicated %>% any){
        names(groups)[groups %>% names %>% duplicated] %<>% paste0(' ')
    }
    
    if (order == 'A-Z'){
        groups = groups[groups %>% names %>% order]
    }
    expression = t(mouseExpr[mouseGene[,field] %in% gene,])
    frame = groups %>% melt
    colors = toColor(mouseDes$ShinyNames[frame$value],coloring)$col
    frame %<>% mutate(GSM = mouseDes$sampleName[value], 
                      gene = expression[value,], 
                      prop= L1,  
                      color = colors ,
                      reference = mouseDes$Reference[value],
                      rnaSeq = mouseDes$Platform[value] %in%  'RNAseq',
                      PMID = mouseDes$PMID[value]) %>% 
        select(GSM,gene,prop,color, reference,rnaSeq,PMID)
    
    # amygdala fix. if a region doesnt exist returns an empty matrix
    if (nrow(frame)==0){
        frame[,2] = integer(0)
    }
    
    frame$color = apply(col2rgb(frame$color),2,function(x){
        x = x/255
        rgb(x[1],x[2],x[3])
    })
    # if color is false, set all to black.
    if (!color){
        frame$color = "#000000"
    }
    
    frame$rnaSeq %<>%
        replaceElement(c("FALSE" = 'Microarray', 'TRUE' = 'RNAseq')) %$% 
        newVector %>% 
        factor()
    frame$`Data Source` = frame$rnaSeq
    
    
    
    # amygdala fix again
    if(nrow(frame)==0){
        frame = cbind(frame,data.frame(id=character(0)))
        return(frame)
    }
    return(frame)
}

load('memoReg.rda')

regionGroups = designs %>% lapply(function(x){
    out = memoReg(x,regionNames,prop,
                  regionHierarchy = regionHierarchy
    )
    
    names(out) = sapply(names(out),function(x){
        strsplit(x,split = '_')[[1]][1]
    })
    return(out)
})

hierarchyNames = list(NeuronTypes = c('MajorType','Neurotransmitter','ShinyNames'),
                      Methodology = c('Method', 'Reference'))

hierarchize = function(levels,design){
    out = vector(mode = 'list', length = len(unique(design[levels[1]]) %>% trimNAs))
    
    out = lapply(out,function(x){structure('',stselected = TRUE)})
    names(out) = unique(design[levels[1]]) %>% trimNAs %>% sort
    
    if ((len(levels)>1) & (nrow(design)>0)){
        out = lapply(names(out),function(x){
            hierarchize(levels[-1] ,design[design[,levels[1]] %in% x,])
        })
        names(out) = unique(design[levels[1]]) %>% trimNAs %>% sort
        for(i in 1:len(out)){
            if (len(out[[i]])==1 && names(out[[i]]) == names(out[i])){
                out[[i]] = structure('',stselected = TRUE)}
        }
    }
    return(out)
}

hierarchies = lapply(hierarchyNames, function(levels){
    hierarchize(levels,designs$GPL339[!is.na(designs$GPL339[,levels[len(levels)]]),])
})

ensemblIDs = gemmaAPI::readDataFile('Data/Generic_mouse_ensembl')


#' @get /nx_plot
#' @png (width = 600*1.5, height = 500*1.5)
nx_plot = function(gene = 'Ogn',ncbi = NULL, ensembl = NULL, region = 'Cortex',dataset = 'GPL339', fixed = FALSE){
    gc()
    
    if(is.null(ncbi) & is.null(ensembl)){
        geneRow = which(genes[[dataset]]$Gene.Symbol == gene)
    } else if(!is.null(ncbi)){
        geneRow = which(genes[[dataset]]$NCBIids == ncbi)
    } else if(!is.null(ensembl)){
        ncbi = ensemblIDs$NCBIids[ensemblIDs$ProbeName %in% ensembl]
        geneRow = which(genes[[dataset]]$NCBIids == ncbi)
    }
    
    if(is.null(ncbi) & is.null(ensembl)){
        field = 'Gene.Symbol'
    } else if(!is.null(ncbi)){
        gene = ncbi
        field = 'NCBIids'
    }
    
    expression = exprs[[dataset]][geneRow,]
    regionSelect = regionGroups[[dataset]][[region]]
    
    createFrame(gene = gene,
                geneList =  genes[[dataset]],
                expression = exprs[[dataset]],
                design = designs[[dataset]],
                prop = "ShinyNames",
                reference = "Reference",
                pmid = "PMID",
                coloring = coloring,
                field = field,
                regionSelect =  regionGroups[[dataset]][[region]],color = TRUE,order = 'Cell type',treeChoice = 'NeuronTypes',
                treeSelected = character(0)
    ) -> frame
    frame$prop %<>% factor(levels = unique(frame$prop))
    
    p = ggplot(frame, aes(y = gene, x = prop, fill = color,shape = rnaSeq)) + geom_point(size = 5) + scale_fill_identity() + 
        theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 16),
              axis.text.y = element_text(size = 16),
              axis.title.y = element_text(size= 16),
              legend.title = element_blank()) + xlab('') +
        scale_shape_manual(values=c(21,24)) + 
        ylab(paste(genes[[dataset]]$Gene.Symbol[geneRow],'log2 expression')) 
    
    if(fixed){
        p = p + coord_cartesian(ylim = c(minValue, maxValue))
    }
    
    print(p)
    
}


#' @get /nx_data
#' @serializer contentType list(type='text/json')
nx_data = function(gene = 'Ogn',ncbi = NULL, ensembl = NULL, region = 'Cortex',dataset = 'GPL339'){
    gc()
    
    
    if(is.null(ncbi) & is.null(ensembl)){
        geneRow = which(genes[[dataset]]$Gene.Symbol == gene)
    } else if(!is.null(ncbi)){
        geneRow = which(genes[[dataset]]$NCBIids == ncbi)
    } else if(!is.null(ensembl)){
        ncbi = ensemblIDs$NCBIids[ensemblIDs$ProbeName %in% ensembl]
        geneRow = which(genes[[dataset]]$NCBIids == ncbi)
    }
    
    if(is.null(ncbi) & is.null(ensembl)){
        field = 'Gene.Symbol'
    } else if(!is.null(ncbi)){
        gene = ncbi
        field = 'NCBIids'
    }
    
    expression = exprs[[dataset]][geneRow,]
    regionSelect = regionGroups[[dataset]][[region]]
    
    createFrame(gene = gene,
                geneList =  genes[[dataset]],
                expression = exprs[[dataset]],
                design = designs[[dataset]],
                prop = "ShinyNames",
                reference = "Reference",
                pmid = "PMID",
                coloring = coloring,
                field = field,
                regionSelect =  regionGroups[[dataset]][[region]],color = TRUE,order = 'Cell type',treeChoice = 'NeuronTypes',
                treeSelected = character(0)
    ) -> frame
    
    frame %<>% select(GSM, gene, prop, reference,PMID,`Data Source`)
    names(frame) = c('Sample','Expression','Cell Type','Reference','PMID','Data Source')
    frame %>% jsonlite::toJSON(pretty = TRUE)
}