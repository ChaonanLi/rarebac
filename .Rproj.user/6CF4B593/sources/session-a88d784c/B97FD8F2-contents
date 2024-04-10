# Create a directory
createDir <- function(rootPath, dir){
    path <- paste(rootPath, dir, sep = '/')
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    return(path)
}

# Save PDF
savePDF <- function(object, path, filename, width = 8, height = 6){
    pdf(paste(path, filename, sep = '/'), width = width, height = height)
    print(object)
    dev.off()
}

# Save xlsx with single sheet
saveSingleXlsx <- function(object, path, filename, row.names = T, col.names = T){
    if (nrow(object) == 0){
        object1 <- as.data.frame(matrix(nrow = 2, ncol = ncol(object)))
        colnames(object1) <- colnames(object)
        object <- object1
    }
    xlsx::write.xlsx(x = object, file = paste(path, filename, sep = '/'), 
                     col.names = col.names, row.names = row.names, append = F)
} 

# Save xlsx with mutiple sheets
saveMutipleXlsx <- function(objectList, path, filename, row.names = T, col.names = T){
    for (i in seq(length(objectList))){
        sheetname <- names(objectList[i])
        dataframe <- objectList[[i]]
        if (nrow(dataframe) == 0){
            dataframe1 <- as.data.frame(matrix(nrow = 2, ncol = ncol(dataframe)))
            colnames(dataframe1) <- colnames(dataframe)
            dataframe <- dataframe1
        }
        append <- TRUE
        if (i == 1) append <- FALSE
        xlsx::write.xlsx(x = dataframe, file = paste(path, filename, sep = '/'), 
                         col.names = col.names, row.names = row.names, 
                         append = append, sheetName = sheetname)
    }
}

# Save stdout
saveStdout <- function(StdoutVar, Outdir, prefix){
    file <- file(paste(Outdir, paste(prefix, '.txt', sep = ''), sep = '/'))
    sink(file)
    print(StdoutVar)
    sink()
    close(file)
}
