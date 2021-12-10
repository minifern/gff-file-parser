file_parse <- function(gffFile) {

  ##read in file

  Z01 <- readLines(gffFile)

  Z02 <- grep(pattern = "(^#)",

              x = Z01)

  Z03 <- Z01[-Z02]

  Z04 <- strsplit(Z03,

                  split = "\t",

                  fixed = TRUE)

  Z05 <- do.call(rbind, 

                 Z04)

  

  ##data to a data frame format

  gffData <- as.data.frame(Z05, 

                           stringsAsFactors = FALSE)

  

  ##coersion

  gffData[,4] <- as.integer(gffData[,4])

  gffData[,5] <- as.integer(gffData[,5])

  gffData[,7] <- ifelse(test = gffData[,7] == "+",

                        yes = 0L,

                        no = 1L)

  

  ##split by feature 

  geneData <- subset(gffData, gffData[,3] == 'gene')

  cdsData <- subset(gffData, gffData[,3] == 'CDS')

  

  ##isolating column data to use later

  cdsStart <- cdsData[,4]

  cdsStop <- cdsData[,5]

  geneStart <- geneData[,4]

  geneStop <- geneData[,5]

  cdsContig <- cdsData[,1]

  geneContig <- geneData[,1]

  cdsStrand <- cdsData[,7]

  geneStrand <- geneData[,7]

  cdsAttrib <- strsplit(cdsData[,9], split = ";", fixed = TRUE)

  

  ##isolate proteinID, product and geneID

  #product

  cdsProduct <- sapply(cdsAttrib, 

                       function(x)if(any(grepl("product=", x, fixed = TRUE))) {

                         str_extract(x[which(grepl("product=", 

                                                   x, 

                                                   fixed = TRUE))], 

                                     pattern = "(?<=product=)(.*)")

                       } else { 

                         NA_character_

                       },

                       USE.NAMES = FALSE,

                       simplify = TRUE)

  #proteinID

  cdsProteinID <- sapply(cdsAttrib,

                         function(x)if(any(grepl("protein_id=", x, fixed = TRUE))) {

                           str_extract(x[which(grepl("protein_id=", 

                                                     x, 

                                                     fixed = TRUE))], 

                                       pattern = "(?<=protein_id=)(.*)")

                         } else { 

                           NA_character_

                         },

                         USE.NAMES = FALSE, 

                         simplify = TRUE)

  #geneID

  geneID <- sapply(cdsAttrib, 

                   function(x)if(any(grepl("ID=", 

                                           x, 

                                           fixed = TRUE))) {

                     str_extract(x[which(grepl("ID=", 

                                               x, 

                                               fixed = TRUE))], 

                                 pattern = "(?<=ID=)(.*)")

                   } else { 

                     NA_character_

                   },

                   USE.NAMES = FALSE, 

                   simplify = TRUE)



  ##find matching CDS indices for each gene

  cdsMatchedInd <- vector(mode = "list", 

                          length = nrow(geneData))

  for(i in 1:nrow(geneData)) {

    cdsMatchedInd[[i]] <- which(cdsStart >= geneStart[i] &

                                  cdsStop <= geneStop[i] &

                                  cdsStart <= geneStop[i] &

                                  cdsStop >= geneStart[i] &

                                  cdsContig == geneContig[i] &

                                  cdsStrand == geneStrand[i])

    

    #ensure IDs are equal within the matches

    if(length(cdsMatchedInd[[i]]) != 0) {

      cdsMatchedInd[[i]] <- sapply(cdsMatchedInd[[i]],

                                   function(x) if(geneID[x] == geneID[cdsMatchedInd[[i]][1]]) {

                                     return(x)

                                   } else {

                                     return(NA_integer_)

                                   },

                                   USE.NAMES = FALSE,

                                   simplify = TRUE)

    } else {

      cdsMatchedInd[[i]] <- NA_integer_

    }

  }

  

  ##collapsed starts and stops

  collapsedLocs <- vector(mode = "list", 

                          length = length(cdsMatchedInd))

  ##matched cds info

  matchedCdsProduct <- vector(mode = "character", 

                              length = nrow(geneData))

  matchedCdsID <- vector(mode = "character", 

                         length = nrow(geneData))

  matchedCdsProteinID <- vector(mode = "character", 

                                length = nrow(geneData))



  for (i in seq_along(cdsMatchedInd)) {

    #collapsing starts and stops

    if(!is.na(cdsData[cdsMatchedInd[[i]], 4]) | !is.na(cdsData[cdsMatchedInd[[i]], 5])){

      collapsedLocs[[i]] <- paste(cdsData[cdsMatchedInd[[i]], 4], 

                                  "X", 

                                  cdsData[cdsMatchedInd[[i]], 5], 

                                  sep = "")

    } else {

      collapsedLocs[[i]] <- NA_character_

    }

    #matched info

    if(length(cdsMatchedInd[[i]] > 1)) {

      matchedCdsProduct[i] <- cdsProduct[cdsMatchedInd[[i]][1]]

      matchedCdsID[i] <- geneID[cdsMatchedInd[[i]][1]]

      matchedCdsProteinID[i] <- cdsProteinID[cdsMatchedInd[[i]][1]]

    } else {

      matchedCdsProduct[i] <- NA_character_

      matchedCdsID[i] <- NA_character_

      matchedCdsProteinID[i] <- NA_character_

    }

  }

  #collapse matches together

  collapsedLocs <- sapply(collapsedLocs, function(x) paste(x,

                                                           collapse = "Y",

                                                           sep = ""))

  

  ##create data frame with strand, starts, starts and stops, matches, geneId, cdsProduct, CDSproteinID

  parsedData <- data.frame("Strand" = geneStrand, 

                           "Start" = geneStart, 

                           "Stop" = geneStop, 

                           "Matches" = collapsedLocs, 

                           "Gene_ID" = matchedCdsID, 

                           "CDS_Product" = matchedCdsProduct, 

                           "CDS_Protein_ID" = matchedCdsProteinID)

  return(c(nrow(geneData), nrow(cdsData)))

  return(parsedData)

}
