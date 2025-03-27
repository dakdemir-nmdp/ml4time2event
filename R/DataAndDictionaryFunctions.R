## Module: Data and Dictionary Functions
## This file contains functions for processing dictionaries and data.
## Sections include:
##   - Dictionary processing functions (e.g. readDict, filllasttoDict, bindDicts, CatVarsDict, etc.)
##   - Data processing functions (e.g. readData, ApplyDictToData, RemoveMonoVarsData, etc.)
########1- Mostly dictionary processing functions
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'





#' readDict
#'
#' Import a dictionary that is in standard format
#' The dictionary should be in a csv, xlsx, or sas7bdat file
#' @param file path to dictionary file
#' @param filetype type of dictionary file (csv, xlsx, sas7bdat)
#' @param ... additional arguments to pass to read.csv, openxlsx::read.xlsx,
#' haven::read_sas
#' @return data frame with variable information
#' @examples
#' \dontrun{
#' readDict(file = "data/dictionary.csv", filetype = "csv")
#' }
#' @export
readDict <- function(file = NULL, filetype="csv", ...) {
## This code reads in a dictionary file and returns a data frame
## with the variable information
## Input: file = path to dictionary file
##        filetype = type of dictionary file (csv, xlsx, sas7bdat)
## Output: DICDATA = data frame with variable information
  if (!is.null(file)) {
    if (filetype=="csv"){
      DICDATA <- as.data.frame(read.csv(file, ...))
    } else if (filetype=="xlsx") {
      DICDATA <- as.data.frame(openxlsx::read.xlsx(file, ...))
    } else if (filetype=="sas7bdat"){
     DICDATA <- as.data.frame(haven::read_sas(file, ...))
    } else {
      stop("Provide filetype (one of csv, xlsx, sas7bdat)!")
    }
    if (sum(
      colnames(DICDATA) %in% c(
        "Variable",
        "Description",
        "Type",
        "Value",
        "Label",
        "Notes"
      ))==6
    ) {
      DICDATA <- DICDATA[, colnames(DICDATA) %in% c("Variable",
                                                    "Description",
                                                    "Type",
                                                    "Value",
                                                    "Label",
                                                    "Notes")]

      return(DICDATA)
    } else {
      stop(
        "Colnames of data should contain \n
           Variable, Description, Type,Value, Label, and Notes. \n
           Check help file for dictionary format"
      )
    }
  } else {
    stop("Provide file location")
  }
}



#' filllasttoDict
#'
#' Utility function to put dictionary in a special format, fill forward
#'
#' @param dict data frame with dictionary information
#' @return data frame with dictionary information
#' @export

filllasttoDict <- function(dict) {
# This code uses the last non-NA value in a data frame column to impute  NA values
# in that column.
filllastvector <- function(x) {
    y <- x
    filli <- NA
    for (i in 1:length(y)) {
      if (!is.na(y[i])) {
        filli <- y[i]
      } else {
        y[i] <- filli
      }
    }
    return(c(unlist(y)))
  }
  DICT <- apply(dict, 2, filllastvector)
  DICT <- as.data.frame(DICT)
  colnames(DICT) <- colnames(dict)
  return(DICT)
}



#' bindDict
#'
#' Utility function to combine two dictionaries in a serial fashion
#'
#' @param dict1 data frame with dictionary information
#' @param dict2 data frame with dictionary information
#' @return data frame with dictionary information
#' @export
bindDicts <- function(Dict1, Dict2) {
  # This function binds two dictionaries together
  # Dict1 is a dictionary of variable names and their corresponding
  #   descriptions
  # Dict2 is a dictionary of variable names and their corresponding
  #   descriptions
  # The function returns a dictionary that contains all of the variable
  #   names and descriptions from Dict1 and Dict2. If a variable name
  #   is in both Dict1 and Dict2, then the description from Dict1 is
  #   used.
  Names1 <- unique(Dict1$Variable)
  Names2 <- unique(Dict2$Variable)
  Dict2int <- Dict2[Dict2$Variable %in% intersect(Names2, Names1), ]

  Dict2rem <- Dict2[Dict2$Variable %in% setdiff(Names2, Names1), ]
  DICT <- dplyr::bind_rows(Dict1, Dict2rem)
  rownames(DICT) <- 1:nrow(DICT)
  return(DICT)
}

bindDicts <- function(Dict1, Dict2) {
  Names1 <- unique(Dict1$Variable)
  Names2 <- unique(Dict2$Variable)
  Dict2int <- Dict2[Dict2$Variable %in% intersect(Names2, Names1), ]

  Dict2rem <- Dict2[Dict2$Variable %in% setdiff(Names2, Names1), ]
  DICT <- dplyr::bind_rows(Dict1, Dict2rem)
  rownames(DICT) <- 1:nrow(DICT)
  return(DICT)
}



#' CatVarsDict
#'
#' Utility function to get dictionary for the categorical variables
#'
#' @param dict data frame with dictionary information
#' @return data frame with dictionary information
#' @export
CatVarsDict <- function(dict) {
  return(dict[dict$Type == "Categorical", ])
}


#' ConVarsDict
#'
#' Utility function to get dictionary for the continuous variables
#' @param dict data frame with dictionary information
#' @return data frame with dictionary information
#'
#' @export
ConVarsDict <- function(dict) {
  return(dict[dict$Type == "Continuous", ])
}


#' CharVarsDict
#'
#' Utility function to get only the Character variables
#' @param dict data frame with dictionary information
#' @return data frame with dictionary information
#'
#' @export
CharVarsDict <- function(dict) {
  return(dict[dict$Type == "Character", ])
}


#' OneValCat2NumDict
#'
#' Utility function to correct the class of Continuous variables if they were
#' misidentified
#'
#' @param dict data frame with dictionary information
#' @param vars character vector of variable names
#' @return data frame with dictionary information
#' @export
OneValCat2NumDict <- function(dict, vars) {
  tabx <- table(dict$Variable)
  numCats <- names(tabx)[which(tabx == 1)]
  dict[dict$Variable %in% intersect(numCats, vars),]$Type <-
    "Continuous"
  return(dict)
}

#' DetectIdVarsDict
#'
#' Utility function to get only the id variables
#' @param dict data frame with dictionary information
#' @return character vector of variable names
#' 
#' @export
DetectIdVarsDict <- function(dict, idvars = NULL) {
  variables <- unique(dict$Variable)
  strout <-
    ((((
      stringr::str_detect(variables, "id") |
        stringr::str_detect(variables, "ccn")
    ) |
      stringr::str_detect(variables, "pseudo")
    ) |
      stringr::str_detect(variables, "team")
    ) | stringr::str_detect(variables, "case"))
  return(variables[unique(c(strout, idvars))])
}

#' NumIdVarstoCharDict
#'
#' Utility function to convert id vars to Character
#'
#' @param dict data frame with dictionary information
#' @param idvars character vector of variable names
#' @return data frame with dictionary information
#' @export
NumIdVarstoCharDict <-
  function(dict,
           idvars = c("`#`", "pseudoid", "pseudoccn")) {
    for (vari in idvars) {
      dict[dict$Variable == vari, "Type"] <- "Character"
    }
    return(dict)
}



#' NumIdVarstoCharDict
#'
#' Utility function to detect NA strings in a dictionary
#' @param dict data frame with dictionary information
#' @param knownNaLabels character vector of known NA labels
#' @param knownNaValues character vector of known NA values
#' @return data frame with dictionary information
#'
#' @export
DetectNAStringsDict <-
  function(dict,
           knownNaLabels = NULL,
           knownNaValues = NULL) {
    nalabels <-
      c(
        "missing",
        "unknown",
        "not answ",
        "not answered",
        "n/a",
        "unspecified",
        "nan")
    navalues <-
      c(".b", ".c", ".f", ".y", "<NA>", "99", "999")
    if (!is.null(knownNaLabels)) {
      nalabels <- union(nalabels, knownNaLabels)
    }
    if (!is.null(knownNaValues)) {
      navalues <- union(navalues, knownNaValues)
    }

    detectNAs <- function(dict) {
      whichNAlabels <-
        sapply(nalabels, function(x) {
          which(stringr::str_detect(tolower(dict$Label), x))
        })
      whichNAvalues <-
        sapply(navalues, function(x) {
          which(stringr::str_detect(tolower(dict$Value), x))
        })
      NACats <- c(whichNAlabels, whichNAvalues)
      return(NACats)
    }
    Nalocs <- detectNAs(dict)
    dict$NACats <- "VALUE"
    for (listi in Nalocs) {
      dict$NACats[listi] <- "NAVAL"
    }
    return(dict)
  }





#' DateVarstoDateDict
#'
#' Utility function to convert date vars to Date
#' @param dict data frame with dictionary information
#' @param datevars character vector of variable names
#' @return data frame with dictionary information
#'
#' @export
DateVarstoDateDict <- function(dict, datevars = NULL) {
  if (!is.null(datevars)) {
    for (vari in datevars) {
      dict[dict$Variable == vari, 'Type'] <- "Date"
    }
  }
  return(dict)
}

#' finddatevars
#'
#' Utility function to find possible date variables
#' @param names character vector of variable names
#' @param knowndatevars character vector of known date variables
#' @return character vector of variable names
#'
#' @export
finddatevars<-function(names, knowndatevars=NULL){
  datevars<-union(names[grepl("dte", names)],names[grepl("date", names)])
  datevars<-union(datevars, knowndatevars)
  datevars
}


#' addDiseaseSpecColsDict
#'
#' process the dictionary to add disease specific / not applicable features
#'
#' @param dict data frame with dictionary information
#' @return data frame with dictionary information
#'  @export
addDiseaseSpecColsDict<-function(dict){

  dict[dict$Value%in%".E", ]$NACats<-"NotApp"

  DiseaseSpecificVars<-unique(dict[dict$Label%in%"N/A, other disease",]$Variable)
  alldiseases<-dict[dict$Variable%in%"disease",]$Label

  DiseaseSpecificVarsDescr<-unique(dict[dict$Label%in%"N/A, other disease",]$Description)


  dict$DisSpec<-NA
  for (i in 1:length(alldiseases)){
    before<-dict$DisSpec[dict$Description%in%DiseaseSpecificVarsDescr[grepl(alldiseases[i], DiseaseSpecificVarsDescr)]]
    toadd<-rep(alldiseases[i], length(dict$DisSpec[dict$Description%in%DiseaseSpecificVarsDescr[grepl(alldiseases[i], DiseaseSpecificVarsDescr)]]))
    after<-apply(cbind(before, toadd), 1, function(x){paste(na.omit(x), collapse=" ")})
    dict$DisSpec[dict$Description%in%DiseaseSpecificVarsDescr[grepl(alldiseases[i], DiseaseSpecificVarsDescr)]]<-after
  }

  dict
}




########2- Mostly data processing functions


#' readData
#'
#' Import a data that is in standard format
#' @param file character vector of file location
#' @param filetype character vector of file type (csv, xlsx, sas7bdat)
#' @return data frame with data
#'
#'  @export
readData <- function(file = NULL, filetype="csv",...) {
  if (!is.null(file)) {
    if (filetype=="csv"){
      DATA <- as.data.frame(read.csv(file, ...))
    } else if (filetype=="xlsx") {
      DATA <- as.data.frame(openxlsx::read.xlsx(file, ...))
    } else if (filetype=="sas7bdat"){
      DATA <- as.data.frame(haven::read_sas(file, ...))
    } else {
    stop("Provide filetype (one of csv, xlsx, sas7bdat)!")
    }
    colnames(DATA) <- tolower(colnames(DATA))
    return(DATA)
  } else {
    stop("Provide data location!")
  }
}



#' ApplyDictToData
#'
#' apply dictionary to data
#' @param dict data frame with dictionary information
#' @param data data frame with data
#' @return data frame with data
#'
#' @export
ApplyDictToData <- function(dict, data) {
  data <- as.data.frame(data)
  for (vari in unique(dict$Variable)) {
    dictx <- dict[dict$Variable == vari, ]
    VarName <- vari
    if (VarName %in% colnames(data)) {
      x <- c(unlist(data[, colnames(data) %in% VarName]))
      x[x %in% dictx$Value[dictx$NACats == "NAVAL"]] <- NA
      data[, colnames(data) %in% VarName] <- x
      comment(data[, colnames(data) %in% VarName]) <-
        dictx$Description[[1]]

      if (dictx$Type[[1]] %in% c("Categorical")) {
        Values <- trimws(dictx$Value[dictx$NACats != "NAVAL"])
        Labels <- trimws(dictx$Label[dictx$NACats != "NAVAL"])
        if (length(Values) > 1) {
            x <- data[, colnames(data) %in% VarName]
            data[, colnames(data) %in% VarName] <- x
            data[, colnames(data) %in% VarName] <-
              tryCatch(
                factor(trimws(data[, colnames(data) %in% VarName]), levels = Values, labels = Labels), # nolint
                error = function(e) {
                  print(VarName)
                  print("Warning: ")
                  print(e)
                  as.character(data[, colnames(data) %in% VarName])
                },
                warning=function(w){
                  print(VarName)
                  print(w)
                  print("Warning")
                }
              )
        } else {
          data[, colnames(data) %in% VarName] <-
            as.character(data[, colnames(data) %in% VarName])
        }


      } else if (dictx$Type[[1]] %in% c("Character")) {
        data[, colnames(data) %in% VarName] <-
          as.character(data[, colnames(data) %in% VarName])
      } else {
        data[, colnames(data) %in% VarName] <-
          as.numeric(data[, colnames(data) %in% VarName])
        x <- data[, colnames(data) %in% VarName]

        x[x == 99] <- NA
        data[, colnames(data) %in% VarName] <- x
      }
      comment(data[, colnames(data) %in% VarName]) <-
        dictx$Description[[1]]
    }
  }
  return(data)
}


#' RemoveMonoVarsData
#'
#' If a variable has one value remove it
#' @param data data frame with data
#' @return data frame with data
#'
#' @export
RemoveMonoVarsData <- function(data) {
  numcats <- apply(data, 2, function(x)
    length(table(x)))
  data[, numcats > 1]
}

#' RemoveMonoVarsData
#'
#' If a variable has only NA's remove it
#'
#'
RemoveAllNAVars <- function(data) {
  toremove <- apply(data, 2, function(x)
    length(x) == sum(is.na(x)))
  data[, !toremove]
}



#' findDiseaseSpecVarsData
#'
#' Utility function to try to detect disease spec variables from data
#' @param data data frame with data
#' @param pin numeric value between 0 and 1 for the proportion of non NA values in the disease group
#' @param pout numeric value between 0 and 1 for the proportion of NA values in the other groups
#' @return list of variables that are disease specific
#'
#'  @export
findDiseaseSpecVarsData <- function(data, pin = .5, pout = .9) {
  varstocheck <- setdiff(colnames(data), "disease")
  dislist <- vector(mode = "list")
  for (disi in na.omit(unique(data$disease))) {
    varlist <- as.character(c())
    for (vari in varstocheck) {
      print(vari)
      x <- data[, colnames(data) %in% vari]
      y <- which(data$disease == disi)
      yother <- which(data$disease != disi)
      isnax <- which(is.na(x))
      notisnax <- which(!is.na(x))
      if (((sum(notisnax %in% y) / length(notisnax)) > pin) &
          ((sum(isnax %in% yother) / length(isnax)) > pout)) {
        varlist <- c(varlist, vari)
      }
    }

    if (length(varlist) > 0) {
      dislist <- c(dislist, list(c(unlist(varlist))))
      names(dislist)[length(dislist)] <- disi
    }
  }

  dislist

}




#' processDiseaseSpecVarsData
#'
#' Utility function to process disease spec variables
#'
#' @param data data frame with data
#' @param dislist list of variables that are disease specific
#' @return data frame with data
#' @export
processDiseaseSpecVarsData <- function(data, dislist) {
  processonedisease <- function(disease) {
    vars <- dislist[[disease]]
    whichd <- which(data$disease == disease)
    whichnotd <- which(data$disease != disease)

    for (vari in vars) {
      data[, vari] <- as.character(data[, vari])
      x <- data[whichnotd, vari]
      isnax <- is.na(x)
      data[intersect(whichnotd, which(isnax)), vari] <- "NotApp"
      data[, vari] <- as.factor(data[, vari])
    }
    return(data)

  }

  for (disease in names(dislist)) {
    data <- processonedisease(disease)
  }
  return(data)
}




#' RemoveMissinVarsData
#'
#' Utility function to remove variables with more than a certain percentage
#' of missing
#' @param data data frame with data
#' @param maxprop numeric value between 0 and 1 for the maximum proportion of NA values
#' @return data frame with data
#'
#' @export
RemoveMissinVarsData <- function(data, maxprop = .2) {
  props <- apply(data, 2, function(x)
    sum(is.na(x)) / length(x))
  data[, props <= maxprop]
}

#' RemoveMissinRecordsData
#'
#' Utility function to remove records with more than a certain percentage
#' of missing
#' @param data data frame with data
#' @param maxprop numeric value between 0 and 1 for the maximum proportion of NA values
#' @return data frame with data
#'
#'  @export
RemoveMissinRecordsData <- function(data, maxprop = .2) {
  props <- apply(data, 1, function(x)
    sum(is.na(x)) / length(x))
  data[props <= maxprop, ]
}



#' getcharcolsData
#'
#' Utility function to det the names of Character columns in data
#' @param data data frame with data
#' @return character vector with names of character columns
#'
#'  @export
getcharcolsData <- function(data) {
  charcols <- sapply(colnames(data), function(x) {
    is.character(data[, colnames(data) %in% x])
  })
  return(names(which(charcols)))
}







#' ImputeMissinRecordsData
#'
#' Imputation of data using missRanger
#'
#' @param data data set to be imputed
#' @param dontuse the variables not to use when imputing,
#' thesev ariables will not be used for imputation or for imputing other
#'  variables.
#' @inheritParams missRanger::missRanger
#' @return Imputed data
#'
#' @seealso [missRanger::missRanger()] which this function wraps.
#' @examples
#' \dontrun{
#' DATAImputed<-DATA%>%ImputeMissinRecordsData()
#' }
#'
#' @export
ImputeMissinRecordsData <- function(data, dontuse = NULL, ...) {
  namesindata<-colnames(data)
  charcols <- getcharcolsData(data)
  charcols <- union(charcols, dontuse)
  datatoimp<-data[, !(colnames(data) %in% charcols)]
  datatoimpnames<-colnames(datatoimp)
  colnames(datatoimp)<-make.names(colnames(datatoimp))
  datatoimp <-
    missRanger::missRanger(datatoimp, ...)
  colnames(datatoimp)<-datatoimpnames
  data[,!(colnames(data) %in% charcols)]<-datatoimp
  return(data)
}




#' RemoveRareCategoriesData
#'
#' Put NA's to rare (freq less than a specified percentage) categories of
#' factor (Categorical) variables
#' @param data data frame with data
#' @param minfreq numeric value between 0 and 1 for the minimum frequency of a category
#' @return data frame with data
#'
#'
#'
#' @export

RemoveRareCategoriesData <- function(data, minfreq = .01) {
  for (vari in colnames(data)) {
    if (is.factor(data[, colnames(data) %in% vari])) {
      tabx <- table(data[, colnames(data) %in% vari])
      ftabx <- tabx / sum(tabx)
      removelevels <- names(tabx)[ftabx <= minfreq]
      x <- data[, colnames(data) %in% vari]
      x[x %in% removelevels] <- NA
      data[, colnames(data) %in% vari] <- x
      data[, colnames(data) %in% vari] <-
        droplevels(data[, colnames(data) %in% vari])
    }
  }

  return(data)
}


#' RemoveRareBinaryVarsData
#'
#' Put NA's to rare (freq less than a specified percentage) categories of
#' binary variables
#' @param data data frame with data
#' @param minfreq numeric value between 0 and 1 for the minimum frequency of a category
#' @return data frame with data
#'
#'  @export
RemoveRareBinaryVarsData <- function(data, minfreq = .01) {
  for (vari in colnames(data)) {
    if (length(unique(data[, colnames(data) %in% vari]))==2) {
      tabx <- table(data[, colnames(data) %in% vari])
      ftabx <- tabx / sum(tabx)
      removelevels <- names(tabx)[ftabx <= minfreq]
      x <- data[, colnames(data) %in% vari]
      x[x %in% removelevels] <- NA
      data[, colnames(data) %in% vari] <- x
    }
  }

  return(data)
}






#' CollapseRareCategoriestoOtherData
#'
#' Put 'Other' to rare (freq less than a specified percentage) categories of
#' factor (Categorical) variables
#' @param data data frame with data
#' @param minfreq numeric value between 0 and 1 for the minimum frequency of a category
#' @return data frame with data
#'
#'  @export
CollapseRareCategoriestoOtherData <- function(data, minfreq = .01) {
  for (vari in colnames(data)) {
    if ((is.factor(data[, colnames(data) %in% vari]) | is.character(data[, colnames(data) %in% vari]))) {
      x<-as.character(data[, colnames(data) %in% vari])
      tabx<-table(x)
      if (length(tabx)>3){
      namestocollapse<-names(tabx[tabx<ceiling(minfreq*nrow(data))])
      x[x%in%namestocollapse]<-"Other"
      x<-as.factor(x)
      data[, colnames(data) %in% vari]<-x
      }
    }
  }
  return(data)
}



#' droplevelsoffactorsData
#'
#' Drop unused factor levels from data
#' @param data data frame with data
#' @return data frame with data
#'
#'
#'  @export
droplevelsoffactorsData <- function(data) {
  for (vari in colnames(data)) {
    if (is.factor(data[, vari])) {
      data[, vari] <- droplevels(data[, vari])
    }
  }
  data
}

#' findvarsnamesthatrepeatData
#'
#' Find variable names that are similar to other variable names
#' @param data data frame with data
#' @return a vector with variable names that are similar to other variable names
#'
#'
#'  @export
findvarsnamesthatrepeatData <- function(data) {
  outmat <- c()
  for (vari in colnames(data)) {
    if (length(grep(vari, setdiff(colnames(data), vari))) > 0) {
      outmat <-
        c(outmat, paste(setdiff(colnames(data), vari)[grep(vari, setdiff(colnames(data), vari))], collapse =
                          ", "))
      names(outmat)[length(outmat)] <- vari
    }
  }

  outmat
}




#' ZeroOneScalerData
#'
#' Scale numeric variables to the range of 0 and 1
#' @param data data frame with data
#' @return a list with the data frame with data and the min and max values of each variable
#'
#'
#'  @export
ZeroOneScalerData<-function(data){
  datanew<-data
  minxvec<-c()
  maxxvec<-c()
  for (i in 1:ncol(datanew)){
    if (is.numeric(datanew[,i])){
      minx<-min(datanew[,i])
      maxx<-max(datanew[,i])
      datanew[,i]<-(datanew[,i]-minx)/(max(maxx-minx,1))
      minxvec<-c(minxvec,minx)
      maxxvec<-c(maxxvec,maxx)
    } else {
      minxvec<-c(minxvec,NA)
      maxxvec<-c(maxxvec,NA)

    }
  }
  names(minxvec)<-names(maxxvec)<-colnames(datanew)
  return(list(data=datanew, minxvec=minxvec,maxxvec=maxxvec ))
}



#' ZeroOneScalerApplierData
#'
#' Scale numeric variables to the range of 0 and 1 using given mins and maxs
#' @param data data frame with data
#' @param mins a vector with the minimum values of each variable
#' @param maxs a vector with the maximum values of each variable
#' @return a data frame with the data
#'
#'
#'  @export
ZeroOneScalerApplierData<-function(data, mins, maxs){
  datanew<-data
  for (i in colnames(datanew)){
    if (is.numeric(datanew[,i])){
      minx<-mins[i]
      maxx<-maxs[i]
      datanew[,i]<-(datanew[,i]-minx)/(max(maxx-minx,1))
    }
  }
  return(datanew)
}

#' UndoZeroOneScalerApplierData
#'
#' Un-scale numeric variables using given mins and maxs
#'
#'
#'  @export
UndoZeroOneScalerApplierData<-function(data, mins, maxs){
  datanew<-data

  for (i in colnames(datanew)){
    if (is.numeric(datanew[,i])){
      minx<-mins[i]
      maxx<-maxs[i]
      datanew[,i]<-max(maxx-minx,1)*datanew[,i]+minx
    }
  }
  return(datanew)
}





#' NumVarstCatsData
#'
#' Make categorcal variables from numeric variables
#' @param data data frame with data
#' @param numgroups a vector with the number of groups for each variable
#' @param cuts a vector with the cut points for each variable
#' @return a data frame with the data
#'
#'
#'  @export
NumVarstCatsData<-function(data,numgroups=NULL, cuts=NULL){
  quantcat<-function (x, q = 4, na.rm = TRUE, ...)
  {
    if (length(q) == 1) {
      q <- seq(0, 1, length.out = q + 1)
    }
    quant <- round(quantile(x, q, na.rm = na.rm))
    dups <- duplicated(quant)
    if (any(dups)) {
      flag <- x %in% unique(quant[dups])
      retval <- ifelse(flag, paste("[", as.character(x), "]",
                                   sep = ""), NA)
      uniqs <- unique(quant)
      reposition <- function(cut) {
        flag <- x >= cut
        if (sum(flag, na.rm = na.rm) == 0) {
          return(cut)
        }
        else {
          return(min(x[flag], na.rm = na.rm))
        }
      }
      newquant <- sapply(uniqs, reposition)
      retval[!flag] <- as.character(cut(x[!flag], breaks = newquant,
                                        include.lowest = TRUE, ...))
      levs <- unique(retval[order(x)])
      retval <- factor(retval, levels = levs)
      mkpairs <- function(x) {
        sapply(x, function(y) if (length(y) == 2)
          y[c(2, 2)]
          else y[2:3])
      }
      pairs <- mkpairs(strsplit(levs, "[^0-9+\\.\\-]+"))
      rownames(pairs) <- c("lower.bound", "upper.bound")
      colnames(pairs) <- levs
      closed.lower <- rep(F, ncol(pairs))
      closed.upper <- rep(T, ncol(pairs))
      closed.lower[1] <- TRUE
      for (i in 2:ncol(pairs)) {
        if (pairs[1, i] == pairs[1, i - 1] && pairs[1, i] ==
            pairs[2, i - 1]) {
          closed.lower[i] <- FALSE
        }
      }
      for (i in 1:(ncol(pairs) - 1)) {
        if (pairs[2, i] == pairs[1, i + 1] && pairs[2, i] ==
            pairs[2, i + 1]) {
          closed.upper[i] <- FALSE
        }
      }
      levs <- ifelse(pairs[1, ] == pairs[2, ], pairs[1, ],
                     paste(ifelse(closed.lower, "[", "("), pairs[1, ],
                           ",", pairs[2, ], ifelse(closed.upper, "]", ")"),
                           sep = ""))
      levels(retval) <- levs
    }
    else {
      retval <- cut(x, quant, include.lowest = TRUE, ...)
    }
    return(retval)
  }

  cutcat<-function (x, cuts, na.rm = TRUE, ...)
  {
    quant <- c(min(x),cuts, max(x))
    retval <- cut(x, quant, include.lowest = TRUE, ...)

    return(retval)
  }


  if (is.null(cuts)){
  datanew<-data

  for (i in 1:ncol(datanew)){
    if (is.numeric(datanew[,i])){
      if (length(unique(datanew[,i]))>(numgroups+3)){
      datanew[,i]<-quantcat(datanew[,i],q = numgroups-1)
      }
    }
  }
  } else {

    datanew<-data

    for (i in 1:ncol(datanew)){
      if (is.numeric(datanew[,i])){
        if (length(unique(datanew[,i]))>(length(cuts)+3)){
          datanew[,i]<-cutcat(datanew[,i],cuts = cuts)
        }
      }
    }
  }



  return(datanew)
}





########3- Mostly data summary functions


#' pairwiserelationshipsDataSummmary
#'
#' calculate pairwise correlations among variables using the BCD correlation
#' @param data data frame with data
#' @return a matrix with variable names and correlations
#'
#'  @export
pairwiserelationshipsDataSummmary <- function(data) {
  Cormat <- matrix(0, ncol(data), ncol(data))
  for (i in 1:(ncol(data) - 1)) {
    xi <- model.matrix( ~ data[, i] - 1)
    for (j in (i + 1):ncol(data)) {
      xj <- model.matrix( ~ data[, j] - 1)
      dcorij <- energy::bcdcor(xi, xj)
      Cormat[i, j] <- dcorij
      pbapply::pbupdate(i/ncol(data))

    }
  }
  Cormat <- Cormat + t(Cormat) + diag(ncol(data))
  rownames(Cormat) <- colnames(Cormat) <- colnames(data)
  return(Cormat)
}




#' gethighcorvarsDataSummmary
#'
#' Filter correlations based on a threshold and list names
#' @param pmat correlation matrix
#' @param corcutoff correlation cutoff
#' @return a matrix with variable names and correlations
#'
#'
#'  @export
gethighcorvarsDataSummmary <- function(pmat, corcutoff = .8) {
  corvarsmat <- NULL
  ExpVars<-rownames(pmat)
  for (vari in 1:(length(ExpVars) - 1)) {
    for (varj in (vari + 1):length(ExpVars)) {
      if (pmat[vari, varj] > corcutoff) {
        corvarsmat <-
          rbind(corvarsmat, c(ExpVars[vari], ExpVars[varj], pmat[vari, varj]))
      }
    }
  }
  return(corvarsmat)
}


#' OneAgainstRestCorDataSummmary
#'
#' Find variables correlated most with the rest
#' @param data data frame with data
#' @return a vector with correlations of each variable with the rest
#'
#'  @export
OneAgainstRestCorDataSummmary <- function(data) {
  Corvec <- rep(0,ncol(data))
  for (i in 1:(ncol(data))) {
    xi <- model.matrix( ~ data[, i] - 1)
      xrest <- model.matrix( ~ . - 1, data=data[,-i])
      dcorij <- energy::bcdcor(xi, xrest)
      Corvec[i] <- dcorij
      pbapply::pbupdate(i/ncol(data))
  }
  names(Corvec) <-  colnames(data)
  return(Corvec)
}



#' SummaryTableDataSummmary
#'
#' Make a data summary of a list of variables in the data
#' @param data data frame with data
#' @param UseVars a vector of variables to include in the summary
#' @return a gtsummary table
#'
#'
#'  @export
SummaryTableDataSummmary <-
  function(data,
           UseVars) {
    table2 <-data[, UseVars]%>%
    gtsummary::tbl_summary() %>%
      gtsummary::modify_header(label = "**Variable**") %>% # update the column header
      gtsummary::bold_labels()%>%
      gtsummary::italicize_labels()
    table2
  }



#' MakelabelledSASDict
#'
#' Make a dictionary from labelled data
#' @param fileloc location of the sas file
#' @return a labelled dictionary
#'
#'
#'  @export
MakelabelledSASDict<-function(fileloc){
  dataraw<-haven::read_sas(fileloc)
  dictionary <- labelled::generate_dictionary(dataraw)
  dictionary
}


#' A prepared dictionary for HSCT data from CIBMTR
#' 
#'
#' @docType data
#'
#'
"DICT"





#' ReplaceOutlierNumValsData
#'
#' Replace the outlier values for numeric vars
#' 
#' @param data data frame with data
#' @param multIQR multiplier for IQR
#' @param minnumgroup minimum number of groups to consider a variable
#' @return a data frame with outliers replaced
#'
#'
#'
#'  @export
ReplaceOutlierNumValsData<-function(data,multIQR=1.5,minnumgroup=10){
  datanew<-data
  minxvec<-c()
  maxxvec<-c()
  for (i in 1:ncol(datanew)){
    if (is.numeric(datanew[,i])){
      if (length(unique(datanew[,i]))>minnumgroup){
      print(colnames(datanew)[i])
      x<-datanew[,i]
      Q75<-quantile(x, .75, na.rm=TRUE)
      Q25<-quantile(x, .25, na.rm=TRUE)
      IQR=Q75-Q25
     UL<- Q75+IQR*multIQR
     LL<- Q25-IQR*multIQR
     x[x>UL]<-UL
     x[x<LL]<-LL
     datanew[,i]<-x
      }
    }
    }
datanew
}


####################################


#' MakeTestDataConfWithTrainData
#'
#' Make sure test data does not have levels that are not observed in the training
#' @param traindata data frame with training data
#' @param testdata data frame with test data
#' @return a data frame with test data
#'
#'
#'  @export
MakeTestDataConfWithTrainData<-function(traindata,testdata){
  testdata<-testdata[,colnames(testdata)%in%colnames(traindata)]
  for (vari in colnames(testdata)){
    if (is.factor(testdata[, vari])){
      valsintrain<-unique(traindata[, vari])
      x<-testdata[, vari]
      x[!(x%in%valsintrain)]<-NA
      x<-factor(as.character(x), levels=levels(traindata[, vari]))
      testdata[, vari]<-x
    }
  }
  testdata
}



#' RemoveEmptySpacesData
#'
#' Remove empty spaces from factors
#' @param DATA data frame with data
#' @return a data frame with data
#'
#'
#' @export
RemoveEmptySpacesData<-function(DATA){
for (vari in colnames(DATA)){
  if (is.factor(DATA[,vari])){
    x<-trimws(as.character(x), "both")
    DATA[,vari]<-as.factor(DATA[,vari])
  }
}
  DATA
}
