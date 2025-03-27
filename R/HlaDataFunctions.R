######## Mostly hla data processing
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' PyArd
#'
#' Load reticulate library in R.
#'
#' library(reticulate)
#' If using virtual environments, call use_virtualenv function.
#' use_virtualenv('~/py-ard/venv/')
#' Verify your Python environment
#' py_config()
#' This should show you the Python installation that reticulate is using.
#' Install py-ard in the Python environment
#'
#' py_install('py-ard')
#'
#'
PyArd<-function(...){ #3290 for latest imgt hla alleles
  pyd<-reticulate::import("pyard")
  ard = pyd$ARD(...)
  return(ard)
}



#' ConverDist2Sim
#'
#' Convert a distance matrix to a similarity matrix
#'
#'
ConverDist2Sim<-function(DistMat){
  n<-nrow(DistMat)
  CenteringMat<-matrix(-1/n, n,n)
  diag(CenteringMat)<-1-1/n
  OutMat<- -.5*CenteringMat%*%DistMat%*%CenteringMat
  rownames(OutMat)<-colnames(OutMat)<-rownames(DistMat)
  return(OutMat)
}


#' ConverSim2PCs
#'
#' Get the principal components from a similarity matrix
#'
#'
ConverSim2PCs<-function(SimMat){
  svdSimMat<-svd(SimMat)
  PCs<-SimMat%*%svdSimMat$v
  return(PCs)
}





######Note: The desmat rows could be adjusted using allele frequencies,

#' GetAlleleFeatures
#'
#' Creates a list that includes an allelesindata list for alleles as they match
#' to allelenames and a design matrix for these matched alleles
#'
#'
#'
GetAlleleFeatures<-function(alleles,field="A", allelenames){
  ard<-PyArd()
  FieldforARD<-paste(c("HLA-",field,"*"),collapse="")
  FieldforARD2<-paste(c(field,"\\*"),collapse="")
  FieldforARD3<-paste(c("HLA-",field,"\\*"),collapse="")
  ardres<-sapply(paste(FieldforARD,alleles, sep=""),function(x){
    x<-c(unlist(strsplit(x,split="/",fixed=TRUE)))[1]
    ard$redux(x,"W")
  })
  allelesindata<-lapply(ardres,function(x){
    allelespyard<-strsplit(gsub(FieldforARD2,"",gsub(FieldforARD3,"",x)), "/",fixed=TRUE)[[1]]
    allelespyard[allelespyard%in%allelenames]
  })
  namestoberesolved<-unique(names(allelesindata)[which(sapply(allelesindata, length)==0)])
  print(namestoberesolved)
  namesresolved<-lapply(namestoberesolved, function(x){

    y<-gsub("HLA-","",x)
    outvec<-c()
    xexp<-tryCatch(ard$expand_mac(y), error=function(e){
      print(y)
      print("failed")
      return(y)
      })
    if (length(xexp)>0){
    for (rd in xexp){
    out<-ard$redux(rd, "W")
    outvec<-c(outvec, strsplit(gsub(FieldforARD2,"",gsub(FieldforARD3,"",out)), "/",fixed=TRUE)[[1]])
    }
    } else {
      print("error2")
      outvec<-x
      }
    allelespyard<-lapply(ardres,function(x){
      allelespyard<-strsplit(gsub(FieldforARD2,"",gsub(FieldforARD3,"",x)), "/",fixed=TRUE)[[1]]
    })
    out<-unique(c(unlist(allelespyard[allelespyard%in%outvec])))
    if (length(out)>0){return(out)}else {return(x)}
    })
  names(namesresolved)<-namestoberesolved
  print(namesresolved)

  for (namei in namestoberesolved){
    for (nameii in which(names(allelesindata)%in%namei)){
    allelesindata[[nameii]]<-namesresolved[names(namesresolved)%in%namei][[1]]
  }
  }
  DesMat<-sapply(allelesindata, function(x){
    y<-rep(0,length(allelenames))
    y[allelenames%in%x]<-1/max(c(sum(allelenames%in%x),1))
    y
  })
  DesMat<-t(DesMat)
  colnames(DesMat)<-allelenames
  return(list(allelesindata=allelesindata,DesMat=DesMat))
}


#' Calc6MatchingDists
#'
#' Gets 6 sequence based distance matrices for a given set of donor
#' recipient alleles
#'
#'
#'
Calc6MatchingDists<-function(d_typ1,d_typ2,r_typ1,r_typ2,field="A", allelenames,DistAlleles){
  df1<-GetAlleleFeatures(alleles=d_typ1,field=field, allelenames=allelenames)
  df2<-GetAlleleFeatures(alleles=d_typ2,field=field, allelenames=allelenames)
  rf1<-GetAlleleFeatures(alleles=r_typ1,field=field, allelenames=allelenames)
  rf2<-GetAlleleFeatures(alleles=r_typ2,field=field, allelenames=allelenames)
  r12<-sapply(1:length(df1$allelesindata), function(ri){
    print(ri)
    print(rf1$allelesindata[[ri]])
    print(rf2$allelesindata[[ri]])
    aaf<-function(ri){mean(DistAlleles[rf1$allelesindata[[ri]],rf2$allelesindata[[ri]]])}
    tryCatch(aaf(ri), error=function(e){NA})
  })
  d12<-sapply(1:length(df1$allelesindata), function(ri){

    aaf<-function(ri){mean(DistAlleles[df1$allelesindata[[ri]],df2$allelesindata[[ri]]])}
    tryCatch(aaf(ri), error=function(e){NA})
  })

  dr11<-sapply(1:length(df1$allelesindata), function(ri){
    aaf<-function(ri){mean(DistAlleles[df1$allelesindata[[ri]],rf1$allelesindata[[ri]]])}
    tryCatch(aaf(ri), error=function(e){NA})

  })
  dr21<-sapply(1:length(df1$allelesindata), function(ri){
    aaf<-function(ri){mean(DistAlleles[df1$allelesindata[[ri]],rf2$allelesindata[[ri]]])}
    tryCatch(aaf(ri), error=function(e){NA})

  })
  dr12<-sapply(1:length(df1$allelesindata), function(ri){
    aaf<-function(ri){mean(DistAlleles[df2$allelesindata[[ri]],rf1$allelesindata[[ri]]])}
    tryCatch(aaf(ri), error=function(e){NA})


  })
  dr22<-sapply(1:length(df1$allelesindata), function(ri){
    aaf<-function(ri){mean(DistAlleles[df2$allelesindata[[ri]],rf2$allelesindata[[ri]]])}
    tryCatch(aaf(ri), error=function(e){NA})


  })
  outdf<-data.frame(d12=d12,r12=r12,dr11=dr11,dr21=dr21,dr12=dr12,dr22=dr22)
  return(outdf)
}




#' CalcHLAGenotypeFeatures
#'
#' Given sequence based features  calculate features for a genotype
#' at a HLA field
#'
#'
CalcHLAGenotypeFeatures<-function(typ1,typ2, field="A", allelenames,FeatAlleles){
  f1<-GetAlleleFeatures(alleles=typ1,field=field, allelenames=allelenames)
  f2<-GetAlleleFeatures(alleles=typ2,field=field, allelenames=allelenames)

  F1<-f1$DesMat%*%FeatAlleles
  F2<-f2$DesMat%*%FeatAlleles

  outMat<-cbind(F1,F2)
  return(outMat)
}

