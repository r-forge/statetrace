################################################################################
#################  SETUP sta OBJECT ############################################
################################################################################
# sta object (made by calls to staMake) is a list of lists
# Top level: all= results summarizing over participants
#            ss = list with one entry per participant
# ss is a list with one entry for each participant in the experiment with
#   attributes "slevs" and "dlevs" and "nt" giving names of state and dimension
#   factor levels and number of trace levels, common to all participants
#
#   Each entry of ss is named `1`,`2`,...  and has an attribute "datasource"
#   giving the data file name from which the data was obtained (or the
#   participant column level corresponding to that participant when data all
#   obtained from one combined column.
#
# Each entry in ss is a list with components (only the first is filled in on
# initial creation, including prior probability calculation):
#   d=augmented data representation (see getdesign in ssMake)
#   m=statistcs for each model: t=trace, m=monotonic
#   p=plots for each model: e=encompassing, t=trace, m=monotonic
# Each element in m has: control=ancilliary info used in getting stats
################################################################################



################################################################################
#######  READ IN .txt DATA FILES AND MAKE UP sta OBJECT USING staMake ##########
################################################################################

# Data are read in from rectangular (i.e., same number of columsn in each row)
# tab delimited text files. Rows with NA in the final column are stripped.
#
# Two types are supported:
#   "sum", containing summed response frequencies with one line for each design
#          cell
#   "bin" containing binary reponse indicators with one line for each response.
#
# Two designs (with correspondingly different data files) are supproted, a no
# baseline ("B0" or "2AFC") or state baseline ("B2" or "Yes/No") design.
#
# There may be one file for data from all participants (see example files:
# sum2AFC.txt, bin2AFC.txt, sumYN.txt, binYN.txt) or one file for each
# partcipant (e.g., seperate files for each set of lines corresponding to one
# participant in the example files) or set of participants.
#
# For all file types/designs a single input file containing one or more
#   participants data sets is structured as follows:
# column 1=participant identifier,
# column 2=state factor     (must be 2 levels, can be character or numeric)
# column 3=dimension factor (must be 2 levels, can be character or numeric),
# column 4=trace factor     (must be integer, 1:n levels)
# NOTE: Trace level codes must be bigger (in the sense of a sort order, i.e.,
#       can be alpha as well as numeric) for condtions that should be associated
#       with greater accuracy, but otherwise may be different for different
#       participants/states/dimesnions (but must have the same number (T) of
#       unique values. These will then be converted to 1:T
#       e.g., in the 2AFC example data, trace=1/2/3=study for 66/200/800ms for
#             dimension="U" but 200/800/1800ms for dimesnion="I"
#       e.g., in the YN example data, trace=1/2/3=study for 33/100/267ms for
#             dimension="U" but 267/800/2048ms for dimesnion="I"
#
# Each of the four file types has a different structure for the remaining
# columns, with "sum" files having 6 columns and "bin" files 5 columns.
#
# For bin2AFC files the last column (column 5) has a 1 for a CORRECT choice, and
#     a 0 for an INCORRECT choice.
# For sum2AFC files the second last column (column 5) contains the number of
#     CORRECT responses in a cell and the last column (column 6) has the TOTAL
#     number of responses in a cell.
#
# For binYN files the last column (column 5) has a 1 for a YES choice, and
#     a 0 for a NO choice.
# For sum2AFC files the second last column (column 5) contains the number of
#     YES responses in a cell and the last column (column 6) has the TOTAL
#     number of responses in a cell.
#
# NOTE where multiple input files are specified to staMake, which indicates that
#      each corresponds to one particpant, there is no participant identifier
#      column so subtract one from all of the above column number specifications
#
# NOTE files for YN (B0) designs have extra lines corresponding to responses to
#      "noise" or "new" stimuli. There can only be TWO types of such stimuli,
#      one for each state level. THESE LINES ARE INDICATED BY "NA" ENTRIES
#      IN THE DIMENSION AND TRACE FACTOR COLUMNS.
#
# NOTE no provision is made for a one-baseline (B1) design as that can be
#      analyzed by treating it as a B0 (no-baseline) design, e.g., analyze only
#      the "Yes" data.
#

staMake <- function(staname="",fnams=NULL,folder="",extension="txt",
                    usecols=NULL,header=TRUE,sep="\t",na.strings="NA",
                    multiparticipant=FALSE,
                    a=1){

  # staname: string name of object to make in global environment
  # If it names an existing sta object data will be appended to that, if
  #   not a completely new sta object is created

  # fnams: A character vector with either a single element specifying
  #   the name of a file in the working directory or the path + name of a file
  #   in another directory, or with multiple elements, containing such
  #   specifications.
  # folder, extension: fills fnams with "folder//*.extension"
  # NOTE: If fnams is null and folder="" (default) on Windows fnams is
  #       filled with a call to "choose.files". On other systems, or if
  #       folders!="", fnams is filled with "folder//*.extension". Use
  #       folders="." to force this method in the current directory
  # NOTE: Where more than one element in fnams, each is assumed to specify data
  #         for one participant with no column specifying participants
  #         (participant numbering added automatically).
  # a: first and second beta shape parameters for prior calculation
  # usecols: vector of integers or names specifying which coloumns to use
  #          from each data file
  # header, sep, na.strings = file format arguements with read.delim defaults
  # NOTE that all lines with NA in the final column are stripped
  #
  # output is an sta object

  ssMake=function(fdats,snums,snams,a=1) {
  # fdats:  list of arrays of frequency data for each participant
  #         Array Format: tdsn=tx2x2x2 array
  #                       t=trace, d=dimension, s=state, n=n/N
  #                       attr(tdsn,"b"): NULL (B0 design) or 2x2 = state x n/N
  # snumnam: vector of sequinteger participant names, P1, P2 ...
  # snams: data source (or levels from single file) used to set attr datasource
  # a: first and second beta shape parameters for prior calculation

    getdesign <- function(fdat,snum,a=1) {
    # Sets up dat object for one element of fdats
    # INCLUDING CALCULATION OF PRIOR PROBABILITIES
    # Output list has entries:
    # 1) m,n,s1,s2: sucesses, num. trial, and two shape parameter matrices for
    #      each chain (one row per chain)
    # 2) odr: list of lists specifying order constraints to apply to each chain,
    # 3) priorp: prior probability of each chain in set
    #            NB: product is prior prob. for set
    # 4) chain.names: design cell(s) for each chain, D=dimension, S=state
    # 5) dsn: string with design name (B2, B0)

      if (!all(dim(fdat)[-1]==c(2,2,2)))
        stop("Data must be trace x state(2) x dimension(2) x m/n(2) array")
      if (any(is.na(fdat))) stop("NAs in data not supported")
      mold=fdat[,,,1]; nold=fdat[,,,2]
      nt=dim(mold)[1] # number of trace levels
      if (all(is.null(attr(fdat,"b")))) {           # B0 design
        s1=t(matrix(as.vector(fdat[,,,1]),ncol=4))
        s2=t(matrix(as.vector(fdat[,,,2]),ncol=4))
        odr=list(c1=list(1:nt),c2=list(1:nt),c3=list(1:nt),c4=list(1:nt))
        pp=rep(1/factorial(nt),4); dsn="B0"
        chain.names=c("D1S1","D2S1","D1S2","D2S2")
      } else {                                      # B2 design
        b=attr(fdat,"b")
        if (!all(dim(b)==c(2,2)))
          stop("Baseline must be state(2) x m/n(2) matrix")
        s1=cbind(as.vector(b[,1]),t(matrix(as.vector(fdat[,,,1]),ncol=2)))
        s2=cbind(as.vector(b[,2]),t(matrix(as.vector(fdat[,,,2]),ncol=2)))
        odr=list(c(1,2:(nt+1)),c(1,(nt+2):(2*nt+1))); odr=list(c1=odr,c2=odr)
        pp=rep(((factorial(nt))^-2)/(2*nt+1),2); dsn="B2"
        chain.names=c("S1","S2")
      }
      m=s1; n=s2
      s2=s2-s1+a; s1=s1+a  # first and second beta shape parameters
      list(sname=snum,m=m,n=n,s1=s1,s2=s2,odr=odr,priorp=pp,dsn=dsn,
           chain.names=chain.names)
    }

    nss=length(fdats); out=vector(mode="list",length=nss)
    j=0
    for (i in 1:length(fdats)) {
      dat=getdesign(fdats[[i]],snum=snums[i],a=a)
      out[[i]]=list(
        d=dat,                                      # dat in dsn form
        m=list(t=list(control=NULL),     # trace model
               m=list(control=NULL)      # monotonic model
        ),
        # plots for e=encompasing, trace and monotonic
        p=list(p=list(e=NULL,t=NULL,m=NULL),z=list(e=NULL,t=NULL,m=NULL))
      )
      attr(out[[i]],"datasource")=snams[i]
    }
    out
  }

  makeTraceInterger=function(sdt) {
    tmp=tapply(sdt[,3],list(sdt[,1],sdt[,2]),function(x){sort(unique(x))})
    nt=length(tmp[1,1][[1]])
    for (i in dimnames(tmp)[[1]])
      for (j in dimnames(tmp)[[2]]) {
        if (length(tmp[i,j][[1]])!=nt)
          stop("Number of trace factor levels is not the same in all cells")
        for (k in 1:nt)
          sdt[sdt[,1]==i & sdt[,2]==j & sdt[,3]==tmp[i,j][[1]][k],3]=k
    }
    as.numeric(sdt[,3])
  }

  getfdat=function(rdat,slevs,dlevs,nt) {
    bdat=rdat[is.na(rdat[,2]),-c(2,3)]
    rdat=rdat[!is.na(rdat[,2]),]
    rdat[,3]=makeTraceInterger(rdat[,1:3])
    if ( max(rdat[,3])!=nt )
      stop("Number of trace levels does not match previous files")
    sl=levels(rdat[,1])
    if ( is.null(sl) ) sl=sort(unique(rdat[,1]))
    dl=levels(rdat[,2])
    if ( is.null(dl) ) dl=sort(unique(rdat[,2]))
    if ( !all(sl==slevs) )
      stop("State levels do not match levels in previous files")
    if ( !all(dl==dlevs) )
      stop("Dimension levels do not match levels in previous files")
    if ( dim(rdat)[2]==5 ) { # summary
      n=tapply(rdat[,4],list(rdat[,3],rdat[,2],rdat[,1]),sum)
      N=tapply(rdat[,5],list(rdat[,3],rdat[,2],rdat[,1]),sum)
      if ( dim(bdat)[2]!=0 ) {
        bn=tapply(bdat[,2],bdat[,1],sum)
        bN=tapply(bdat[,3],bdat[,1],sum)
      }
    } else {                 # binary
      if ( !all(sort(unique(rdat[,4]))) %in% c(0,1) )
        stop("Data column must contain only 0 or 1 values")
      n=tapply(rdat[,4],list(rdat[,3],rdat[,2],rdat[,1]),sum)
      N=tapply(rdat[,4],list(rdat[,3],rdat[,2],rdat[,1]),length)
      if ( dim(bdat)[2]!=0 ) {
        if ( !all(sort(unique(bdat[,2]))) %in% c(0,1) )
        stop("Data column must contain only 0 or 1 values")
        bn=tapply(bdat[,2],bdat[,1],sum)
        bN=tapply(bdat[,2],bdat[,1],length)
      }
    }
    dn=list(paste("t",1:dim(n)[1],sep=""),c("d1","d2"),c("s1","s2"),c("n","N"))
    out=array(c(n,N),dim=c(dim(n),2),dimnames=dn)
    if ( dim(bdat)[1]!=0 ) {
      dn=list(c("s1","s2"),c("n","N"))
      attr(out,"b")=array(c(bn,bN),dim=c(2,2),dimnames=dn)
    }
    out
  }

  if (!is.character(staname) || nchar(staname)==0)
    stop("Name of sta object must be be given as quoted text")
  if (exists(staname,envir=.GlobalEnv)) {
    sta=get(staname,envir=.GlobalEnv)
    cat("\nADDING DATA TO AN EXISTING sta OBJECT!\n")
    if (class(sta)!="sta")
      stop(paste("Existing sta object",staname,"does not have class sta"))
  }

  if(sep=="tab")   {sep="\t"}
  if(sep=="space")   {sep=""}
  if(sep=="comma")     {sep=","}

  # fill fnams
  if (is.null(fnams))
    if ( Sys.info()["sysname"]=="Windows" && folder=="" ) {
      fnams=choose.files(
        caption="Select .txt data file(s) for state-trace analysis",
        filters = Filters["txt",]
      )
      fnams=c(fnams[-1],fnams[1]) # corrects for weird choose.files ordering
    } else fnams=dir(folder,full.names=TRUE,
                     pattern=glob2rx(paste("*",extension,sep=".")))

  # read in the .txt files, check and format data and place in a list
  fdats=vector(mode="list",length=0)
  snum=0; snams=character()
  for ( i in fnams ) {
    cat(paste("Reading in data file",i,"\n"))
    rawdats=read.delim(i,header=header,sep=sep,na.strings=na.strings)
    if ( !is.null(usecols) ) {
      if ( is.numeric(usecols) && max(usecols) > dim(rawdats)[2])
        stop("Column numbers specified in \"usecols\" too large")
      if ( is.character(usecols) && !all(usecols %in% names(rawdats)) )
        stop("Column names not in data file")
      rawdats=rawdats[,usecols]
    }
    ncols=dim(rawdats)[2]
    rawdats=rawdats[!is.na(rawdats[,ncols]),]
    if (multiparticipant) {
      ncols=ncols-1
      rawdat=rawdats[,-1]
    } else rawdat=rawdats
    bindat=switch(ncols-3,TRUE,FALSE,NA)
    if ( is.na(bindat) )
      stop(paste("Data file does not have right number of columns\n",
          "You may also have tried to read in a multiparticipant file",
          "without specifying multiparticipant=T."))
    # State factor
    if ( any(is.na(rawdat[,1])) )
      stop("NAs not allowed in state factor")
    slevs=levels(rawdat[,1])
    if ( is.null(slevs) ) slevs=sort(unique(rawdat[,1]))
    if (i==fnams[1]) slevs1=slevs
    if ( length(slevs)!=2 )
      stop(paste("The data must have exactly 2 state factor levels.\n",
          "You may also have tried to read in a multiparticipant file",
          "without specifying multiparticipant=T."))
    if ( !all(slevs==slevs1) )
      stop("Participants have different state levels.\n")
    # Dimension Factor
    nad=is.na(rawdat[,2])
    if ( any(nad) ) isB2=TRUE
    dlevs=levels(rawdat[,2])
    if ( is.null(dlevs) ) dlevs=sort(unique(rawdat[,2]))
    if (i==fnams[1]) dlevs1=dlevs
    if ( length(dlevs)!=2 )
      stop(paste("The data must have exactly 2 dimension factor levels.\n",
        "You may also have tried to read in a multiparticipant file",
        "without specifying multiparticipant=T."))
    if ( !all(dlevs==dlevs1) )
      stop("Participants have different dimension levels.\n")
    # Trace Factor
    nat=is.na(rawdat[,3])
    if ( !all(nat==nad) )
     stop(paste("Dimension and Trace Factors must have NAs in the same rows.",
        "You may also have tried to read in a multiparticipant file",
        "without specifying multiparticipant=T."))
    nt=max(makeTraceInterger(rawdat[!is.na(rawdat[,2]),1:3]))
    if (i==fnams[1]) nt1=nt
    if ( nt1!=nt )
      stop("Participants have different numbers of trace levels.\n")
    # If adding to previous sta object check if consistent
    if ( (exists("sta",inherits=FALSE)) ) {
      if ( !all(attr(sta$ss,"slevs")==slevs) )
        stop("State levels do not match levels in old sta object")
      if ( !all(attr(sta$ss,"dlevs")==dlevs) )
        stop("Dimension levels do not match levels in old sta object")
      if ( attr(sta$ss,"nt")!=nt )
        stop("Number of trace levels does not match old sta object")
    }
    if ( multiparticipant ) {
      snamsi=unique(as.character(rawdats[,1]))
      snams=c(snams,snamsi)
      for ( j in snamsi ) {
        snum=snum+1
        rawdat=rawdats[rawdats[,1]==j,-1]
        fdats[[as.character(snum)]]=getfdat(rawdat,slevs,dlevs,nt)
      }
    } else {
      snams=c(snams,i)
      snum=snum+1
      fdats[[as.character(snum)]]=getfdat(rawdat,slevs,dlevs,nt)
    }
  }
  snum = 1:snum
  if ( exists("sta",inherits=FALSE) ) snum = snum+length(sta$ss)
  ss=ssMake(fdats,snum,snams,a=a)
  dsns=unlist(lapply(ss,function(x){x$d$dsn}))
  if ( !all(dsns==dsns[1]) )
    stop("Participants have different baseline designs.")
  names(ss)=snum
  dsn=ss[[1]]$d$dsn
  if ( !exists("sta",inherits=FALSE) ) {
    attr(ss,"slevs")=slevs; attr(ss,"dlevs")=dlevs; attr(ss,"nt")=nt
    names(ss)=snum
    sta=list(all=list(p=list(e=NULL,t=NULL,m=NULL)),ss=ss)
    class(sta) <- "sta"
    cat(paste("\nObject",staname,
              "has state levels",paste(slevs,collapse=","),
              "dimension levels",paste(dlevs,collapse=","),
              "and",nt,"trace levels, with baseline design",dsn,"\n\n"))
  } else {
    slevsi=attr(sta$ss,"slevs"); dlevsi=attr(sta$ss,"dlevs")
    nti=attr(sta$ss,"nt"); dsni=sta$ss[[1]]$d$dsn
    cat(paste("New data has state levels",paste(slevs,collapse=","),
              "dimension levels",paste(dlevsi,collapse=","),
              "and",nt,"trace levels, with baseline design",dsn,"\n"))
    if (!all(slevs==slevsi) || !all(dlevs==dlevsi) || nt!=nti || dsn!=dsni)
      stop("New data doesn't have design matching object existing object")
    for (i in snum) sta$ss[[as.character(i)]]=ss[[as.character(i)]]
  }
  dss=unlist(lapply(sta$ss,function(x){attr(x,"datasource")}))
  dups=duplicated(dss)
  if ( any(dups) ) {
    cat("\nData sources are duplicated for the following participants:\n")
    print(dss[dups]); cat("\n")
  }
  assign(staname,sta,envir=.GlobalEnv)
  paste("Object",staname,"has been created/updated")
}


# stanames is single character name, in which case stored samples are managed
# or a vector of character names in which case the corresponding objects
# are joined and saved to an object with name given by staname (by default the
# first name in stanames).
# With defult settings nothing is done to stored samples but their size is
# reported. nkeep parameters allow their size to be reduced, and keepbestm=T
# first removes all but the best monotonic samples (as recoreded in
# $m$control$nm rather than in the actual stored sample)
staManage <- function(stanames="",staname=stanames[1],
                      checkConverge=FALSE,nmcmc=5000,mcmcSave="",
                      nkeepe=Inf,nkeept=Inf,
                      nkeepm=Inf,keepbestm=FALSE) {

  getmcmcl=function(ssi,nmcmc=5000) {
    tFC=ssi$p$t$s
    nt= dim(tFC)[1]; ntnd=nt*dim(tFC)[2]; nsim=dim(tFC)[3]
    tFC=aperm(tFC,c(3,2,1))
    nlist=nsim%/%nmcmc
    if (nlist>0) {
      out=vector(mode="list",length=nlist)
      tFC=array(tFC[1:(nmcmc*nlist),,],dim=c(nmcmc,nlist,ntnd))
      dimnames(tFC)[[3]]=
        as.vector(t(outer(ssi$d$chain.names,1:nt,paste,sep="T")))
      for (i in 1:nlist) out[[i]]=mcmc(tFC[,i,])
      out=mcmc.list(out)
      if (nlist<=1) attr(out,"mR")=NA else
        attr(out,"mR")=gelman.diag(out,autoburnin=FALSE,transform=TRUE)[[2]]
      attr(out,"badnmcmc")=nsim%%nmcmc
    } else {out=NA; attr(out,"badnmcmc")=nsim%%nmcmc; attr(out,"mR")=NA}
    out
  }

  if ( !is.character(stanames) || !is.vector(stanames) )
    stop("Name(s) of sta object(s) must be be given as a character vector")
  if ( !is.character(staname) || !is.vector(staname) || length(staname)!=1)
    stop("Save name must be be given as a character vector of length one")
  if ( any(sapply(c(stanames,staname),nchar)==0) )
    stop("Names of sta objects cannot be empty strings")
  ok=sapply(stanames,exists,envir=.GlobalEnv)
  if (any(!ok))
    stop(paste("Object(s)",paste(stanames[!ok],collapse=","),
               "not in the global environment"))
  stas=mget(stanames,envir=.GlobalEnv)
  ok=unlist(lapply(stas,function(x){class(x)=="sta"}))
  if (any(!ok))
    stop(paste("Object(s)",paste(stanames[!ok],collapse=","),
               "not of class sta"))
  if (length(stas)>1) { # join existing sta objects
    ss=stas[[1]]$ss
    slevs=attr(ss,"slevs"); dlevs=attr(ss,"dlevs")
    nt=attr(ss,"nt"); dsn=ss[[1]]$d$dsn
    cat(paste("Object",stanames[1],
              "has state levels",paste(slevs,collapse=","),
              "dimension levels",paste(dlevs,collapse=","),
              "and",nt,"trace levels, with baseline design",dsn,"\n"))
    snum=length(ss)
    for (i in 2:length(stas)) {
      slevsi=attr(stas[[i]]$ss,"slevs"); dlevsi=attr(stas[[i]]$ss,"dlevs")
      nti=attr(stas[[i]]$ss,"nt"); dsni=stas[[i]]$ss[[1]]$d$dsn
      cat(paste("Object",stanames[i],
                "has state levels",paste(slevsi,collapse=","),
                "dimension levels",paste(dlevsi,collapse=","),
                "and",nti,"trace levels, with baseline design",dsni,"\n"))
      if (!all(slevs==slevsi) || !all(dlevs==dlevsi) || nt!=nti || dsn!=dsni)
          stop("Design of this sta object does not match previous designs")
      for (j in 1:length(stas[[i]]$ss)) {
        snum=snum+1
        ss[[as.character(snum)]]=stas[[i]]$ss[[j]]
      }
    }
    dss=unlist(lapply(ss,function(x){attr(x,"datasource")}))
    dups=duplicated(dss)
    if ( any(dups) ) {
      cat("Data sources are duplicated for the following participants:\n")
      print(dss[dups])
    }
    sta=list(all=list(p=list(e=NULL,t=NULL,m=NULL)),ss=ss)
    class(sta)="sta"
  } else if (checkConverge) {
    sta=stas[[1]]
    require(coda)
    lmcmcl=lapply(sta$ss,getmcmcl,nmcmc=nmcmc)
    cat("\nGelman Diagnostic for Convergence (multivariate)\n")
    print(round(unlist(lapply(lmcmcl,function(x){attr(x,"mR")})),2))
    badnmcmc=unlist(lapply(lmcmcl,function(x){attr(x,"badnmcmc")}))
    isbad=badnmcmc!=0
    if (any(isbad)) {
      cat("\nSamples ignored (nmcmc not a multiple of stored samples):\n")
      print(badnmcmc[isbad])
    }
    if (is.character(mcmcSave) && nchar(mcmcSave[1])>0) {
      if (all(unlist(lapply(lmcmcl,function(x){any(is.na(unlist(x)))}))))
        cat(paste("\nTrace samples not saved as insufficient stored samples",
                       "for any participant\n"))  else
      {
        cat(paste("\nList of trace sample mcmc.list objects for each",
                  "participant saved as",mcmcSave[1],"\n"))
        assign(mcmcSave,lmcmcl,envir=.GlobalEnv)
      }
    }
  } else { # manage stored samples
    sta=stas[[1]]
    slevs=attr(sta$ss,"slevs"); dlevs=attr(sta$ss,"dlevs")
    nt=attr(sta$ss,"nt"); dsn=sta$ss[[1]]$d$dsn
    cat(paste("Object",stanames[1],
              "has state levels",paste(slevs,collapse=","),
              "dimension levels",paste(dlevs,collapse=","),
              "and",nt,"trace levels, with baseline design",dsn,"\n\n"))
    if ( !is.integer(nkeepe) && nkeepe>=0 ) {
      for (i in names(sta$ss)) {
        if (nkeepe==0 || is.null(sta$ss[[i]]$p$e$s) )
          sta$ss[[i]]$p$e$s=NULL else
          if (nkeepe<dim(sta$ss[[i]]$p$e$s[[1]])[2]) {
            sta$ss[[i]]$p$e$s[[1]]=sta$ss[[i]]$p$e$s[[1]][,1:nkeepe]
            sta$ss[[i]]$p$e$s[[2]]=sta$ss[[i]]$p$e$s[[2]][,1:nkeepe]
          }
      }
      cat("Stored encompassing samples now of sizes:\n")
      print(unlist(lapply(sta$ss,function(x){dim(x$p$e$s[[1]])[2]})))
    } else stop("nkeepe must be a non-negative integer")
    if ( !is.integer(nkeept) && nkeept>=0 ) {
      for (i in names(sta$ss)) {
        if (nkeept==0 || is.null(sta$ss[[i]]$p$t$s) )
          sta$ss[[i]]$p$t$s=NULL else
          if (nkeept<dim(sta$ss[[i]]$p$t$s)[3])
            sta$ss[[i]]$p$t$s=sta$ss[[i]]$p$t$s[,,1:nkeept]
      }
      cat("\nStored trace samples now of sizes:\n")
      print(unlist(lapply(sta$ss,function(x){dim(x$p$t$s)[3]})))
    } else stop("nkeept must be a non-negative integer")
    if ( !is.integer(nkeepm) && nkeepm>=0 ) {
      for (i in names(sta$ss)) {
        if (nkeepm==0 || is.null(sta$ss[[i]]$p$m$s) )
          sta$ss[[i]]$p$m$s=NULL else {
          if (keepbestm) {
            nm=sta$ss[[i]]$m$m$control$nm
            bestm=dimnames(sta$ss[[i]]$p$m$s)[[3]]==names(nm)[which.max(nm)]
            sta$ss[[i]]$p$m$s=sta$ss[[i]]$p$m$s[,,bestm]
          }
          if (nkeepm<dim(sta$ss[[i]]$p$m$s)[3])
            sta$ss[[i]]$p$m$s=sta$ss[[i]]$p$m$s[,,1:nkeepm]
        }
      }
      cat("\nStored monotonic samples now of sizes:\n")
      print(unlist(lapply(sta$ss,function(x){dim(x$p$m$s)[3]})))
    } else stop("nkeepm must be a non-negative integer")
  }
  assign(staname,sta,envir=.GlobalEnv)
  paste("Object",staname,"has been created/updated.")
}

################################################################################
############################### SAMPLING  ######################################
################################################################################

# verbose:
#   2 gives per participant timing,
#   1 gives model timing
#   0 gives no timing info
stSample=function(bosname="",
                  refresh=FALSE,maxt=0,ci=95,
                  BFe=TRUE,sampe=1e5,ed=.0005,nkeepe=1e4,
                  BFt=TRUE,sampt=5e3,burn=100,td=.005,nkeept=1e4,nkeepm=Inf,
                  verbose=1) {

  ##############################################################################
  ########################   ENCOMPASSING MODEL  ###############################
  ##############################################################################

  bosBFt <- function(bos,maxt=0,refresh=FALSE,ci=.95,d=5e-3,
                     sampe=1e5,nkeepe=1e4,
                     verbose=2) {
  # Input=sta object for one experiment, "bos"
  # Output=updated version of input
  #
  # sampe=sample size per pass
  # d=precision of BF estimate with confidence ci
  #
  # Two modes
  # refresh=T: used to update control parameters
  #            given d, ci (used when these changed to update timing)
  #            Fast (i.e., no sampling) UNLESS no sampling done yet, in which
  #            case does one pass to get necessary timing info
  # refresh=FALSE: Runs at most for time=maxt (hours). If maxt=0 just one pass
  #
  # verbose: 0=silent, 1=no output for each participant,
  #          2=output for each participant
  # NOTE: d and ci: can all be specified as a scalar or vector with one
  #       entry per participant

    # Does main compute for one participant
    # Input=sta object for one participant, bo, ouput is updated version
    boBFt <- function(bo,refresh=FALSE,maxt=0,
                      ci=.95,d=5e-3,sampe=1e5,nkeepe=1e4,
                      verbose=2) {

      # Get encompasing samples
      encompsamp <- function(o,s1,s2,nci,sampe) {
        d=matrix(nrow=nci,ncol=sampe)
        d[o[[1]][1],]=rbeta(sampe,s1[o[[1]][1]],s2[o[[1]][1]])
        for (j in 1:length(o))
          for (i in 2:length(o[[j]]))
            d[o[[j]][i],]=rbeta(sampe,s1[o[[j]][i]],s2[o[[j]][i]])
        d
      }

      # Returns number of columns in esamp sample that fulfill
      # conditions specified in o (local version of inorder)
      ino <- function(d,o) {
        di=1:dim(d)[2]
        for (j in 1:length(o)) {
          for (i in 2:length(o[[j]])) {
            di=di[d[o[[j]][i-1],di] <= d[o[[j]][i],di]]
            if (length(di)==0) break
          }
          if (length(di)==0) break
        }
        length(di)
      }

      # credible interval of width ci for posterior proportion
      calcci=function(ci,n,N) {
        c(qbeta(1-(1-ci)/2,n+1,N-n+1),
          qbeta((1-ci)/2,n+1,N-n+1))
      }

      # Updates control confidence intervals
      updateci=function(ctrl,dsn) {
        nchain=length(ctrl$doned)
        for (i in 1:nchain) { # Update CIs
          ctrl$cis[i,]=calcci(ctrl$ci,ctrl$ntraces[i],ctrl$Ntraces[i])[2:1]
        }
        ctrl$dcis=apply(ctrl$cis,1,diff)
        ctrl$doned=(ctrl$dcis<ctrl$d)
        ctrl
      }

      # Updates control timing estimates
      updatetime=function(ctrl,dsn,sampe) {
        # finds N necessary to get a width ci credible inverval = d
        getn=function(ci,n,d) {
          # divergence from desired ci=d, used in line search
          sfun=function(N,n,ci,d) {
            N=round(N)
            abs(log(d/-diff(calcci(ci,n,N))))
          }
          optimize(sfun,c(1e5,1e14),n=n,ci=ci,d=d)$minimum
        }
        nchain=length(ctrl$doned)
        for (i in 1:nchain) { # Figure out timing
          if (ctrl$doned[i]) {ctrl$moren[i]=0; ctrl$moret[i]=0}
          else {
            n=ctrl$ntraces[i]
            dval=ctrl$d
            ctrl$moren[i]=max(getn(ci,n,dval)-ctrl$Ntraces[i],sampe)
            ctrl$moret[i]=ctrl$tpers*ctrl$moren[i]/60
          }
        }
        ctrl
      }

      # MAIN BODY of boBFt
      #Initialize
      dsn=bo$d; nchain=nrow(dsn$s1)
      nc=apply(dsn$s1,1,function(x){sum(!is.na(x))})
      ctrl=bo$m$t$control
      if (is.null(ctrl)) { # First time
         ctrl=list(models="M1=trace, M2=encompassing",
                   d=d, doned=rep(FALSE,nchain),ci=ci,
                   tpers=NA,moret=NA,moren=NA,
                   cis=matrix(nrow=nchain,ncol=2),dcis=rep(NA,nchain),
                   ntraces=rep(0,nchain),Ntraces=rep(0,nchain))
         if (refresh) {maxt=0; refresh=FALSE; onepass=TRUE}
      } else {             # Load in new ci, badp and d
          ctrl$d=d; ctrl$ci=ci; onepass=FALSE
          ctrl=updateci(ctrl,dsn)
      }
      if ( !refresh && !all(ctrl$doned)) {   # Havent given up as BF<badp
        # Initial loop to get timing for currant run and/or add graph
        ctrl$tpers=system.time({
          if ( is.null(bo$p$e$s) )
            bo$p$e$s=vector(mode="list",length=nchain)
          for (i in 1:nchain) {
            if ( !ctrl$doned[i] ) {
              s1=dsn$s1[i,]; s2=dsn$s2[i,]
              esamp=encompsamp(dsn$odr[[i]],s1,s2,nc[i],sampe)
              nget=min(c(nkeepe-attr(bo$p$e$s,"nsamp"),dim(esamp)[2]))
              if ( is.null(attr(bo$p$e$s,"nsamp")) ) {
                bo$p$e$s[[i]]=esamp[,1:nget]
              } else if (nget>0) {
                bo$p$e$s[[i]]=cbind(bo$p$e$s[[i]],esamp[,1:nget])
              }
              ctrl$ntraces[i]=ctrl$ntraces[i]+ino(esamp,dsn$odr[[i]])
              ctrl$Ntraces[i]=ctrl$Ntraces[i]+sampe
            }
          }
          attr(bo$p$e$s,"nsamp")=dim(bo$p$e$s[[1]])[2]
        })[3]/(sum(!ctrl$doned)*sampe)
        ctrl=updateci(ctrl,dsn)
        if  (onepass || all(ctrl$doned)) maxpass=0
          else maxpass=floor(maxt/(sampe*ctrl$tpers/60))
        pass=0
        repeat {
          pass=pass+1
          if ((pass>maxpass) | all(ctrl$doned)) break
          else {
            for (i in 1:nchain) {
              if (!ctrl$doned[i]) {
                s1=dsn$s1[i,]; s2=dsn$s2[i,]
                esamp=encompsamp(dsn$odr[[i]],s1,s2,nc[i],sampe)
                nget=min(c(nkeepe-attr(bo$p$e$s,"nsamp"),dim(esamp)[2]))
                if (nget>0) {
                  bo$p$e$s[[i]]=cbind(bo$p$e$s[[i]],esamp[,1:nget])
                }
                ctrl$ntraces[i]=ctrl$ntraces[i]+ino(esamp,dsn$odr[[i]])
                ctrl$Ntraces[i]=ctrl$Ntraces[i]+sampe
              }
            }
            ctrl=updateci(ctrl,dsn)
          }
        }
      }
      ctrl=updatetime(ctrl,dsn,sampe)
      pvals=" (Chain CI-Size: "
      for (i in 1:nchain)
        pvals=paste(pvals,i,"=",round(ctrl$dcis[i],4)," ",sep="")
      if (verbose>1) {
        cat(paste("Participant ",bo$d$sname,pvals,"): ",
                  round(sum(ctrl$moret),2),"\n",sep=""))
      }
      bo$m$t$control=ctrl
      bo
    } # END of boBFt

    # MAIN BODY of bosBFt
    ns=length(bos$ss)
    if (length(d)==1) d=rep(d,ns)
    if (length(ci)==1) ci=rep(ci,ns)
    if (verbose>0) {
     cat("\nUPDATING REMAINING TIME FOR ENCOMPASSING MODEL SAMPLING (minutes)\n")
    }
    for (i in 1:ns) { # update to refresh tmod$control ok, doned, moret and moren
      firsttime=is.null(bos$ss[[i]]$m$t$control)
      bos$ss[[i]]=boBFt(bos$ss[[i]],refresh=TRUE,d=d[i],ci=ci[i],
                        sampe=sampe,verbose=verbose,nkeepe=nkeepe)
      if (firsttime) bos$ss[[i]]$m$t$control$doned=FALSE
    }
    if (verbose>0)  cat(paste("TOTAL TIME REMAINING FOR ENCOMPASSING SAMPLING:",
                        round(sum(unlist(lapply(bos$ss,function(x){
                          x$m$t$control$moret
                        }))),2),"\n")
                    )
    if ( !refresh ) {
      todo=unlist(lapply(bos$ss,function(x){any(!x$m$t$control$doned)}))
      if (any(todo)) {
        tmax=maxt/sum(todo)
        if (verbose>0)
           cat("\nWORKING ON ENCOMPASSING MODEL FOR REMAINING PARTICIPANTS\n")
         repeat {
            tdone=system.time({
               for (i in 1:length(bos$ss)) {
                  if (todo[i])
                     bos$ss[[i]]=boBFt(bos$ss[[i]],maxt=tmax,d=d[i],
                        ci=ci[i],sampe=sampe,verbose=verbose,nkeepe=nkeepe)
               }
            })[3]/60
            tmax=tmax-tdone
            todo=unlist(lapply(bos$ss,function(x){any(!x$m$t$control$doned)}))
            if (!any(todo) || (tmax<=0)) break
         }
      }
    }
    bos
  } # END of bosBFt


  ##############################################################################
  ############################  TRACE MODEL ####################################
  ##############################################################################


  bosBFm <- function(bos,maxt=0,refresh=FALSE,
                     ci=.95,d=5e-3,burn=100,
                     sampt=1e4,nkeept=1e4,nkeepm=1e4,
                     verbose=2) {
  # Subset of bosBFt parameters, except
  #   sampt=sample size per pass
  #   burn=number of burnin samples

    # Input=BOAST object for one participant, bo
    boBFm <- function(bo,sampt=1e4,burn=100,d=1e-3,ci=.95,firsttime=FALSE,
                      refresh=FALSE,maxt=0,verbose=2,nkeept=1e4,nkeepm=Inf) {

      table1d=function(tab,dsn,sampt=1e4,burn=100) {
      # creates table of counts for each monotonic model
      # Used by bosBFm and by grid version

        mtab=function(res,tab,dsn) {
        # keeps only monotinic samples in res (adding names giving their ranks)
        # and updates the count of each type in tab

          # updates tab with counts of monotonic model orders in res
          addtable=function(tab,res) {
            tab1=table(dimnames(res)[[1]])
            tnam1=names(tab1)
            for (i in 1:length(tab1)) {
              isin=tnam1[i]==names(tab)
              if (any(isin)) {
                n=tab[isin]+tab1[i]
                tab[isin]=tab[isin]+tab1[i]
              } else {
                tab=c(tab1[i],tab)
              }
            }
            tab
          }

          # MAIN BODY of mtab
          idx=1:dim(res)[1]; bmat=matrix(nrow=dim(res)[1],ncol=2)
          if (dsn$dsn=="B0") {            # B0 design
            nt=dim(res)[2]
            resa=array(res,dim=c(dim(res)[1],2*nt,2))
          } else {                        # B2 design
            nt=(dim(res)[2]-1)/2
            resa=res[,-1,]                # remove FA
          }
          for (j in 1:nt) for (k in (nt+1):(2*nt)){
            bmat[idx,1]=resa[idx,j,1]<resa[idx,k,1]
            bmat[idx,2]=resa[idx,j,2]<resa[idx,k,2]
            idx=idx[bmat[idx,1]==bmat[idx,2]]
            if (length(idx)==0) break
          }
          if (length(idx)==0)
            list(mtab=tab,mres=array(dim=c(0,dim(res)[2:3]))) else {
            if (dsn$dsn=="B0") res=resa
            if (length(idx)==1) {
              rdim=dim(res)
              res=res[idx,,]
              idx=paste(rank(as.vector(res[,1]),ties.method="random"),
                        collapse="")
              res=array(res,dim=c(1,rdim[2:3]))
            } else {
              res=res[idx,,]
              idx=apply(res[,,1],1,function(x){
                paste(rank(as.vector(x),ties.method="random"),collapse="")
              })
            }
            dimnames(res)=list(idx,dimnames(res)[[2]],dimnames(res)[[3]])
            list(mtab=addtable(tab,res),mres=res)
          }
        } # END of mtab

        # Samples one trace at a time (with or without single baseline)
        gibbst <- function(tab,dsn,nsamp=1e4,burn=100) {
          s1s=dsn$s1; s2s=dsn$s2
          ntrace=dim(s1s)[1]
          nsamp=nsamp+burn+1; n=length(s1s[1,])
          out <- array(dim=c(nsamp,n,ntrace))
          for (t in 1:ntrace) {
            out[1,,t] <- c(1:n+1)/(n+2) # monotonic initial sample
            s1=s1s[t,]; s2=s2s[t,]
            for (i in 2:(nsamp)) {
              out[i,,t]=out[i-1,,t]
              ki <- trunc(runif(n)*n)+1
              u <- runif(n)
              for (k in 1:n) {
                if (ki[k]==1) plow=0 else
                  plow=pbeta(out[i,ki[k]-1,t],s1[ki[k]],s2[ki[k]]) # lower bound
                if (ki[k]==n) phi=1  else
                  phi =pbeta(out[i,ki[k]+1,t],s1[ki[k]],s2[ki[k]]) # upper bound
                # inverse probability sample
                out[i,ki[k],t] <-
                  qbeta(plow+u[k]*(phi-plow),s1[ki[k]],s2[ki[k]])
              }
            }
          }
          res=out[-c(1:(burn+1)),,]
          outs=mtab(res,tab,dsn)
          outs$res=res
          outs
        }

        # Samples 2 traces with same baseline less than 1st value of each trace
        gibbs2t <- function(tab,dsn,nsamp=1e4,burn=100) {
          s1s=dsn$s1; s2s=dsn$s2
          ntrace=dim(s1s)[1]; nsamp=nsamp+burn+1; n=length(s1s[1,])
          out <- array(dim=c(nsamp,n,ntrace)); n1=(n-1)/2; n2=n1
          for (t in 1:ntrace) {
            # arbitary initial sample
            out[1,,t] <- c(c(1:(n1+1))/(n1+2),c(2:(n2+1))/(n2+2))
            s1=s1s[t,]; s2=s2s[t,]
            for (i in 2:(nsamp)) {
              out[i,,t]=out[i-1,,t]
              ki <- trunc(runif(n)*n)+1
              u <- runif(n)
              for (k in 1:n) {
                # lower bound
                if (ki[k]==1) plow=0 else
                  if (ki[k]==(n1+2))
                    plow=pbeta(out[i,1,t],s1[ki[k]],s2[ki[k]]) else
                    plow=pbeta(out[i,ki[k]-1,t],s1[ki[k]],s2[ki[k]])
                # upper bound
                if ((ki[k]==(n1+1)) || (ki[k]==n)) phi=1 else
                  if (ki[k]==1) {
                    phi = min(pbeta(out[i,2,t],s1[ki[k]],s2[ki[k]]),
                              pbeta(out[i,(n1+2),t],s1[ki[k]],s2[ki[k]]))
                  } else phi = pbeta(out[i,ki[k]+1,t],s1[ki[k]],s2[ki[k]])
                # inverse probability sample
                out[i,ki[k],t] <-
                  qbeta(plow+u[k]*(phi-plow),s1[ki[k]],s2[ki[k]])
              }
            }
          }
          res=out[-c(1:(burn+1)),,]
          outs=mtab(res,tab,dsn)
          outs$res=res
          outs
        }

        # MAIN BODY of table1D
        if (dsn$dsn=="B0") tmp=gibbst(tab,dsn,nsamp=sampt,burn=burn)
        if (dsn$dsn=="B2") tmp=gibbs2t(tab,dsn,nsamp=sampt,burn=burn)
        tmp
      } # END of table1D


      # credible interval of width ci proportions of no overlap,
      # monotonic and non-monotonic
      calccim=function(ci,n,N) {
        ciu=1-(1-ci)/2; cil=(1-ci)/2
        rbind(c(qbeta(ciu,n[1]+1,N-n[1]+1),qbeta(cil,n[1]+1,N-n[1]+1)),
                c(qbeta(ciu,n[2]+1,N-n[2]+1),qbeta(cil,n[2]+1,N-n[2]+1)),
                c(qbeta(ciu,n[3]+1,N-n[3]+1),qbeta(cil,n[3]+1,N-n[3]+1))
        )
      }

      updatecim=function(ctrl) {
        n=c(ctrl$n1d,ctrl$Nm-sum(ctrl$n1d))
        ctrl$cis=calccim(ctrl$ci,n,ctrl$Nm)[,2:1]
        ctrl$dcis=apply(ctrl$cis,1,diff)
        ctrl$doned=all(ctrl$dcis<ctrl$d)
        ctrl
      }

      # time to get to required posterior proportion accuracy
      updatetimem=function(ctrl,sampt) {

        # finds N necessary to get a width ci credible inverval = d
        getnm=function(ci,n,d,sampt) {
          # divergence from desired ci=d, used in line search
          sfun=function(N,n,ci,pp1d,d,cinum) {
            N=round(N)
            abs(log(d/-diff(calccim(ci,n,N)[cinum,])))
          }
          moren1=optimize(sfun,c(sampt,1e14),n=n,ci=ci,d=d,cinum=1)$minimum
          moren2=optimize(sfun,c(sampt,1e14),n=n,ci=ci,d=d,cinum=2)$minimum
          moren3=optimize(sfun,c(sampt,1e14),n=n,ci=ci,d=d,cinum=3)$minimum
          max(c(moren1,moren2,moren3))
        }

        if (ctrl$doned) {ctrl$moren=0; ctrl$moret=0}
        else {
            n=c(ctrl$n1d,ctrl$Nm-sum(ctrl$n1d))
            dval=ctrl$d
            ctrl$moren=max(getnm(ctrl$ci,n,dval,sampt)-ctrl$Nm,sampt)
            ctrl$moret=ctrl$tpers*ctrl$moren/60
        }
        ctrl
      }

      lap=function(tab) {
      # returns vector of counts for No Overlap and Overlap models
      # Used by boBFm
        tnams=names(tab)
        if (is.null(tnams)) tmp=c(0,0) else {
          n=nchar(tnams[1]); n2=round(n/2)
          nn1=paste(1:n,collapse="")
          nn2=paste(c((n2+1):n,1:n2),collapse="")
          is1=tnams==nn1; is2=tnams==nn2
          if (any(is1)) nl=tab[is1] else nl=0
          if (any(is2)) nl=nl+tab[is2]
          tmp=c(nl,sum(tab)-nl)
        }
        names(tmp)=c("NoLap","OLap")
        tmp
      }

      # MAIN BODY of boBFm
      dsn=bo$d
      ctrl=bo$m$m$control
      if ( is.null(ctrl) ) { # First time
          nt=dim(dsn$s1)[2]; if (dsn$dsn=="B2") {nt=nt-1; nt=round(nt/2)}
          pp1d=((factorial(nt))^-2)*factorial(2*nt); pp1d=c(2,pp1d-2)/(pp1d^2)
          ctrl=list(models="M1=1D, M2=MultiD",
                    d=d,doned=FALSE,ci=ci,
                    tpers=NA,moret=NA,moren=NA,
                    cis=NA,dcis=NA,
                    pp1d=pp1d,
                    nm=NULL,n1d=c(0,0),Nm=0)
          if (refresh) {maxt=0; refresh=FALSE; onepass=TRUE}
      } else {
        ctrl$d=d; ctrl$ci=ci; onepass=FALSE
        ctrl=updatecim(ctrl)
        if (firsttime) ctrl$doned=FALSE
      }
      if (!refresh && !ctrl$doned) {
        ctrl$tpers=system.time({
          mtrr=table1d(ctrl$nm,dsn,sampt=sampt,burn=burn)
          ctrl$nm=mtrr[[1]]
          ctrl$n1d=lap(ctrl$nm)
          ctrl$Nm=ctrl$Nm+sampt
          # update stored trace samples
          if ( is.null(bo$p$t$s) ) {
            nget=min(c(nkeept,dim(mtrr[[3]])[1]))
            if (nget>0) {
              bo$p$t$s=aperm(mtrr[[3]][1:nget,,,drop=FALSE],c(2,3,1))
              attr(bo$p$t$s,"nsamp")=nget
            }
          } else {
            nget=min(c(nkeept-attr(bo$p$t$s,"nsamp"),dim(mtrr[[3]])[1]))
            if (nget>0) {
              if (nget==1) {
                tmp=mtrr[[3]][1,,]
                tmp=array(tmp,dim=c(1,dim(tmp)))
              } else tmp= mtrr[[3]][1:nget,,]
              nin=dim(bo$p$t$s)[3]+nget
              bo$p$t$s=array(c(bo$p$t$s,aperm(tmp,c(2,3,1))),
                             dim=c(dim(bo$p$t$s)[1:2],nin))
              attr(bo$p$t$s,"nsamp")=nin
            }
          }
          # update stored monotonic samples
          if ( is.null(bo$p$m$s) ) {
            nget=min(c(nkeepm,dim(mtrr[[2]])[1]))
            if (nget>0) {
            	bo$p$m$s=aperm(mtrr[[2]][1:nget,,,drop=FALSE],c(2,3,1))
                attr(bo$p$m$s,"nsamp")=nget
            }
          } else {
            nget=min(c(nkeepm-attr(bo$p$m$s,"nsamp"),dim(mtrr[[2]])[1]))
            if (nget>0) {
              if (nget==1) {
                tmp=mtrr[[2]][1,,]
                tmp=array(tmp,dim=c(1,dim(tmp)))
                dimnames(tmp)=list(dimnames(mtrr[[2]])[[1]][1],NULL,NULL)
              } else tmp=mtrr[[2]][1:nget,,]
              nin=dim(bo$p$m$s)[3]+nget
              newnams=c(dimnames(bo$p$m$s)[[3]],dimnames(tmp)[[1]])
              bo$p$m$s=array(c(bo$p$m$s,aperm(tmp,c(2,3,1))),
                             dim=c(dim(bo$p$m$s)[1:2],nin))
              dimnames(bo$p$m$s)[[3]]=newnams
              attr(bo$p$m$s,"nsamp")=nin
            }
          }
        })[3]/sampt
        ctrl=updatecim(ctrl)
        if (onepass || ctrl$doned) maxpass=0 else
                                   maxpass=floor(maxt/(ctrl$tpers*sampt/60))
        pass=0
        repeat {
          pass=pass+1
          if ((pass>maxpass) | !ctrl$doned ) break else {
            mtrr=table1d(ctrl$nm,dsn,sampt=sampt,burn=burn)
            ctrl$nm=mtrr[[1]]
            ctrl$n1d=lap(ctrl$nm)
            ctrl$Nm=ctrl$Nm+sampt
            nget=min(c(nkeept-attr(bo$p$t$s,"nsamp"),dim(mtrr[[3]])[1]))
            if (nget>0) {
              if (nget==1) {
                tmp=mtrr[[3]][1,,]
                tmp=array(tmp,dim=c(1,dim(tmp)))
              } else tmp= mtrr[[3]][1:nget,,]
              nin=dim(bo$p$t$s)[3]+nget
              bo$p$t$s=array(c(bo$p$t$s,aperm(tmp,c(2,3,1))),
                             dim=c(dim(bo$p$t$s)[1:2],nin))
              attr(bo$p$t$s,"nsamp")=nin
            }
            nget=min(c(nkeepm-attr(bo$p$m$s,"nsamp"),dim(mtrr[[2]])[1]))
            if (nget>0) {
              if (nget==1) {
                tmp=mtrr[[2]][1,,]
                tmp=array(tmp,dim=c(1,dim(tmp)))
              } else tmp= mtrr[[2]][1:nget,,]
              newnams=c(dimnames(bo$p$m$s)[[3]],dimnames(tmp)[[1]])
              nin=dim(bo$p$tms)[3]+nget
              bo$p$m$s=array(c(bo$p$m$s,aperm(tmp,c(2,3,1))),
                             dim=c(dim(bo$p$m$s)[1:2],nin))
              dimnames(bo$p$m$s)[[3]]=newnams
              attr(bo$p$m$s,"nsamp")=nin
            }
            ctrl=updateci(ctrl)
          }
        }
      }
      ctrl=updatetimem(ctrl,sampt)
      bo$m$m$control=ctrl
      pvals=paste("(CI-Size: NoOverlap=",round(ctrl$dcis[1],4),
             " 1D=",round(ctrl$dcis[2],4),
             " MD=",round(ctrl$dcis[3],4),")",sep="")
      if (verbose>1)
        cat(paste("Participant ",bo$d$sname,pvals,": ",
                  round(sum(ctrl$moret),2),"\n",sep=""))
      bo
    } # END of boBFm

    # MAIN BODY of bosBFm
    firsttime=unlist(lapply(bos$ss,function(x){is.null(x$m$m$control)}))
    ns=length(bos$ss)
    if (length(d)==1) d=rep(d,ns)
    if (length(ci)==1) ci=rep(ci,ns)
    if (verbose>0)
      cat("\nUPDATING REMAINING TIME FOR TRACE MODEL SAMPLING (minutes)\n")
    for (i in 1:ns)  # update to refresh tmod$control: doned, moret and moren
      bos$ss[[i]]=boBFm(bos$ss[[i]],refresh=TRUE,d=d[i],ci=ci[i],sampt=sampt,
                        burn=burn,verbose=verbose,nkeept=nkeept,nkeepm=nkeepm)
    if (verbose>0) cat(paste("TOTAL TIME REMAINING FOR TRACE SAMPLING:",
                        round(sum(unlist(lapply(bos$ss,function(x){
                          x$m$m$control$moret
                        }))),2),"\n"))
    if ( !refresh ) {
      todo=unlist(lapply(bos$ss,function(x){!x$m$m$control$doned}))
      todo = todo | firsttime
      if (any(todo)) {
        tmax=maxt/sum(todo)
        if (verbose>0)
           cat("\nWORKING ON TRACE MODEL FOR REMAINING PARTICIPANTS\n")
        repeat {
           tdone=system.time({
              for (i in 1:ns) {
                if (todo[i])
                  bos$ss[[i]]=boBFm(bos$ss[[i]],maxt=tmax,d=d[i],
                     ci=ci[i],sampt=sampt,burn=burn,firsttime=firsttime[i],
                     nkeept=nkeept,nkeepm=nkeepm,verbose=verbose)
              }
           })[3]/60
           tmax=tmax-tdone
           todo=unlist(lapply(bos$ss,function(x){!x$m$m$control$doned}))
           if (!any(todo) || (tmax<=0)) break
        }
      }
    }
    bos
  } # END of bosBFm

  # main body of stSample
  if (!is.character(bosname) || nchar(bosname)==0)
    stop("Name of sta object must be be given as quoted text")
  if (!exists(bosname,envir=.GlobalEnv))
    stop(paste("Object",bosname,"not present in global environment"))
  bos=get(bosname,envir=.GlobalEnv)
  if (class(bos)!="sta")
    stop(paste("Object",bosname,"does not have class \"sta\""))
  # calls to optimize in calccim estimating time remaining can cause warnings,
  # optimize handles sensibly so hide from user but be careful in debugging!
  op=options(); options(warn=-1)
  ci=ci/100    # convert ci to probability
  maxt=maxt*60 #convert maxt to minutes
  if (refresh) maxt=0
  repeat {
    t1=system.time({
      tt=maxt*BFe/(BFe+BFt); tm=maxt*BFt/(BFe+BFt)
      if (BFe) bos=bosBFt(bos,maxt=tt,refresh=refresh,
                          ci=ci,d=ed,sampe=sampe,nkeepe=nkeepe,
                          verbose=verbose)
      if (BFt) bos=bosBFm(bos,maxt=tm,refresh=refresh,ci=ci,d=td,
                          burn=burn,sampt=sampt,nkeept=nkeept,nkeepm=nkeepm,
                          verbose=verbose)
    })[3]/60
    if (verbose>0) {
      etime=sum(unlist(lapply(bos$ss,function(x){x$m$t$control$moret})))
      ttime=sum(unlist(lapply(bos$ss,function(x){x$m$m$control$moret})))
      cat(paste("\nTOTAL TIME REMAINING FOR ENCOMPASSING SAMPLING (minutes):",
                round(etime,2),"\n"))
      cat(paste("TOTAL TIME REMAINING FOR TRACE SAMPLING (minutes):",
                round(ttime,2),"\n"))
      cat(paste("TOTAL TIME REMAINING FOR ALL SAMPLING (minutes):",
                round(etime+ttime,2),"\n\n"))
    }
    maxt=maxt-t1
    alldone=all(unlist(lapply(bos$ss,function(x){
       c(x$m$m$control$doned,x$m$t$control$doned)
    })))
    if (maxt<=0 || (alldone)) break
  }
  options(op)
  assign(bosname,bos,envir=.GlobalEnv)
  paste("Object",bosname,"has been updated.")
}


################################################################################
############################### SUMMARY   ######################################
################################################################################

#  pmpBootav=function(pmp,nbsamp=1e4) {
#    ns=dim(pmp)[1]
#    samps=matrix(nrow=nbsamp,ncol=4)
#    for (i in 1:nbsamp) {
#       dat=
#       samps[i,]=apply(apply(
#        pmp[sample.int(ns,replace=T),],1,function(x){rmultinom(1,1,x)}),1,sum
#       )/ns
#    }
#    mn=apply(samps,2,mean); names(mn)=dimnames(pmp)[[2]]
#    ci95=apply(samps,2,quantile,probs=c(.025,.975))
#    dimnames(ci95)[[2]]=names(mn)
#    list(mean=mn,ci95=ci95)
#  }

# summarised information from et
# exclude is a boolean vector indicating participants to be included
stSummary=function(bosname="",BF=FALSE,traceTrue=FALSE,eachs=TRUE,nsig=3,
                   exclude=NULL,sorts="",extras="",
                   guiarg=NULL) {   # GUI related

  et=function(bos,exclude=NULL) {
  ##############################################################################
  ####### et - statistics relating to encompassing and trace model samples #####
  ##############################################################################
  # Extracts prior and posterior statistics encompassing and trace samples
  # in sta object
  # LEVEL 1: e=results based on encompassing samples, t=results based on trace
  #            samples
  # LEVEL 2: prior, counts, bf and p
  # Level 3: *For bf and p only, prior and counts always for each participant*:
  #          all=combination over participants, s=for each participant,
  #
  # prior:  For e, trace model prior
  #         For t, mono=all monotonic, overlap=no-overlap and overlap (matrix).
  #         NB: Any one mono is overlap[1,]/2
  # counts: For e, n=trace samples, N=number of e samples
  #         For t, n=mono samples,  N=number of t samples, noverlap=no-overlap
  #         and overlap monos (matrix)
  #         nmono=each specific mono samples (participant name list of order
  #               name vectors)
  # For e: bf and p for trace model (e.g., bf$all and bf$s)
  # For t, bf and p (e.g., bf$all and bf$s) for:
  #        m=monotonic model list, m$m=mono rel to trace,
  #        m$nm=non-mono rel to trace, m$mnm=mono vs. non-mono
  #        o=overlap list, o$on=overlapping vs. mono non-overlap,
  #        o$ono=momoe overlapping vs. mono-nonoverlapping,
  #        o$t, relative to trace model: o$t$o=overlap, o$t$no=no overlap,
  #        o$m, relative to mono model:  o$m$o=overlap, o$m$no=no overlap
  #        d=dimensionality
  #
  # PRIMARY MODEL COMPARISONS (e.g., for p$$all)
  # 1) TRACE TEST  (trace vs. non-trace): e$p$all (trace comparison)
  # 2) OVERLAP TEST (mono overlap vs. no lap):
  #       t$p$all$o$ono (monotonic overlap comparison)
  # 3) DIMENSIONALITY TEST  (mono overlap vs. nom-mono) :
  #       t$p$all$d (dimensionality comparison)

    test4=function(bos,traceTrue=FALSE) {
      cnams=bos[[1]]$d$chain.names; nchain=length(cnams)
      # Encompassing samples
      prpsE=apply(matrix(unlist(lapply(bos,function(x){x$d$priorp})),
                  nrow=nchain),2,prod)
      nE=matrix(unlist(lapply(bos,function(x){x$m$t$control$ntraces})),
                nrow=nchain)
      NE=matrix(unlist(lapply(bos,function(x){x$m$t$control$Ntraces})),
                nrow=nchain)
      popsE=apply((nE+1)/(NE+2),2,prod)
      nt=1-popsE
      # Trace based priors
      prpsT=matrix(unlist(lapply(bos,function(x){x$m$m$control$pp1d})),nrow=2)
      prpsT=t(rbind(prpsT,1-apply(prpsT,2,sum)))
      dimnames(prpsT)=list(names(bos),c("NoLap","OLap","MD"))
      # Trace samples
      nT=matrix(unlist(lapply(bos,function(x){x$m$m$control$n1d})),nrow=2)
      NT=unlist(lapply(bos,function(x){x$m$m$control$Nm}))
      nT=t(rbind(nT,NT-apply(nT,2,sum)))
      popsT=(nT+1)/(NT+2); dimnames(popsT)=dimnames(prpsT)
      popsT=popsT/apply(popsT,1,sum)
      # Multiply by trace probability so all four on same scale
      prpsT=prpsT*prpsE; popsT=popsT*popsE
      prpsT=cbind(prpsT,NT=1-prpsE); popsT=cbind(popsT,NT=nt)
      # For each model among 4 models
      bfs=popsT/prpsT; dimnames(bfs)[[2]][2]="UD"
      bf=apply(log(bfs),2,sum)
      if ( traceTrue ) {
       	pmps=bfs[,-4]/apply(bfs[,-4],1,sum)
      	pmp=exp(bf[-4]-max(bf[-4]))/sum(exp(bf[-4]-max(bf[-4])))
      } else { # remove non-trace from calculation
     	pmps=bfs/apply(bfs,1,sum)
      	pmp=exp(bf-max(bf))/sum(exp(bf-max(bf)))
      }
      dimnames(pmps)[[2]][2]="UD"
      list(lbfs=log(bfs),pmps=pmps,lbf=bf,pmp=pmp)
    }

  # Extract results from encompassing samples

    # FOR EACH CHAIN gets trace model priorp, n, N & BF or for each participant
    eps=function(bos) {
      cnams=bos[[1]]$d$chain.names; nchain=length(cnams)
      pp=matrix(unlist(lapply(bos,function(x){x$d$priorp})),nrow=nchain)
      dimnames(pp)=list(cnams,names(bos))
      ns=matrix(unlist(lapply(bos,function(x){x$m$t$control$ntraces})),
                nrow=nchain)
      dimnames(ns)=dimnames(pp)
      Ns=matrix(unlist(lapply(bos,function(x){x$m$t$control$Ntraces})),
                nrow=nchain)
      dimnames(Ns)=dimnames(pp)
      list(pp=pp,ns=ns,Ns=Ns)
    }

  # Extract results from trace samples

    # monotonic model N, overall priorp (pp), n and N and priorp (pps) & n (ns)
    # for overlap and no-overlap and number of monotonic model samples of each
    # type for each participant
    tps=function(bos) {
      pps=matrix(unlist(lapply(bos,function(x){x$m$m$control$pp1d})),nrow=2)
      dimnames(pps)=list(c("NoLap","OLap"),names(bos))
      pp=apply(pps,2,sum)
      ns=matrix(unlist(lapply(bos,function(x){x$m$m$control$n1d})),nrow=2)
      dimnames(ns)=dimnames(pps)
      n=apply(ns,2,sum); N=unlist(lapply(bos,function(x){x$m$m$control$Nm}))
      names(N)=names(pp); names(n)=names(pp)
      nmono=lapply(bos,function(x){x$m$m$control$nm})
      list(pps=pps,pp=pp,ns=ns,n=n,N=N,nmono=nmono)
    }

    # main body of et
    bos=bos$ss
    if (!is.null(exclude)) {
      oknams=vector(length=0)
      ns=length(bos); nok=ns-length(exclude)
      ok=vector(mode="list",length=nok); oki=1 # participant names with exclude
      for (i in 1:ns) if (!any(i==exclude)) {
        ok[[oki]]=bos[[i]]; oki=oki+1; oknams=c(oknams,i)
      }
      names(ok)=oknams
      bos=ok
    }
    eraw=eps(bos); traw=tps(bos); bfs=t(((eraw$ns+1)/(eraw$Ns+2))/eraw$pp)
    list(e=list(prior=eraw$pp,count=list(n=eraw$ns,N=eraw$Ns),
                bfs=bfs,pmps=bfs/(1+bfs)),
         t=list(prior=list(mono=traw$pp,overlap=traw$pps),count=list(n=traw$n,
                  N=traw$N,noverlap=traw$ns,nmono=traw$nmono)),
         test=test4(bos,traceTrue))
  } # end of et

  # Main body of sumbos
  if (!is.character(bosname) || nchar(bosname)==0)
    stop("Name of sta object must be be given as quoted text")
  if (!exists(bosname,envir=.GlobalEnv))
    stop(paste("Object",bosname,"not present in global environment"))
  bos=get(bosname,envir=.GlobalEnv)
  if (class(bos)!="sta")
    stop(paste("Object",bosname,"does not have class \"sta\""))
  dsn=bos$ss[[1]]$d$dsn
  nt=nt=attr(bos$ss,"nt")
  dlevs=attr(bos$ss,"dlevs")
  # details for precision/ci and timings
  noe=unlist(lapply(bos$ss,function(x){is.null(x$m$t$control)}))
  not=unlist(lapply(bos$ss,function(x){is.null(x$m$m$control)}))
  if (any(noe))
    cat("One or more participants with no encompassing model samples.\n")
  if (any(not))
    cat("One or more participants with no trace model samples.\n")
  if ( any(noe) || any(not) )
    stop("Run stSample with refresh=T to fix this problem!")
  dnoe=unique(unlist(lapply(bos$ss,function(x){x$m$t$control$d})))
  dnot=unique(unlist(lapply(bos$ss,function(x){x$m$m$control$d})))
  cino=unique(unlist(lapply(bos$ss,function(x){x$m$t$control$ci})))*100
  if( all(unlist(lapply(bos$ss,function(x){x$m$t$control$doned}))) &&
      all(unlist(lapply(bos$ss,function(x){x$m$m$control$doned}))) ) {# finished
     cat(paste("\nSampling complete at ",cino,
      "% credible interval with precision based on:\nEncompassing samples=",
      dnoe," and Trace samples=",dnot,"\n",sep=""))
  } else {                                                        # not finished
     cat(paste("\nSampling not complete at ",cino,
      "% credible interval with precision based on:\nEncompassing samples=",
      dnoe," and Trace samples=",dnot,"\n",sep=""))
     cat(paste("Estimated time remaining: Encompassing model - ",
               signif(sum(unlist(lapply(bos$ss,function(x){
                  x$m$t$control$moret
               }))),3)," minutes; ",
               "Trace model - ",
               signif(sum(unlist(lapply(bos$ss,function(x){
                 x$m$m$control$moret
               }))),3)," minutes","\n",sep=""))
  }
  if(any(is.null(exclude))) exclude="" else
    exclude=as.numeric(unlist(strsplit(as.character(exclude)," ")))
  etbos=et(bos,exclude)
  if(any(is.na(extras))) extras=""
  if ( is.logical(extras) && all(extras) ) extras=c(
    "Data sources",
    "Trace model prior probability",
    "No-Overlap model prior probability",
    "Uni-dimensional model prior probability",
    "Multi-dimensional model prior probability",
    "Non-Trace model samples from encompassing model",
    "Total samples from encompassing model",
    "Non-overlapping monotonic samples from trace model",
    "Overlapping monotonic samples from trace model",
    "Non-monotonic samples from trace model",
    "Total samples from trace model",
    "Monotonic model counts"
  )
  cat(paste("\nNumber of participants=",length(etbos$e$bfs[,1]),"\n",sep=""))

  if ( !BF & traceTrue )
  	cat("\nTrace model assumed to be true\n")

  cat(paste("\nGroup level analysis",sep=""))
  if(!BF) {
    if ( traceTrue ) {
       tmp=etbos$test$pmp[c(1:3)]
       names(tmp)=c("No-Overlap","Uni-dimensional","Multi-dimensional")
          cat(paste(" (posterior model probabilities)\n",sep=""))
    } else {
       tmp=etbos$test$pmp[c(4,1:3)]
       names(tmp)=c("Non-Trace","No-Overlap","Uni-dimensional","Multi-dimensional")
          cat(paste(" (posterior model probabilities)\n",sep=""))
    }
    print(round(tmp,nsig)) ; cat(paste("\n",sep=""))
  } else {
    tmp=exp(etbos$test$lbf[c(4,1:3)])
    names(tmp)=c("Non-Trace","No-Overlap","Uni-dimensional","Multi-dimensional")
    cat(paste(" (Bayes Factors)\n",sep=""))
    print(signif(tmp,(nsig+1))) ; cat(paste("\n",sep=""))
  }
  if(eachs) {
    if (!BF) {  #p
      stype="(posterior model probabilities)"
      if (traceTrue)
         tmp=round(t(etbos$test$pmps)[c(1,2,3),],nsig) else
         tmp=round(t(etbos$test$pmps)[c(4,1,2,3),],nsig)
    } else {  #BF
      stype="(Bayes Factors)"
      tmp=round(exp(t(etbos$test$lbfs)[c(4,1,2,3),]),nsig)
    }

    if(any(is.na(sorts)) || (sorts=="Non-Trace model" & !BF & traceTrue))
       sorts="" else {
       if(!BF & traceTrue) {
          if(sorts=="No-Overlap model") tmp=tmp[,order(tmp[1,])]
          if(sorts=="Uni-dimensional model") tmp=tmp[,order(tmp[2,])]
          if(sorts=="Multi-dimensional model") tmp=tmp[,order(tmp[3,])]
       } else {
          if(sorts=="Non-Trace model") tmp=tmp[,order(tmp[1,])]
          if(sorts=="No-Overlap model") tmp=tmp[,order(tmp[2,])]
          if(sorts=="Uni-dimensional model") tmp=tmp[,order(tmp[3,])]
          if(sorts=="Multi-dimensional model") tmp=tmp[,order(tmp[4,])]
       }
    }
    if ( !BF & traceTrue)
       dimnames(tmp)[[1]]=
          c("No-Overlap","Uni-dimensional","Multi-dimensional") else
       dimnames(tmp)[[1]]=
          c("Non-Trace","No-Overlap","Uni-dimensional","Multi-dimensional")
    cat(paste("\nIndividual participant analysis ",stype,"\n",sep=""))
    print(tmp) ; cat(paste("\n",sep=""))
  }

  if (any(extras=="Data sources")) {
    cat(paste("Data sources\n",sep=""))
    print(unlist(lapply(bos$ss,function(x){attr(x,"datasource")})))
    cat("\n")
  }
  # non-trace = 1-prior for trace model = 1-etbos$e$prior
  if(any(extras=="Trace model prior probability")) {
    cat(paste("Trace model prior probability\n",sep=""))
    print(prod(etbos$e$prior[,1])); cat("\n")
  }
  # no-overlap = NoLap monotonic samples = etbos$t$prior$overlap[1,]
  if(any(extras=="No-Overlap model prior probability")) {
    cat(paste("No-Overlap model prior probability within Trace Model\n",sep=""))
    print(etbos$t$prior$overlap[1,1]); cat("\n")
  }
  # Uni-dimensional = OLap monotonic samples = etbos$t$prior$overlap[2,]
  if(any(extras=="Uni-dimensional model prior probability")) {
    cat(paste("Uni-dimensional model prior probability within Trace Model\n",
              sep=""))
    print(etbos$t$prior$overlap[2,1]); cat("\n")
  }
  # Multi-dimensional = non-monotonic samples = 1-etbos$t$prior$mono
  if(any(extras=="Multi-dimensional model prior probability")) {
    cat(paste("Multi-dimensional model prior probability within Trace Model\n",
              sep=""))
    print(1-etbos$t$prior$mono[1]); cat("\n")
  }
  cat("\n")
  if(any(extras=="Non-Trace model samples from encompassing model")) {
    cat(paste("Non-Trace model samples from encompassing model\n",sep=""))
    print(etbos$e$count$N-etbos$e$count$n); cat("\n")
  }
  if(any(extras=="Total samples from encompassing model")) {
    cat(paste("Total samples from encompassing model\n",sep=""))
    print(etbos$e$count$N); cat("\n")
  }
  cat("\n")
  if(any(extras=="Non-overlapping monotonic samples from trace model")) {
    cat(paste("Non-overlapping monotonic samples from trace model\n",sep=""))
    print(etbos$t$count$noverlap[1,]); cat("\n")
  }
  if(any(extras=="Overlapping monotonic samples from trace model")) {
    cat(paste("Overlapping monotonic samples from trace model\n",sep=""))
    print(etbos$t$count$noverlap[2,]); cat("\n")
  }
  if(any(extras=="Non-monotonic samples from trace model")) {
    cat(paste("Non-monotonic samples from trace model\n",sep=""))
    print(etbos$t$count$N-etbos$t$count$n); cat("\n")
  }
  if(any(extras=="Total samples from trace model")) {
    cat(paste("Total samples from trace model\n",sep=""))
    print(etbos$t$count$N); cat("\n")
  }
  if(any(extras=="Monotonic model counts")) {
    cat(paste("Monotonic model counts for each participant\n",sep=""))
    r1=1:nt+(dsn=="B2")
    r2=(nt+1):(2*nt)+(dsn=="B2")
    if (dsn=="B2") cat(paste("1 indicates the baseline condition\n"))
    cat(paste(paste(r1,collapse=","),
       "indicate trace levels of dimension level",dlevs[1],"\n"))
    cat(paste(paste(r2,collapse=","),
       "indicate trace levels of dimension level",dlevs[2],"\n\n"))
    tmp=etbos$t$count$nmono
    for (i in names(tmp)) tmp[[i]]=sort(tmp[[i]],decreasing=T)
    print(tmp); cat("\n")
  }
}

################################################################################
############################### plotting #######################################
################################################################################

# posterior model probability plot
stProbplot=function(bosname="",
                    exclude=NULL,maint="",maino="",mainu="",mainm="",
                    ylabt="p(Non-Trace)",ylabo="p(No-Overlap)",
                    ylabu="p(Uni-dimensional)",ylabm="p(Multi-dimensional)",
                    xlab="Sorted Participants",
                    symb=1,weakl=TRUE,strongl=FALSE,plines=FALSE,pnames=TRUE,
                    ymin=0,ymax=1,guiarg=NULL) {

  et=function(bos,exclude=NULL) {
  ##############################################################################
  ####### et - statistics relating to encompassing and trace model samples #####
  ##############################################################################
  # Extracts prior and posterior statistics encompassing and trace samples
  # in sta object
  # LEVEL 1: e=results based on encompassing samples, t=results based on trace
  #            samples
  # LEVEL 2: prior, counts, bf and p
  # Level 3: *For bf and p only, prior and counts always for each participant*:
  #          all=combination over participants, s=for each participant,
  #
  # prior:  For e, trace model prior
  #         For t, mono=all monotonic, overlap=no-overlap and overlap (matrix).
  #         NB: Any one mono is overlap[1,]/2
  # counts: For e, n=trace samples, N=number of e samples
  #         For t, n=mono samples,  N=number of t samples, noverlap=no-overlap
  #         and overlap monos (matrix)
  #         nmono=each specific mono samples (participant name list of order
  #               name vectors)
  # For e: bf and p for trace model (e.g., bf$all and bf$s)
  # For t, bf and p (e.g., bf$all and bf$s) for:
  #        m=monotonic model list, m$m=mono rel to trace,
  #        m$nm=non-mono rel to trace, m$mnm=mono vs. non-mono
  #        o=overlap list, o$on=overlapping vs. mono non-overlap,
  #        o$ono=momoe overlapping vs. mono-nonoverlapping,
  #        o$t, relative to trace model: o$t$o=overlap, o$t$no=no overlap,
  #        o$m, relative to mono model:  o$m$o=overlap, o$m$no=no overlap
  #        d=dimensionality
  #
  # PRIMARY MODEL COMPARISONS (e.g., for p$$all)
  # 1) TRACE TEST  (trace vs. non-trace): e$p$all (trace comparison)
  # 2) OVERLAP TEST (mono overlap vs. no lap):
  #       t$p$all$o$ono (monotonic overlap comparison)
  # 3) DIMENSIONALITY TEST  (mono overlap vs. nom-mono) :
  #       t$p$all$d (dimensionality comparison)

    test4=function(bos) {
      cnams=bos[[1]]$d$chain.names; nchain=length(cnams)
      # Encompassing samples
      prpsE=apply(matrix(unlist(lapply(bos,function(x){x$d$priorp})),
                  nrow=nchain),2,prod)
      nE=matrix(unlist(lapply(bos,function(x){x$m$t$control$ntraces})),
                nrow=nchain)
      NE=matrix(unlist(lapply(bos,function(x){x$m$t$control$Ntraces})),
                nrow=nchain)
      popsE=apply((nE+1)/(NE+2),2,prod)
      nt=1-popsE
      # Trace based priors
      prpsT=matrix(unlist(lapply(bos,function(x){x$m$m$control$pp1d})),nrow=2)
      prpsT=t(rbind(prpsT,1-apply(prpsT,2,sum)))
      dimnames(prpsT)=list(names(bos),c("NoLap","OLap","MD"))
      # Trace samples
      nT=matrix(unlist(lapply(bos,function(x){x$m$m$control$n1d})),nrow=2)
      NT=unlist(lapply(bos,function(x){x$m$m$control$Nm}))
      nT=t(rbind(nT,NT-apply(nT,2,sum)))
      popsT=(nT+1)/(NT+2); dimnames(popsT)=dimnames(prpsT)
      popsT=popsT/apply(popsT,1,sum)
      # Multiply by trace probability so all four on same scale
      prpsT=prpsT*prpsE; popsT=popsT*popsE
      prpsT=cbind(prpsT,NT=1-prpsE); popsT=cbind(popsT,NT=nt)
      # For each model among 4 models
      bfs=popsT/prpsT; dimnames(bfs)[[2]][2]="UD"
      pmps=bfs/apply(bfs,1,sum); dimnames(pmps)[[2]][2]="UD"
      bf=apply(log(bfs),2,sum)
      pmp=exp(bf-max(bf))/sum(exp(bf-max(bf)))
      list(lbfs=log(bfs),pmps=pmps,lbf=bf,pmp=pmp)
    }

  # Extract results from encompassing samples

    # FOR EACH CHAIN gets trace model priorp, n, N & BF or for each participant
    eps=function(bos) {
      cnams=bos[[1]]$d$chain.names; nchain=length(cnams)
      pp=matrix(unlist(lapply(bos,function(x){x$d$priorp})),nrow=nchain)
      dimnames(pp)=list(cnams,names(bos))
      ns=matrix(unlist(lapply(bos,function(x){x$m$t$control$ntraces})),
                nrow=nchain)
      dimnames(ns)=dimnames(pp)
      Ns=matrix(unlist(lapply(bos,function(x){x$m$t$control$Ntraces})),
                nrow=nchain)
      dimnames(Ns)=dimnames(pp)
      list(pp=pp,ns=ns,Ns=Ns)
    }

  # Extract results from trace samples

    # monotonic model N, overall priorp (pp), n and N and priorp (pps) & n (ns)
    # for overlap and no-overlap and number of monotonic model samples of each
    # type for each participant
    tps=function(bos) {
      pps=matrix(unlist(lapply(bos,function(x){x$m$m$control$pp1d})),nrow=2)
      dimnames(pps)=list(c("NoLap","OLap"),names(bos))
      pp=apply(pps,2,sum)
      ns=matrix(unlist(lapply(bos,function(x){x$m$m$control$n1d})),nrow=2)
      dimnames(ns)=dimnames(pps)
      n=apply(ns,2,sum); N=unlist(lapply(bos,function(x){x$m$m$control$Nm}))
      names(N)=names(pp); names(n)=names(pp)
      nmono=lapply(bos,function(x){x$m$m$control$nm})
      list(pps=pps,pp=pp,ns=ns,n=n,N=N,nmono=nmono)
    }

    # main body of et
    bos=bos$ss
    if (!is.null(exclude)) {
      oknams=vector(length=0)
      ns=length(bos); nok=ns-length(exclude)
      ok=vector(mode="list",length=nok); oki=1 # participant names with exclude
      for (i in 1:ns) if (!any(i==exclude)) {
        ok[[oki]]=bos[[i]]; oki=oki+1; oknams=c(oknams,i)
      }
      names(ok)=oknams
      bos=ok
    }
    eraw=eps(bos); traw=tps(bos); bfs=t(((eraw$ns+1)/(eraw$Ns+2))/eraw$pp)
    list(e=list(prior=eraw$pp,count=list(n=eraw$ns,N=eraw$Ns),
                bfs=bfs,pmps=bfs/(1+bfs)),
         t=list(prior=list(mono=traw$pp,overlap=traw$pps),count=list(n=traw$n,
                  N=traw$N,noverlap=traw$ns,nmono=traw$nmono)),
         test=test4(bos))
  } # end of et


  if (!is.character(bosname) || nchar(bosname)==0)
    stop("Name of sta object must be be given as quoted text")
  if (!exists(bosname,envir=.GlobalEnv))
    stop(paste("Object",bosname,"not present in global environment"))
  bos=get(bosname,envir=.GlobalEnv)
  if (class(bos)!="sta")
    stop(paste("Object",bosname,"does not have class \"sta\""))
  options(warn=-1)
  par(mfrow=c(2,2))
  if(any(is.null(exclude))) exclude="" else
    exclude=as.numeric(unlist(strsplit(as.character(exclude)," ")))
  etbos=et(bos,exclude)

  if(is.na(symb)) symb=16
  if(symb=="Unfilled circles")   {symb=1}
  if(symb=="Filled circles")     {symb=16}
  if(symb=="Unfilled triangles") {symb=2}
  if(symb=="Filled triangles")   {symb=17}
  if(symb=="Unfilled squares")   {symb=22}
  if(symb=="Filled squares")     {symb=15}
  if(symb=="Lower case letters") {
    symb=letters[as.numeric(names(etbos$t$prior$mono))]}

  if(symb=="Upper Case LETTERS") {
    symb=LETTERS[as.numeric(names(etbos$t$prior$mono))]}

  if(symb=="Participant numbers") {symb=names(etbos$t$prior$mono)}

  # Non-Trace plot using letters, LETTERS or participants numbers as symbols
  tps=etbos$test$pmps[,"NT"]
  par(mar=c(2,4,3,2)+.1)
  if(!is.numeric(symb)) {
    names(tps)=symb; tps=sort(tps)
    plot.new(); plot.window(xlim=c(1,length(tps)),ylim=c(ymin,ymax))
    text(x=(1:length(tps)),y=as.numeric(tps),labels=names(tps))
    box()
  } else {
    # Non-Trace plot using pch as symbols
    tps=sort(tps)
    plot(1:length(tps),tps,xaxt="n",yaxt="n",
         ylim=c(ymin,ymax),pch=symb,ylab="")
  }
  axis(1,labels=FALSE,tick=FALSE); axis(2,at=c(.05,.25,.5,.75,.95)); box()
  title(xlab=xlab,line=1)
  title(ylab=ylabt,line=3)
  # draw horizontal lines
  if (weakl) {
    abline(h=.05,lty=3,col=8); abline(h=.25,lty=3,col=8)
    abline(h=.5,lty=3,col=8);  abline(h=.75,lty=3,col=8)
    abline(h=.95,lty=3,col=8)
  }
  if (strongl) {
    abline(h=.01,lty=3,col=8); abline(h=.99,lty=3,col=8)
  }
  # legend for trace plot
  pall=etbos$test$pmp["NT"]
  if(plines) abline(h=pall,lty=3,lwd=2,col=1)
  if(pnames)  title(main=paste(maint,"\n","gp(Non-Trace) = ",
                    round(pall,2),sep="")) else
              title(main=maint)
  # No-Overlap comparison plot using letters,
  # LETTERS or participants numbers as symbols
  ops=etbos$test$pmps[,"NoLap"]
  if(!is.numeric(symb)) {
    names(ops)=symb; ops=sort(ops)
    plot.new(); plot.window(xlim=c(1,length(ops)),ylim=c(ymin,ymax))
    text(x=(1:length(ops)),y=as.numeric(ops),labels=names(ops))
    box()
  } else {
    # No-Overlap plot using pch as symbols
    ops=sort(ops)
    plot(1:length(ops),ops,xlab="",ylab="",yaxt="n",xaxt="n",
         ylim=c(ymin,ymax),pch=symb)
  }
  axis(1,labels=FALSE,tick=FALSE); axis(2,at=c(.05,.25,.5,.75,.95))
  title(xlab=xlab,line=1)
  title(ylab=ylabo,line=3)
  # draw horizontal lines
  if (weakl) {
    abline(h=.05,lty=3,col=8); abline(h=.25,lty=3,col=8)
    abline(h=.5,lty=3,col=8); abline(h=.75,lty=3,col=8)
    abline(h=.95,lty=3,col=8)
  }
  if (strongl) {
    abline(h=.01,lty=3,col=8); abline(h=.99,lty=3,col=8)
  }
  # legend for trace plot
  pall=etbos$test$pmp["NoLap"]
  if(plines) abline(h=pall,lty=3,lwd=2,col=1)
  if(pnames)  title(main=paste(maino,"\n","gp(No-Overlap) = ",
                    round(pall,2),sep="")) else
              title(main=maino)
  # Uni-dimensional plot using letters, LETTERS or participants numbers
  ud=etbos$test$pmps[,"UD"]
  if(!is.numeric(symb)) {
    names(ud)=symb; ud=sort(ud)
    plot.new(); plot.window(xlim=c(1,length(ud)),ylim=c(ymin,ymax))
    text(x=(1:length(ud)),y=as.numeric(ud),labels=names(ud))
    box()
  } else {
    # Uni-dimensional plot using pch as symbols
    ud=sort(ud)
    plot(1:length(ud),ud,xlab="",ylab="",yaxt="n",xaxt="n",
          ylim=c(ymin,ymax),pch=symb)
    axis(2,c(.05,.25,.5,.75,.95))
  }
  axis(1,labels=FALSE,tick=FALSE); axis(2,at=c(.05,.25,.5,.75,.95))
  title(xlab=xlab,line=1)
  title(ylab=ylabu,line=3)
  # draw horizontal lines
  if (weakl) {
    abline(h=.05,lty=3,col=8); abline(h=.25,lty=3,col=8)
    abline(h=.5,lty=3,col=8); abline(h=.75,lty=3,col=8)
    abline(h=.95,lty=3,col=8)
  }
  if (strongl) {
    abline(h=.01,lty=3,col=8); abline(h=.99,lty=3,col=8)
  }
  # legend for trace plot
  pall=etbos$test$pmp["UD"]
  if(plines) abline(h=pall,lty=3,lwd=2,col=1)
  if(pnames)  title(main=paste(mainu,"\n","gp(Uni-dimensional) = ",
                    round(pall,2),sep="")) else
              title(main=mainu)
  # Multi-dimensional plot using letters, LETTERS or participants numbers
  md=etbos$test$pmps[,"MD"]
  if(!is.numeric(symb)) {
    names(md)=symb; md=sort(md)
    plot.new(); plot.window(xlim=c(1,length(md)),ylim=c(ymin,ymax))
    text(x=(1:length(md)),y=as.numeric(md),labels=names(md))
    box()
  } else {
    # Multi-dimensional plot using pch as symbols
    md=sort(md)
    plot(1:length(md),md,xlab="",ylab="",yaxt="n",xaxt="n",
          ylim=c(ymin,ymax),pch=symb)
  }
  axis(1,labels=FALSE,tick=FALSE); axis(2,at=c(.05,.25,.5,.75,.95))
  title(xlab=xlab,line=1)
  title(ylab=ylabm,line=3)
  # draw horizontal lines
  if (weakl) {
    abline(h=.05,lty=3,col=8); abline(h=.25,lty=3,col=8)
    abline(h=.5,lty=3,col=8); abline(h=.75,lty=3,col=8)
    abline(h=.95,lty=3,col=8)
  }
  if (strongl) {
    abline(h=.01,lty=3,col=8); abline(h=.99,lty=3,col=8)
  }
  # legend for trace plot
  pall=etbos$test$pmp["MD"]
  if(plines) abline(h=pall,lty=3,lwd=2,col=1)
  if(pnames)  title(main=paste(mainm,"\n","gp(Multi-dimensional) = ",
                    round(pall,2),sep="")) else
              title(main=mainm)
}

################################################################################
### state-trace plot across all participants or for each participant ###########
################################################################################


stBootav=function(staname="",exclude=NULL,eav=TRUE,tav=TRUE,mav=TRUE,
                  nbsamp=1e4,acc=TRUE,bootparticipants=FALSE) {
  # Bootstrap over participants and posterior accuracy measures to get
  # distribution of average over participants, output suitable for plotsamps
  # stored in $all$e (encompassing) and $all$t (trace).
  # acc=T => "p" scale, F => "z" scale

  accsamp=function(bo,acc="p",stype="e") {

    accfn=function(hr,far=NA,acctype="p") {
      if (!all(is.na(far))) {
        if (acctype=="p") hr-far
        else qnorm(hr)-qnorm(far)
      } else {
        if (acctype=="p") hr
        else qnorm(hr)
      }
    }

    dsn=bo$d$dsn
    if (stype=="e") {
      if ( is.null(bo$p$e) )
        stop("Encompasing samples not present, run stSample with nkeepe>0")
      slst=bo$p$e$s
    } else if ( stype=="m" ) {
      if ( is.null(bo$p$m) )
        stop("Monotonic samples not present, run stSample with nkeept>0")
      best=names(bo$m$m$control$nm)[which.max(bo$m$m$control$nm)]
      tmp=bo$p$m$s[,,dimnames(bo$p$m$s)[[3]]==best]
      adim=dim(tmp)
      if ( length(adim)==2) {
      	adim=c(adim,1)
      	if (dsn=="B2") slst=list(tmp[,1,drop=FALSE],tmp[,2,drop=FALSE]) else {
        	tmp=array(tmp,dim=c(adim[1]/2,adim[2]*2,adim[3]))
        	slst=list(tmp[,1,drop=FALSE],tmp[,2,drop=FALSE],
                  	  tmp[,3,drop=FALSE],tmp[,4,drop=FALSE])
        }
      } else {
      	if (dsn=="B2") slst=list(tmp[,1,],tmp[,2,]) else {
        	tmp=array(tmp,dim=c(adim[1]/2,adim[2]*2,adim[3]))
        	slst=list(tmp[,1,],tmp[,2,],
                  	  tmp[,3,],tmp[,4,])
        }
      }
    } else {
      if ( is.null(bo$p$t) )
        stop("Trace samples not present, run stSample with nkeept>0")
      if (dsn=="B2")
        slst=list(bo$p$t$s[,1,],bo$p$t$s[,2,]) else
        slst=list(bo$p$t$s[,1,],bo$p$t$s[,2,],bo$p$t$s[,3,],bo$p$t$s[,4,])
    }
    nsamp=dim(slst[[1]])[2]
    if (dsn=="B2") {
      nt=(dim(slst[[1]])[1]-1)/2
      samp=array(dim=c(nt,2,2,nsamp))
      # i=chains=state, j=dimension, k=trace
      for (i in 1:2) for (j in 1:2) for (k in 1:nt)
        samp[k,j,i,]=accfn(slst[[i]][1+(j-1)*nt+k,],slst[[i]][1,],acc)
    } else { #chains = d1s1 d2s2 d1s2 s2s2
      nt=dim(slst[[1]])[1]
      samp=array(dim=c(nt,2,2,nsamp))
      for (i in 1:2) for (j in 1:2) for (k in 1:nt)
          samp[k,j,i,]=accfn(slst[[(i-1)*2+j]][k,],acc=acc)
    }
    samp
  }

  bootav=function(bos,nss,stype="e",acc="p",nbsamp=1e4,bootparticipants=FALSE) {
    ns=length(nss)                                    # number of participants
    ssamps=lapply(bos$ss,accsamp,acc=acc,stype=stype) # posterior samples
    n=unlist(lapply(ssamps,function(x){dim(x)[4]}))
    stypenam=switch(stype,e="encompassing",t="trace",m="monotonic")
    if (any(n==0))
      stop(paste("No",stypenam,"samples for following participants:\n",
                 nss[n==0]))
    if (any(n<100)) {
      cat(paste("\nLess than 100",stypenam,
                "samples available for participants:\n"))
      print(n[n<100]); cat("\n")
    }
    ddim=dim(ssamps[[1]])
    dat=array(dim=c(ddim[-4],ns))                     # one resample
    out=array(dim=c(ddim[-4],nbsamp))                 # output
    tij=matrix(runif(ns*nbsamp),nrow=nbsamp)
    if ( bootparticipants ) # resample sets of participants
       ssamp=matrix(sample(nss,nbsamp*ns,replace=TRUE),nrow=nbsamp) else
       ssamp=matrix(rep(nss,each=nbsamp),nrow=nbsamp)
    for (i in 1:nbsamp) {
      for (j in 1:ns)
        dat[,,,j] = ssamps[[ssamp[i,j]]][,,,
          min(n[ssamp[i,j]],floor(tij[i,j]*n[ssamp[i,j]])+1)]
      out[,,,i]=apply(dat,1:3,mean)
    }
    attr(out,"nss")=nss
    out
  }

  # main body of stBootav
  acc = ifelse(acc,"p","z")
  if (!is.character(staname) || nchar(staname)==0)
    stop("Name of sta object must be be given as quoted text")
  if (!exists(staname,envir=.GlobalEnv))
    stop(paste("Object",staname,"does not exist in global environment"))
  bos=get(staname,envir=.GlobalEnv)
  if (class(bos)!="sta")
    stop(paste("sta object",staname,"does not have class \"sta\""))
  nss=1:length(bos$ss)
  if( !any(is.null(exclude)) ) {
    exclude=unique(as.numeric(unlist(strsplit(as.character(exclude)," "))))
    if (!all(exclude %in% nss))
      stop("Non-existant participant excluded")
    nss=nss[-exclude]
  }
 if ( length(nss)==1 )
    stop("There is no point boostrap averaging over one participant!")
  cat("Calculating bootstrap average over participants:\n")
  print(nss)
  if ( eav ) {
    cat("\nCalculating encompassing bootstrap average\n")
    if (acc=="p")
      bos$all$p$e=bootav(bos,nss,stype="e",acc=acc,nbsamp=nbsamp,
         bootparticipants=,bootparticipants) else
      bos$all$z$e=bootav(bos,nss,stype="e",acc=acc,nbsamp=nbsamp,
         bootparticipants=,bootparticipants)
  }
  if ( tav ) {
    cat("\nCalculating trace bootstrap average\n")
    if (acc=="p")
      bos$all$p$t=bootav(bos,nss,stype="t",acc=acc,nbsamp=nbsamp,
         bootparticipants=,bootparticipants) else
      bos$all$z$t=bootav(bos,nss,stype="t",acc=acc,nbsamp=nbsamp,
         bootparticipants=,bootparticipants)
  }
  if ( mav ) {
    cat("\nCalculating best monotonic bootstrap average\n")
    if (acc=="p")
      bos$all$p$m=bootav(bos,nss,stype="m",acc=acc,nbsamp=nbsamp,
         bootparticipants=,bootparticipants) else
      bos$all$z$m=bootav(bos,nss,stype="m",acc=acc,nbsamp=nbsamp,
         bootparticipants=,bootparticipants)
  }
  assign(staname,bos,envir=.GlobalEnv)
  paste("Object",staname,"has been updated.")
}


# alls: T=bootstrap average,
#       FALSE=all participants, integer vector=corresponding participants
# line: "n"=no lines, "e"=joins data (i.e., encompassing), "t"=trace lines
#       "m"=most common monotonic order
# stat: plot either "mode", "median", "mean"
# lpt: put points on t and m lines


stPlot=function(bosname="",
                exclude=NULL,main=NA,xlab=NA,ylab=NA,
                legnamd1=NA,legnamd2=NA,acc=TRUE,alls=TRUE,symbd1=1,symbd2=2,
                stat="mode",line="e",lpts=TRUE,preg=TRUE,p=.68,
                smoothfac=5,linew=4,xmin=NA,xmax=NA,ymin=NA,ymax=NA,
                guiarg=NULL       # GUI related
                ) {

  accsamp=function(bo,acc="p",stype="e") {

    accfn=function(hr,far=NA,acctype="p") {
      if (!all(is.na(far))) {
        if (acctype=="p") hr-far
        else qnorm(hr)-qnorm(far)
      } else {
        if (acctype=="p") hr
        else qnorm(hr)
      }
    }

    dsn=bo$d$dsn
    if (stype=="e") {
      if ( is.null(bo$p$e) )
        stop("Encompasing samples not present, run stSample with nkeepe>0")
      slst=bo$p$e$s
    } else if ( stype=="m" ) {
      if ( is.null(bo$p$m) )
        stop("Monotonic samples not present, run stSample with nkeept>0")
      best=names(bo$m$m$control$nm)[which.max(bo$m$m$control$nm)]
      tmp=bo$p$m$s[,,dimnames(bo$p$m$s)[[3]]==best]
      adim=dim(tmp)
      if (dsn=="B2") slst=list(tmp[,1,],tmp[,2,]) else {
        tmp=array(tmp,dim=c(adim[1]/2,adim[2]*2,adim[3]))
        slst=list(tmp[,1,],tmp[,2,],tmp[,3,],tmp[,4,])
      }
    } else {
      if ( is.null(bo$p$t) )
        stop("Trace samples not present, run stSample with nkeept>0")
      if (dsn=="B2")
        slst=list(bo$p$t$s[,1,],bo$p$t$s[,2,]) else
        slst=list(bo$p$t$s[,1,],bo$p$t$s[,2,],bo$p$t$s[,3,],bo$p$t$s[,4,])
    }
    nsamp=dim(slst[[1]])[2]
    if (dsn=="B2") {
      nt=(dim(slst[[1]])[1]-1)/2
      samp=array(dim=c(nt,2,2,nsamp))
      # i=chains=state, j=dimension, k=trace
      for (i in 1:2) for (j in 1:2) for (k in 1:nt)
        samp[k,j,i,]=accfn(slst[[i]][1+(j-1)*nt+k,],slst[[i]][1,],acc)
    } else { #chains = d1s1 d2s2 d1s2 s2s2
      nt=dim(slst[[1]])[1]
      samp=array(dim=c(nt,2,2,nsamp))
      for (i in 1:2) for (j in 1:2) for (k in 1:nt)
          samp[k,j,i,]=accfn(slst[[(i-1)*2+j]][k,],acc=acc)
    }
    samp
  }

  mode2d=function(samp,smoothfac=5) {

    getmode2d=function(dns) {
      n=length(dns$x1)
      x=which.max(dns$fhat)
      div=x%/%n; mod=x%%n
      if(mod==0) c(dns$x1[n],dns$x2[div]) else
                 c(dns$x1[mod],dns$x2[div+1])
    }

    require(KernSmooth)
    nt=dim(samp)[1]; nd=dim(samp)[2]
    smth=0
    for (i in 1:nt) for (j in 1:nd) for (k in 1:2) {
      tmp=dpik(samp[i,j,k,])
      if (tmp>smth) smth=tmp
    }
    out=array(dim=c(nt,nd,2))
    for (i in 1:nt) for (j in 1:nd) {
      dns=bkde2D(t(samp[i,j,,]),smoothfac*smth)
      out[i,j,]=getmode2d(dns)
    }
    out
  }

  getline=function(lsamp,stat="mode",smoothfac=5) {
    if (stat=="mean")  lins=apply(lsamp,1:3,mean)   else
    if (stat=="median") lins=apply(lsamp,1:3,median) else {
                       lins=mode2d(lsamp,smoothfac)
    }
    lins
  }

  plotsamp=function(samp,pts,lins,line,linew=2,main="",preg=TRUE,printacc=TRUE,
                    p=.5,xmin=NA,xmax=NA,ymin=NA,ymax=NA,acc="p") {

    findp=function(p,x) {
      obj=function(t,p,x) abs((sum(x[x>t])/sum(x))-p)
      optimise(f=obj,interval=range(x),p=p,x=x)$minimum
    }

    getlims=function(samp,preg,mns,xmin,xmax,ymin,ymax,p,sfsmth) {

      nt=dim(samp)[1]; nd=dim(samp)[2]
      contourarray=array(dim=c(nt,nd))
      mode(contourarray)="list"
      xlim=c(Inf,-Inf); ylim=c(Inf,-Inf)
      if (preg) {
        for (i in 1:nt) for (j in 1:nd) {
          dns=bkde2D(t(samp[i,j,,]),sfsmth)
          xy=contourLines(x=dns$x1,y=dns$x2,z=dns$fhat,
                          levels=findp(p=p,x=dns$fhat))
          xlim[1]=min(c(xlim[1],xy[[1]]$x))
          xlim[2]=max(c(xlim[2],xy[[1]]$x))
          ylim[1]=min(c(ylim[1],xy[[1]]$y))
          ylim[2]=max(c(ylim[2],xy[[1]]$y))
          contourarray[i,j][[1]]=xy[[1]]
        }
      } else {
        xlim[1]=min(mns[,,1])
        xlim[2]=max(mns[,,1])
        ylim[1]=min(mns[,,2])
        ylim[2]=max(mns[,,2])
      }
      xlim[1]=floor(xlim[1]*20)/20 - .05
      xlim[2]=ceiling(xlim[2]*20)/20 + .05
      ylim[1]=floor(ylim[1]*20)/20 - .05
      ylim[2]=ceiling(ylim[2]*20)/20 + .05
      list(xmin=ifelse(is.na(xmin),xlim[1],xmin),
           xmax=ifelse(is.na(xmax),xlim[2],xmax),
           ymin=ifelse(is.na(ymin),ylim[1],ymin),
           ymax=ifelse(is.na(ymax),ylim[2],ymax),
           contourarray=contourarray)
    }

    # Main body of plotsamp
    nt=dim(samp)[1]; nd=dim(samp)[2]
    if ( preg ) {
      require(KernSmooth)
      smth=0
      for (i in 1:nt) for (j in 1:nd) for (k in 1:2) {
        tmp=dpik(samp[i,j,k,])
        if (tmp>smth) smth=tmp
      }
    }
    lims=getlims(samp,preg,mns=pts,xmin,xmax,ymin,ymax,p,sfsmth=smoothfac*smth)
    plot( as.vector(pts[,1,1]),as.vector(pts[,1,2]),
          xlim=c(lims$xmin,lims$xmax),ylim=c(lims$ymin,lims$ymax),
          xlab=xlab,ylab=ylab,pch=symbd1,cex=3,main=main)
    points(as.vector(pts[,2,1]),as.vector(pts[,2,2]),pch=symbd2,cex=3)
    points(as.vector(pts[,1,1]),as.vector(pts[,1,2]),pch=paste(1:nt),cex=1)
    points(as.vector(pts[,2,1]),as.vector(pts[,2,2]),pch=paste(1:nt),cex=1)

    dimnames(pts)=list(Trace=1:nt,
                       Dimension=c(c(legnamd1,legnamd2)),
                       State=c(xlab,ylab))
    cat(paste("\n",main,"\n\n"))
    print(pts)
    print(c(Average=mean(pts)))

    if (!is.null(lins)) {
      if (line=="e" || line=="t") {
        lines(as.vector(lins[,1,1]),as.vector(lins[,1,2]),lwd=linew)
        lines(as.vector(lins[,2,1]),as.vector(lins[,2,2]),lty=2,lwd=linew)
        if (line=="t" && lpts) {
          points(as.vector(lins[,1,1]),as.vector(lins[,1,2]),pch=symbd1,cex=1.5)
          points(as.vector(lins[,2,1]),as.vector(lins[,2,2]),pch=symbd2,cex=1.5)
        }
      } else if (line=="m") {
        lines(sort(as.vector(lins[,,1])),sort(as.vector(lins[,,2])),lwd=linew)
        if (lpts) {
          pchs=rep(c(symbd1,symbd2),each=dim(lins)[1])[order(lins[,,1])]
          points(sort(as.vector(lins[,,1])),sort(as.vector(lins[,,2])),
                 pch=pchs,cex=1.5)
        }
      }
    }
    if ( !is.null(c(legnamd1,legnamd2)) ) { # draw legend
      if (length(c(legnamd1,legnamd2))!=2)
        stop("Supply two names for data traces (legnamd1 and legnamd2)")
      legend("bottomright",legend=c(legnamd1,legnamd2),pch=c(symbd1,symbd2),
             cex=1.5,bty="n",title="Data")
      if(line!="n" && !is.null(lins) ) {
        tit=switch(line,
                    e="Data Traces",
                    t="Trace Model",
                    m="Best Monotonic Model"
                  )
        if (line=="m")
          legend("topleft",legend=tit,lty=c(1),lwd=2,cex=1.5,bty="n") else
          legend("topleft",legend=c(legnamd1,legnamd2),lty=c(1,2),lwd=2,
                 cex=1.5,bty="n",title=tit)
      }
    }
    if (preg) {
      for (i in 1:nt) for (j in 1:nd)
        lines(lims$contourarray[i,j][[1]])
    }
  } # End of plotsamp

  # main body of stPlot
  acc = ifelse(acc,"p","z")
  if (!is.character(bosname) || nchar(bosname)==0)
    stop("Name of sta object must be be given as quoted text")
  if (!exists(bosname,envir=.GlobalEnv))
    stop(paste("Object",bosname,"not present in global environment"))
  bos=get(bosname,envir=.GlobalEnv)
  if (class(bos)!="sta")
    stop(paste("Object",bosname,"does not have class \"sta\""))
  dsn=bos$ss[[1]]$d$dsn
  if (is.na(xlab)) xlab=attr(bos$ss,"slevs")[1]
  if (is.na(ylab)) ylab=attr(bos$ss,"slevs")[2]
  if (acc=="p") {
    if (dsn=="B2") {
      xlab=paste(xlab," (HR-FAR)")
      ylab=paste(ylab," (HR-FAR)")
    } else {
      xlab=paste(xlab," (HR)")
      ylab=paste(ylab," (HR)")
    }
  } else {
    if (dsn=="B2") {
      xlab=paste(xlab," (zHR-zFAR)")
      ylab=paste(ylab," (zHR-zFAR)")
    } else {
      xlab=paste(xlab," (zHR)")
      ylab=paste(ylab," (zHR)")
    }
  }

  if (is.na(legnamd1)) legnamd1=attr(bos$ss,"dlevs")[1]
  if (is.na(legnamd2)) legnamd2=attr(bos$ss,"dlevs")[2]
  if (is.na(symbd1)) symbd1=2
  if (is.na(symbd2)) symbd2=6
  if (is.na(line)) line="n"
  if (is.na(stat)) stat="mode"

  if ( is.character(stat) ) {
    if(stat=="Posterior mode") {stat="mode"}
    if(stat=="Posterior mean") {stat="mean"}
    if(stat=="Posterior median") {stat="median"}
  }
  if ( is.character(symbd1) ) {
    if(symbd1=="Unfilled circles") {symbd1=1}
    if(symbd1=="Unfilled upright triangles") {symbd1=2}
    if(symbd1=="Unfilled inverted triangles") {symbd1=6}
    if(symbd1=="Unfilled squares") {symbd1=22}
    if(symbd1=="Unfilled diamonds") {symbd1=5}
  }
  if ( is.character(symbd2) ) {
    if(symbd2=="Unfilled circles") {symbd2=1}
    if(symbd2=="Unfilled upright triangles") {symbd2=2}
    if(symbd2=="Unfilled inverted triangles") {symbd2=6}
    if(symbd2=="Unfilled squares") {symbd2=22}
    if(symbd2=="Unfilled diamonds") {symbd2=5}
  }
  if ( is.character(line) && nchar(line)>1 ) {
    if(line=="No lines") {line="n"}
    if(line=="Data Traces") {line="e"}
    if(line=="Trace Model") {line="t"}
    if(line=="Monotonic Model") {line="m"}
  }
  if(is.na(line) | is.na(symbd1) | is.na(symbd2))
    stop("Select line type and/or plotting symbols")

  nss=1:length(bos$ss)
  if( !any(is.null(exclude)) ) {
    exclude=unique(as.numeric(unlist(strsplit(as.character(exclude)," "))))
    if (!all(exclude %in% nss))
      stop("Non-existant participant excluded")
    nss=nss[-exclude]
  }
  if ( is.logical(alls) ) {
    if (!alls) {
      for (i in nss) {
        if (is.na(main) | length(main)!=length(nss))
           main=paste("Participant ",i," posterior ",stat,
                      "\np(posterior region)=",p,sep="")
        bo=bos$ss[[i]]
        samp=accsamp(bo,acc=acc)
        pts=getline(samp,stat=stat,smoothfac=smoothfac)
        if (line=="n") lins=NULL else
        if (line=="e")  lins=pts  else {
          lsamp=accsamp(bo,acc=acc,stype=line)
          nlsamp=dim(lsamp)[4]
          if (nlsamp<100)
            cat(paste("Line for participant",i,"is based on only",
                      dim(lsamp)[4],"samples\n"))
          if (is.na(nlsamp) || nlsamp<5) {
            cat("Cannot draw line for this participant!\n\n")
            lins=NULL
          } else lins=getline(lsamp,stat=stat,smoothfac=smoothfac)
        }
       .Options$device()
        plotsamp(samp,pts=pts,lins=lins,line=line,linew=linew,main=main,acc=acc,
                 preg=preg,p=p,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
      }
    } else {
      if ( (acc=="p" && is.null(bos$all$p$e)) ||
           (acc!="p" && is.null(bos$all$z$e)) )
          stop(paste("Encompasing boostrap average for this accuracy measure",
                     "not present, run stBootav"))
      if ( line=="t" && ((acc=="p" && is.null(bos$all$p$t)) ||
                         (acc!="p" && is.null(bos$all$z$t))) )
          stop(paste("Trace boostrap average for this accuracy measure",
                     "not present, run stBootav"))
      if ( line=="m" && ((acc=="p" && is.null(bos$all$p$m)) ||
                         (acc!="p" && is.null(bos$all$z$m))) )
          stop(paste("Monotonic boostrap average for this accuracy measure",
                     "not present, run stBootav"))
      if ( acc=="p" ) storednss=attr(bos$all$p$e,"nss") else
                        storednss=attr(bos$all$z$e,"nss")
      if (!all(storednss==nss) )
        stop(paste("Included participants do not match encompassing bootstrap",
                   "average, rerun stBootav"))
      if ( line=="t" ) {
        if ( acc=="p" ) storednss=attr(bos$all$p$t,"nss") else
                          storednss=attr(bos$all$z$t,"nss")
          if ( !all(storednss==nss) )
            stop(paste("Included participants do not match trace bootstrap",
                       "average, rerun stBootav"))
      }
      if (line=="m") {
        if ( acc=="p" ) storednss=attr(bos$all$p$m,"nss") else
                          storednss=attr(bos$all$z$m,"nss")
        if ( !all(storednss==nss) )
          stop(paste("Included participants do not match monotonic bootstrap",
                     "average, rerun stBootav"))
      }
      if (is.na(main))
         main=paste("Participant average posterior ",stat,
            "s\np(posterior region)=",p,sep="")
      pts=getline(switch(acc,p=bos$all$p$e,bos$all$z$e),
                  stat=stat,smoothfac=smoothfac)
      if (line=="n") lins=NULL else
      if (line=="e") lins=pts  else
        lins=getline(switch(acc,p=switch(line,t=bos$all$p$t,m=bos$all$p$m),
                                switch(line,t=bos$all$z$t,m=bos$all$z$m)),
                     stat=stat,smoothfac=smoothfac)
      plotsamp(switch(acc,p=bos$all$p$e,bos$all$z$e),
               pts=pts,lins=lins,line=line,linew=linew,main=main,
               preg=preg,p=p,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,acc=acc)
    }
  } else {
    for (i in alls) {
        if (all(is.na(main)) | length(main) != length(alls))
           main=paste("Participant ",i," posterior ",stat,
                      "s\np(posterior region)=",p,sep="")
        bo=bos$ss[[i]]
        samp=accsamp(bo,acc=acc)
        pts=getline(samp,stat=stat,smoothfac=smoothfac)
        if (line=="n") lins=NULL else
        if (line=="e") lins=pts  else {
          lsamp=accsamp(bo,acc=acc,stype=line)
          nlsamp=dim(lsamp)[4]
          if (nlsamp<100)
            cat(paste("Line for participant",i,"is based on only",
                      dim(lsamp)[4],"samples\n"))
          if (is.na(nlsamp) || nlsamp<5) {
            cat("Cannot draw line for this participant!\n\n")
            lins=NULL
          } else lins=getline(lsamp,stat=stat,smoothfac=smoothfac)
        }
        if (length(alls)>1) .Options$device()
        plotsamp(samp,pts=pts,lins=lins,line=line,linew=linew,main=main,acc=acc,
                 preg=preg,p=p,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
    }
  }
}


stFirst=function(staname="",fnams=NULL,folder="",extension="txt",
                 sep="tab",multiparticipant=FALSE,header=TRUE,
                 usecols=NULL,na.strings="NA",acc=TRUE) {
  staMake(staname=staname,fnams=fnams,folder=folder,extension=extension,
          usecols=usecols,header=header,sep=sep,na.strings=na.strings,
          multiparticipant=multiparticipant)
  stSample(staname,verbose=2)
  staManage(staname,checkConverge=TRUE)
  stBootav(staname,eav=TRUE,tav=FALSE,mav=FALSE,acc=acc)
  stPlot(staname,alls=FALSE,acc=acc)
  .Options$device()
  stPlot(staname,acc=acc)
  cat(paste("\nFirst pass completed!\n",
    "Run stSample to obtain accurate posterior model probability estimates\n",
    "then examine results with stSummary and stProbplot.\n",
    "Run stBootav to obtain participant results for average trace and\n",
    "monotonic plots with stPlot\n",
    "Use staManage to join sta objects and manage stored samples\n"))
}



################################################################################
################################ GUIS  #########################################
################################################################################

guistFirst <- function() {
  require(fgui)
  guiSet("ENTRY_WIDTH",20)
  guiSet("LIST_WIDTH",15)
  guiSet("LIST_HEIGHT",5)
  guiv(stFirst,
        title="Generate sta object and complete the first pass",
        argText=c(staname="* sta object name (character string)",
          fnams="Character string of directory + file name of data file/s",
          folder="Character string of directory containing data file/s",
          extension="Character string of file extension",
          sep="* File delimiter",
          usecols="Columns to use from each data file",
          header="Header row in data files?",
          na.strings="String for empty cells in data file",
          multiparticipant="Selected data files each contain data for multiple participants",
          acc="Accuracy based on probabilities"),
        argList=list(sep=NULL),
        argCommand=list(guiarg=staMake_press()),
        argType=list(guiarg="i"),
        argOption=list(header=c("T","F"),multiparticipant=c("T","F"),acc=c("T","F")),
        argGridSticky=c(rep("e",10)),
        closeOnExec=TRUE,
		helps=NULL
  )
}

staMake_press=function() {
  setListElements("sep",c("tab","space","comma"))
}

# gui for subsequent passes of boast
guistSample <- function() {
  require(fgui)
  guiSet("ENTRY_WIDTH",20)
  guiSet("SLIDER_LENGTH",150)
  guiv(stSample,
        title="Refine sampling from the encompassing and trace models",
        argText=c(bosname="* sta object name (character string)",
          refresh="Refresh sta object calculations",
          maxt="Maximum run time (hours)",
          ci="Credible interval (0-100%)",
          BFe="Run encompassing model",
          sampe="Samples per run for encompassing model",
          ed="Credible interval precision for encompassing samples (0-1)",
          nkeepe="Number of encompassing model samples to keep for plotting",
          BFt="Run trace model",
          sampt="Samples per run for trace model",
          burn="Number of burn-in samples",
          td="Credible interval precision for trace samples (0-1)",
          nkeept="Number of trace model samples to keep for plotting",
          nkeepm="Number of monotonic model samples to keep for plotting",
          verbose="verbose"),
        argSlider=list(verbose=c(0,2,1)),
        argOption=list(refresh=c("T","F"),BFe=c("T","F"),BFt=c("T","F")),
        argGridSticky=c(rep("e",15)),
        closeOnExec=TRUE,
		helps=NULL
  )
}


# gui for summary statistics
guistSummary=function() {
  require(fgui)
  guiSet("ENTRY_WIDTH",20)
  guiSet("LIST_WIDTH",50)
  guiSet("LIST_HEIGHT",12)
  guiv(stSummary,
        title="Extract summary model selection results",
        argText=c(bosname="* sta object name (character string)",
          BF="Report Bayes Factors (T) or probabilities (F)",
          traceTrue="Use trace-true (T) or exhaustive (F) strategy for probabilities",
          eachs="Display values for individual participants",
          nsig="Round to how many decimal places?",
          exclude="Participants to exclude",
          sorts="Sort values for individual participants by model",
          extras="Select additional results to display"
        ),
        argList=list(extras=NULL,sorts=NULL),
        argCommand=list(guiarg=stSummary_press()),
        argOption=list(BF=c("T","F"),traceTrue=c("T","F"),eachs=c("T","F")),
        argType=list(guiarg="i"),
        argGridSticky=c(rep("e",6),rep("n",3)),
        closeOnExec=TRUE,
		helps=NULL
  )
}

stSummary_press=function() {
  setListElements("sorts",c("Non-Trace model","No-Overlap model",
                  "Uni-dimensional model","Multi-dimensional model")
  )
  setListElements("extras",
                    c("Data sources",
                    "Trace model prior probability",
                    "No-Overlap model prior probability",
                    "Uni-dimensional model prior probability",
                    "Multi-dimensional model prior probability",
                    "Non-Trace model samples from encompassing model",
                    "Total samples from encompassing model",
                    "Non-overlapping monotonic samples from trace model",
                    "Overlapping monotonic samples from trace model",
                    "Non-monotonic samples from trace model",
                    "Monotonic model counts",
                    "Total samples from trace model")
  )
}

# gui for stProbplot
guistProbplot=function() {
  require(fgui)
  guiSet("SLIDER_LENGTH",200)
  guiSet("ENTRY_WIDTH",30)
  guiSet("LIST_WIDTH",20)
  guiSet("LIST_HEIGHT",10)
  guiSet("EDIT_WIDTH",45)
  guiSet("EDIT_HEIGHT",1)
  gui(stProbplot,
      title="Plot individual participant model selection probabilities",
      argText=c(bosname="* sta object name (character string)",
        exclude="Participants to exclude",
        maint="Non-Trace title",
        maino="No-Overlap title",
        mainu="Uni-dimensional title",
        mainm="Multi-dimensional title",
        xlab="x axis label",
        ylabt="Non-Trace y axis label",
        ylabo="No-Overlap y axis label",
        ylabu="Uni-dimensional model y axis label",
        ylabm="Multi-dimensional y axis label",
        ymin="Minimum y axis value",
        ymax="Maximum y axis value",
        weakl="Insert dashed lines at p(0.05,0.25,0.5,0.75,0.95)",
        strongl="Insert dashed lines at p(0.01,0.99)",
        symb="Select plotting symbols",
        plines="Insert line at group model selection probability",
        pnames="Insert group model selection probability in title"
      ),
      argList=list(symb=NULL),
      argCommand=list(guiarg=stProbplot_press()),
      argOption=list(
        weakl=c("T","F"),
        strongl=c("T","F"),
        plines=c("T","F"),
        pnames=c("T","F")
      ),
      argType=list(guiarg="i"),argGridOrder=c(1:12,13,14:19,20,21),
      argGridSticky=c(rep("e",23)),
      argSlider=list(ymin=c(0,1,0.05),ymax=c(0,1,0.05)),
	  helps=NULL
  )
}

stProbplot_press=function() {
  setListElements("symb",
                  c("Unfilled circles",
                    "Filled circles",
                    "Unfilled triangles",
                    "Filled triangles",
                    "Unfilled squares",
                    "Filled squares",
                    "Lower case letters",
                    "Upper Case LETTERS",
                    "Participant numbers")
                  )
}

# gui for stBootav
guistBootav=function() {
  require(fgui)
  guiSet("ENTRY_WIDTH",20)
  guiv(stBootav,
        title="Generate bootstrap averages for state-trace plots",
        argText=c(staname="* sta object name (character string)",
          exclude="Participants to exclude",
          eav="Generate bootstrap average for encompassing model",
          tav="Generate bootstrap average for trace model",
          mav="Generate bootstrap average for monotonic model",
          nbsamp="Number of bootstrap samples to draw",
          acc="Accuracy based on probabilities",
          bootparticipants="Resample participants"),
        argOption=list(eav=c("T","F"),tav=c("T","F"),mav=c("T","F"),acc=c("T","F"),bootparticipants=c("T","F")),
        argGridSticky=rep("e",8),
		helps=NULL
  )
}

# gui for stPlot
guistPlot=function() {
  require(fgui)
  guiSet("SLIDER_LENGTH",200)
  guiSet("ENTRY_WIDTH",30)
  guiSet("LIST_WIDTH",25)
  guiSet("LIST_HEIGHT",6)
  guiSet("EDIT_WIDTH",50)
  guiSet("EDIT_HEIGHT",1)
  gui(stPlot,
      title="Generate state-trace plots",
      argText=c(bosname="* sta object name (character string)",
        exclude="Participants to exclude",
        main="Main title",
        xlab="x axis label",
        ylab="y axis label",
        legnamd1="Dimension 1 label",
        legnamd2="Dimension 2 label",
        acc="Accuracy based on probabilities",
        alls="Plot average (T) or each participant (F)",
        symbd1="Plot symbols for dimension 1",
        symbd2="Plot symbols for dimension 2",
        stat="Measure of Central Tendency",
        line="Line type",
        lpts="Place points on plot line/s",
        preg="Plot credible p regions",
        p="Width of credible p regions",
        smoothfac="Smoothing factor for p regions",
        linew="Width of plot line/s",
        xmin="Minimum x axis value",
        xmax="Maximum x axis value",
        ymin="Minimum y axis value",
        ymax="Maximum y axis value"
      ),
      argList=list(stat=NULL,line=NULL,symbd1=NULL,symbd2=NULL),
      argCommand=list(guiarg=stPlot_press()),
      argOption=list(acc=c("T","F"),alls=c("T","F"),preg=c("T","F"),lpts=c("T","F")),
      argType=list(guiarg="i"),
      argGridSticky=c(rep("e",14),rep("w",4),rep("e",4),"n"),
      argGridOrder=c(1,2,3,4,5,6,7,8,9,10,10,11,11,12,12,13:20),
      argSlider=list(
        p=c(0,1,0.02),
        smoothfac=c(0.5,10,0.5),
        linew=c(0,10,1)
      ),
	  helps=NULL
  )
}

stPlot_press=function() {
  setListElements("stat",
    c("Posterior mode",
      "Posterior mean",
      "Posterior median")
  )
  setListElements("line",
    c("No lines",
      "Data Traces",
      "Trace Model",
      "Monotonic Model")
  )
  setListElements("symbd1",
    c("Unfilled circles",
      "Unfilled upright triangles",
      "Unfilled inverted triangles",
      "Unfilled squares",
      "Unfilled diamonds")
  )
  setListElements("symbd2",
    c("Unfilled circles",
      "Unfilled upright triangles",
      "Unfilled inverted triangles",
      "Unfilled squares",
      "Unfilled diamonds")
  )
}

# gui for staManage
guistaManage=function() {
  require(fgui)
  guiSet("ENTRY_WIDTH",20)
  guiv(staManage,
        title="Manage sta object",
        argText=c(stanames="* sta object name/s (character string)",
          staname="New name for sta object (optional, character string)",
          checkConverge="Check convergence of trace model MCMC chains?",
          nmcmc="Length of each MCMC chain",
          mcmcSave="Name of saved MCMC samples",
          nkeepe="Number of encompassing model samples to keep",
          nkeept="Number of trace model samples to keep",
          nkeepm="Number of monotonic model samples to keep",
          keepbestm="Keep only samples for the best monotonic model?"
        ),
        argOption=list(checkConverge=c("T","F"),keepbestm=c("T","F")),
        argGridSticky=rep("e",9),
		helps=NULL
  )
}

sta=function(stFirst,stSample,stSummary,stProbplot,stBootav,stPlot,staManage) {}

stacallback=function(arg) {
  if(arg=="stFirst") guistFirst()
  if(arg=="stSample") guistSample()
  if(arg=="stSummary") guistSummary()
  if(arg=="stProbplot") guistProbplot()
  if(arg=="stBootav") guistBootav()
  if(arg=="stPlot") guistPlot()
  if(arg=="staManage") guistaManage()
}

guista=function() {
  require(fgui)
  guiv(sta,title="Select option",
  argText=c(stFirst="stFirst: Generate sta object and run first pass of sampling",
            stSample="stSample: Refine sampling from the encompassing and trace models",
            stSummary="stSummary: Extract summary model selection results",
            stProbplot="stProbplot: Plot individual participant model selection probabilities",
            stBootav="stBootav: Generate bootstrap averages for state-trace plots",
            stPlot="stPlot: State-Trace plots for the group average and for individual participants",
            staManage="staManage: Manage sta object/s"),
  argCommand=list(stFirst=NULL,stSample=NULL,stSummary=NULL,
                  stProbplot=NULL,stBootav=NULL,stPlot=NULL,
                  staManage=NULL),
  callback=stacallback,
  exec=NULL,
  argGridSticky=c(rep("n",8)),
  argGridOrder=1:8,modal=FALSE,
  helps=NULL
  )
}
