getBreaks_th <- function(x,y,threshold=40,postprocessing=list(post=TRUE,adjacent=2),col="GrayLevel",plot=TRUE,
         shiny=FALSE,...){
  
  if (!is.numeric(threshold)){
    stop("threshold must be a percent strictly between 0 and 100")
  } else if ((threshold<=0)||(threshold>=100)||(length(threshold)!=1)){
    stop("threshold must be a percent strictly between 0 and 100")
  }
  if (!is.list(postprocessing)){
    stop("postprocessing must be a list")
  } else {
    if (!is.logical(postprocessing$post)){
      stop("postprocessing$post must be logical")
    } else if (postprocessing$post){
      if (!is.numeric(postprocessing$adjacent)){
        stop("postprocessing$adjacent must be a positive integer")
      } else if ((postprocessing$adjacent<=0)||(floor(postprocessing$adjacent)!=postprocessing$adjacent)){
        stop("postprocessing$adjacent must be a positive integer")
      }
    }
  }
  if (!(is.matrix(y)||(class(y)=="dgeMatrix"))){
    stop("y must be the observations data (or a transformation)")
  }
  if (!is.character(col)){
    stop("col must be a character")
  } else {
    if (length(col)==1){
      if (col=="GrayLevel"){
        couleur=gray(seq(0,1, length=256))
        couleur=couleur[length(couleur):1]
      }else{
        couleur=c(grDevices::rgb((0:201)/201*200,(0:201)/201*200,255,maxColorValue = 255),"white",
                  grDevices::rgb(255,(0:201)/201*200,(0:201)/201*200,maxColorValue = 255)[201:0])
      }
    }else{
      couleur=col
    }
  }
  if (!is.logical(shiny)){
    stop("shiny must be logical")
  }
  if (shiny){
    ############### Version Shiny
    ui_blockseg = shiny::fluidPage(
      shiny::fluidRow(
        shiny::column(3,shiny::sliderInput("Thres","Threshold (%)",min=0,max=100,step=0.1,
                                           value=floor(threshold*10)/10)),
        shiny::column(3,shiny::checkboxInput("Gray", "Gray version", col=="GrayLevel")),
        shiny::column(3,shiny::checkboxInput("post", "Post-processing?", postprocessing$post)),
        shiny::column(3,shiny::numericInput('adj', 'Number of neighbours', postprocessing$adjacent))),
      shiny::column(12,shiny::plotOutput(outputId="plot", width = "100%",
                                         height = paste0(as.character(par("din")[2]*160),"px")))
    )
    ########### Partie Server
    server_blockseg = function(input, output) {
      output$plot<- shiny::renderPlot({
        if (input$Gray){
          plot(x,y,threshold=input$Thres,postprocessing=list(post=input$post,
                                                             adjacent=input$adj),
               col="GrayLevel",shiny=FALSE)
        }else{
          if (col=="GrayLevel"){
            col2="Color"
          }else{
            col2=col
          }
          plot(x,y,threshold=input$Thres,postprocessing=list(post=input$post,
                                                             adjacent=input$adj),
               col=col2,shiny=FALSE)
        }
      })
    }
    ########### Affichage
    shiny::shinyApp(ui = ui_blockseg,server=server_blockseg)
    
  }else{
    brek=seq(min(y),max(y),length=length(couleur)+1)
    
    ValSeuilline=threshold/100*max(x@RowBreaks)
    ValSeuilCol=threshold/100*max(x@ColBreaks)
    ColLine=x@RowBreaks>ValSeuilline
    ColCol=x@ColBreaks>ValSeuilCol
    RowBreak=which(ColLine)
    ColBreak=which(ColCol)
    
    n=nrow(y)
    d=ncol(y)
    
    if (postprocessing$post){
      nr=length(RowBreak)
      if (nr>1){
        dist=RowBreak[2:nr]-RowBreak[1:(nr-1)]
        Rowbreakbis=c()
        ind=1
        Rowcompl=c()
        nr=nr-1
        while(ind<=nr){
          if (dist[ind]<postprocessing$adjacent){
            inddeb=ind
            while((dist[ind]<postprocessing$adjacent)&(ind<nr)){
              ind=ind+1
            }
            indfin=ind
            indref=inddeb-1+which.max(x@RowBreaks[RowBreak[inddeb:indfin]])
            Rowbreakbis=c(Rowbreakbis,RowBreak[indref])
            suite=inddeb:indfin
            Rowcompl=c(Rowcompl,RowBreak[suite[suite!=indref]])
          }else{
            Rowbreakbis=c(Rowbreakbis,RowBreak[ind])
          }
          ind=ind+1
        }
        if (dist[nr]>=postprocessing$adjacent){
          Rowbreakbis=c(Rowbreakbis,RowBreak[nr+1])
        }else{
          Rowcompl=c(Rowcompl,RowBreak[nr+1])
        }
        RowBreak=Rowbreakbis
        ColLine[Rowcompl]=ColLine[Rowcompl]+2
      }
      nc=length(ColBreak)
      if (nc>1){
        dist=ColBreak[2:nc]-ColBreak[1:(nc-1)]
        Colbreakbis=c()
        ind=1
        Colcompl=c()
        nc=nc-1
        while(ind<=nc){
          if (dist[ind]<postprocessing$adjacent){
            inddeb=ind
            while((dist[ind]<postprocessing$adjacent)&(ind<nc)){
              ind=ind+1
            }
            indfin=ind
            indref=inddeb-1+which.max(x@ColBreaks[ColBreak[inddeb:indfin]])
            Colbreakbis=c(Colbreakbis,ColBreak[indref])
            suite=inddeb:indfin
            Colcompl=c(Colcompl,ColBreak[suite[suite!=indref]])
          }else{
            Colbreakbis=c(Colbreakbis,ColBreak[ind])
          }
          ind=ind+1
        }
        if (dist[nc]>=postprocessing$adjacent){
          Colbreakbis=c(Colbreakbis,ColBreak[nc+1])
        }else{
          Colcompl=c(Colcompl,ColBreak[nc+1])
        }
        ColBreak=Colbreakbis
        ColCol[Colcompl]=ColCol[Colcompl]+2
      }
    }
    
    ## matrice resumee
    emplz=diff(unique(c(1,RowBreak,n+1)))
    z=bdiag(lapply(emplz,function(i) rep(1,i) ))
    nz=length(emplz)
    emplw=diff(unique(c(1,ColBreak,d+1)))
    w=bdiag(lapply(emplw,function(i) rep(1,i)))
    nw=length(emplw)
    mean =as.matrix(t(z)%*%y%*%w/(emplz%*%t(emplw)))
    var = pmax(0,as.matrix(t(z)%*%(y^2)%*%w/(emplz%*%t(emplw)))-mean^2)
    resum = mean
    
    if(plot){
      par(mfrow=c(2,2),oma=c(0,0,3,0))
      
      ### Affichage matrice originale
      
      image(1:d,1:n,t(y)[,n:1],xlab="",ylab="",xaxt="n",yaxt="n",main="Original data",col=couleur,breaks=brek)
      abline(v=ColBreak-0.5,col="purple")
      abline(h=n-RowBreak+1.5,col="purple")
      
      ## Affichage plot ligne
      
      plot(x@RowBreaks,n:1,col=ColLine+1,xlab="",ylab="n", ylim=c(1,n),yaxt="n")
      if (length(RowBreak)<=1){
        title(main=paste(as.character(length(RowBreak))," row break",sep=""))
      }else{
        title(main=paste(as.character(length(RowBreak))," row breaks",sep=""))
      }
      abline(v=ValSeuilline,col="purple")
      nax=floor(6*par("din")[2]/4)
      axis(2,at=floor(seq(0,n,length=nax)),labels=floor(seq(n,0,length=nax)))
      
      
      ## Affichage plot colonne
      
      plot(1:d,x@ColBreaks,col=ColCol+1,xlab="n", ylab="")
      if (length(ColBreak)<=1){
        title(main=paste(as.character(length(ColBreak))," column break",sep=""))
      }else{
        title(main=paste(as.character(length(ColBreak))," column breaks",sep=""))
      }
      abline(h=ValSeuilCol,col="purple")
      
      
      if (is.matrix(resum)){
        if (nz!=1){
          image(x=unique(c(1,ColBreak,d)),y=unique(n-c(n,RowBreak[length(RowBreak):1],1)+1),
                z=t(resum[nz:1,]),
                xlab="",ylab="",xaxt="n",yaxt="n",main="Summarized data",col=couleur,breaks=brek)
        }else{
          image(x=unique(c(1,ColBreak,d)),y=unique(n-c(n,RowBreak[RowBreak:1],1)+1),
                z=t(t(t(resum))),
                xlab="",ylab="",xaxt="n",yaxt="n",main="Summarized data",col=couleur,breaks=brek)
        }
        abline(v=ColBreak,col="purple")
        abline(h=n-RowBreak+1,col="purple")
      }
      title(outer=TRUE,main=paste(as.character(threshold),"%",sep=""))
    }
  }
  return(list(ColBreak=unique(c(1,ColBreak,d)),RowBreak=unique(n-c(n,RowBreak[length(RowBreak):1],1)+1),resum=resum,z=z,w=w,y=y,emplz=emplz,emplw=emplw, mean = mean, var=var))
}
