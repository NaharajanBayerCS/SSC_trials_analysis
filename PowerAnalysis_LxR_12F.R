##########################################################################
#    POWER ANALYSIS for general exp with one variance component
#  ---------------------------------------------------------------
# prpgram name:   power_1way_LocxRep.R  : can be used for 
# function: To calculate the location # needed to achieve a desired power when 
#         exp design is RCBD or CRD at all locations
# 	inputs:
#       expdesign= "RCBD" or "CRD" at locations
# 		  ckmean= mean of the check or a reference  
# 		  cv=   sigma.e/ckmean   (sigma.e= square root of residual variance)
#       ntrt=  number of treatment or entry in trial
#       alpha= significnat level (i.e., 0.05, 0.1)
#       maxrep= the max rep to show in the power plot
#		    trait.name=  a text string for trait name, exp name, or researcher
#
# 	output:
#      y-axis: Power from 0 to 100
#      x-axis: rep number from 2 to maxrep
# 		 each line in graph is the power for each %diff 
#            (%diff= 2.5% to 25% by 2.5%  of ckmean)
#   written by Grace Liu                           Feb 15, 2016
###########################################################################
##  expdesign="F21"; minrep=10;  maxrep=60;  cv=8;  ckmean=100;   rLGvsE=50; nsub=3;
## ntrt1=5;  ntrt2=3; alpha=0.05; mindp=5; maxdp=30;  bydp=5 ; trait.name="YLD"; 
##

library(openxlsx)

power.12F.LocxRep<-function(expdesign, ckmean, cv, rLGvsE, nsub, ntrt1, ntrt2, 
                      alpha, minrep, maxrep, mindp, maxdp, bydp,  trait.name){
## -------  Input data dependent variables  -----------
  today=Sys.Date()
  setDesignName(expdesign)
	nrep=seq(2,maxrep,1	)    
  s2e=(cv*ckmean/100.0)**2
  s2ge=s2e*rLGvsE/100.0
	pcnt=seq(mindp, maxdp, bydp)  # in terms of %
	pcnt.txt=paste0(mindp," to ",maxdp, " by ",bydp)
  diff=pcnt*ckmean/100

  power=matrix(0,length(nrep),length(pcnt))
  
  switch(expdesign,
         ##  for 1-factor RCBD
         F11 = {
           ntrt2=0
           print ("Determine Loc Num:1 factor RCBD")
           ## --------  calculate power for Loc No. 1Factor --------------------
           if (nsub >=2) { 
            for (i in 1: length(nrep)){
                df=nrep[i]*(ntrt1-1)*(nsub-1)
              
               for (j in 1: length(pcnt)){
                 sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                 power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
             
           }else{   ### single rep cases  ###
             s2ge=0.0 ; nsub=1;
             for (i in 1: length(nrep)){
             df=(nrep[i]-1)*(ntrt1-1)
               
               for (j in 1: length(pcnt)){
                 sed=sqrt(2*s2e/nrep[i])
                 power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
           } 
         },
         
         F12 = {
           ntrt2=0
           print ("Determine Loc Num:1 factor CRD")
           ## --------  calculate power for Loc No. 1Factor --------------------
           if (nsub >=2) { 
             power=matrix(0,length(nrep),length(pcnt))	
             for (i in 1: length(nrep)){
                       df=nrep[i]*ntrt1*(nsub-1)
               
               for (j in 1: length(pcnt)){
                 sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                 power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
             
           }else{   ### single rep cases  ###
             s2ge=0.0 ; nsub=1;
             power=matrix(0,length(nrep),length(pcnt))	
             for (i in 1: length(nrep)){
              df=(nrep[i]-1)*ntrt1
               
               for (j in 1: length(pcnt)){
                 sed=sqrt(2*s2e/nrep[i])
                 power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
           } 
         },
         
         
        F21 = {
           print ("Determine Loc Num: 2-factor RCBD")
          ### when rep/1Loc >=1  ###
          if (nsub >=2 ) {
            for (i in 1: length(nrep)){
               df=nrep[i]*(ntrt1*ntrt2-1)*(nsub-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
            
          } else { 
            ### when rep/1Loc =1  ###
            s2ge=0.0; nsub=1;
            for (i in 1: length(nrep)){
             df=(nrep[i]-1)*(ntrt1*ntrt2-1)
             
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
          }
         },  
        
        F22 = {
          print ("Determine Loc Num: 2-factor SPD")
          ### when rep/1Loc >=1  ###
          if (nsub >=2 ) {
            for (i in 1: length(nrep)){
              df=nrep[i]*ntrt1*(ntrt2-1)*(nsub-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
            
          } else { 
            ### when rep/1Loc =1  ###
            s2ge=0.0; nsub=1;
            for (i in 1: length(nrep)){
              df=(nrep[i]-1)*ntrt1*(ntrt2-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
          }
        },  
        
        F23 = {
          print ("Determine Loc Num: 2-factor GUBD")
          ### when rep/1Loc >=1  ###
          if (nsub >=2 ) {
            for (i in 1: length(nrep)){
              df=nrep[i]*ntrt1*(ntrt2-1)*(nsub-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
            
          } else { 
            ### when rep/1Loc =1  ###
            s2ge=0.0; nsub=1;
            for (i in 1: length(nrep)){
              df=(nrep[i]-1)*ntrt1*(ntrt2-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
          }
        },  
        
        
        F24 = {
          print ("Determine Loc Num: 2-factor STRIP Design")
          ### when rep/1Loc >=1  ###
          if (nsub >=2 ) {
            for (i in 1: length(nrep)){
              df=nrep[i]*(ntrt1-1)*(ntrt2-1)*(nsub-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
            
          } else { 
            ### when rep/1Loc =1  ###
            s2ge=0.0; nsub=1;
            for (i in 1: length(nrep)){
              df=(nrep[i]-1)*(ntrt1-1)*(ntrt2-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
          }
        },  
        
        
        F25 = {    ##??? df  checking  ###
          print ("Determine Loc Num: 2-factor B nested in A RCBD")
          ### when rep/1Loc >=1  ###
          if (nsub >=2 ) {
            for (i in 1: length(nrep)){
              df=nrep[i]*(ntrt1*ntrt2-1)*(nsub-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
            
          } else { 
            ### when rep/1Loc =1  ###
            s2ge=0.0; nsub=1;
            for (i in 1: length(nrep)){
              df=(nrep[i]-1)*(ntrt1*ntrt2-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
          }
        }, 
        
        F26 = {
          print ("Determine Loc Num: 2-factor AxB CRD")
          ### when rep/1Loc >=1  ###
          if (nsub >=2 ) {
            for (i in 1: length(nrep)){
              df=nrep[i]*ntrt1*ntrt2*(nsub-1)
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i]/nsub+2*s2ge/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
            
          } else { 
            ### when rep/1Loc =1  ###
            s2ge=0.0; nsub=1;
            for (i in 1: length(nrep)){
              df=(nrep[i]-1)*ntrt1*ntrt2
              
              for (j in 1: length(pcnt)){
                sed=sqrt(2*s2e/nrep[i])
                power[i,j]=1-pt(qt(1-alpha/2.0,df),df, diff[j]/sed)	}  }
          }
        }  
        
  
             # if (expdesign=="RCBD2") {df=nrep[i]*(ntrt1*ntrt2-1)*(nsub-1)}
             # else if (expdesign=="CRD2") {df=nrep[i]*ntrt1*ntrt2*(nsub-1)}
             # else if (expdesign=="SPLIT") {df=nrep[i]*ntrt1*(ntrt2-1)*(nsub-1)}
             # else if (expdesign=="STRIP") {df=nrep[i]*(ntrt1-1)*(ntrt2-1)*(nsub-1)}
             # else if (expdesign=="GUBD2") {df=nrep[i]*ntrt1*(ntrt2-1)*(nsub-1)}
              

            ### when rep/1Loc =1  ###
           # s2ge=0.0; nsub=1;
      
             # if (expdesign=="RCBD2") {df=(nrep[i]-1)*(ntrt1*ntrt2-1)}
             # else if (expdesign=="CRD2") {df=(nrep[i]-1)*ntrt1*ntrt2}
             # else if (expdesign=="SPLIT") {df=(nrep[i]-1)*ntrt1*(ntrt2-1)}
              #else if (expdesign=="STRIP") {df=(nrep[i]-1)*(ntrt1-1)*(ntrt2-1)}
              #else if (expdesign=="GUBD2") {df=(nrep[i]-1)*ntrt1*(ntrt2-1)}
              
  )  
      
difflabel = paste0(diff, " (", pcnt, "%)")


###  ---------------   plot power & rep by diff%   ---------------  ###

par(oma=c(1,1,3,8))
ns = seq(minrep, maxrep, 1	)
power0 = power[(minrep-1):(maxrep-1),]

if (nsub >=2) { 
gp1=matplot(ns, power0, type = "l", xlab = paste0("Location No.  (given ",nsub," reps per Loc)"), col=c(length(pcnt):1),
            lty=1, lwd=2,  ylab = "Power", cex = 1.2, bty = "n", axes=FALSE)
} else { 
gp1=matplot(ns, power0, type = "l", xlab = paste0("Location No. (of single rep)"), col=c(length(pcnt):1),
              lty=1, lwd=2,  ylab = "Power", cex = 1.2, bty = "n", axes=FALSE)  
}

axis(side=1, at=seq(minrep, maxrep, by=2))
axis(side=2, at=seq(0, 1, by=0.1))
title( main = paste0( "Power Curve for ", var=trait.name,
                      " with ", analyname, " Design", "\n CV=", cv, "%, Mean=", ckmean, 
                      ",  ", ntrt1," Factor1,  ", ntrt2, "  Factor2" ),  
              cex = 1.2, col = 4)
abline( h = 0.8, v = nrep, col = 'gray60', lty = 2)

par(fig=c(0,1,0,1),oma=c(0,3,0,0), mar=c(0,0,0,0), new=TRUE)
legend("topright", title="Diff (Diff%)", rev(difflabel),
       col=c(1:length(pcnt)), lty=1, lwd=2,
       bty = "n", y.intersp = 1.5, xpd = TRUE)
#dev.off()
 # gp1


###  ----------------    output power table   ------------------------  ###

out = cbind(nrep, power)
out = data.frame(out)
out1 = reshape(out,
               v.names = "Power",
               timevar = "Diff",
               idvar = "nrep",
               varying = 1:length(pcnt)+1,
               direction = "long",
               sep = "")

out1$Diff = factor(out1$Diff, labels = difflabel)
names(out1)[names(out1)=="nrep"] = "Rep"

out2 = data.frame(Diff = difflabel)
power7 = subset(out1, out1$Power >= 0.7)
power8 = subset(out1, out1$Power >= 0.8)
power9 = subset(out1, out1$Power >= 0.9)
power7 = power7[!duplicated(power7$Diff), 1:2]
power8 = power8[!duplicated(power8$Diff), 1:2]
power9 = power9[!duplicated(power9$Diff), 1:2]

out2$P7 = power7[match(difflabel, power7$Diff),"Rep"]
out2$P8 = power8[match(difflabel, power8$Diff),"Rep"]
out2$P9 = power9[match(difflabel, power9$Diff),"Rep"]
out2[is.na(out2)] = paste0("> ",maxrep)

colnames(out2) = c("Diff_pcntDiff",  "N_Power70", "N_Power80", "N_Power90")
out2=cbind(out2, "CV"=cv, "Mean"=ckmean, "Rep"=nsub)

###  ----   output power table to Excel ------###
###  ----   output power table to Excel ------###
rep=paste0(nrep)
colnames(power) = paste0("Diff_",pcnt,"pcnt")
round(power, digits=4)
power1=cbind(rep, power)

info1=cbind("Exp_Design"= analyname,  "Trait_Name"=trait.name, "CV"=cv, "Ratio_LxT_to_MSE"=rLGvsE, "Rep_per_Loc"=nsub, 
              "Levels_Factor1"=ntrt1, "Levels_Factor2"=ntrt2, "Alpha_TypeI"=alpha)
info2=t(info1)
info1=data.frame("Input_Parameter"=rownames(info2), "Input_Value"=info2[,1])

outtabs=list("Input_Info"=info1, "Rep_table"=data.frame(out2), "Power_table"=data.frame(power1))
xls_name= paste0("Report_Power_N_Table_w",nsub,"Rep_",expdesign,"_",trait.name,"_CV", cv,"_", today,".xlsx") 

write.xlsx( outtabs, file=xls_name, row.names=FALSE, col.names=TRUE)    

outall=list("plot1"=gp1,"tab1"=data.frame(out2),"expInfo"=info1,"xls_name"=xls_name)
return(outall)
 
   }

##################################################################################
####---------------------------------------------------------------------####
###   function  setDesignName   to generate text "expdesign"   for         ##
####       the experiemtal design by analysis_type to the output          ###
#### --------------------------------------------------------------------####
setDesignName <- function(expdesign){
  ###Specify the fixed, sparse and random formulas:
  switch(expdesign,
         #######   1-factor designs    #####

         F11 = {
           assign("analyname","Design: Multi LOC MULTI REP RCBD 1-FACTOR", envir = .GlobalEnv)
         },
         F12 = {
           assign("analyname","Design: Multi LOC MULTI REP CRD 1-FACTOR", envir = .GlobalEnv)
         },
         
#######   2-factor designs    #####
         F21 = {
           assign("analyname","Design: Multi LOC MULTI REP RCBD AxB", envir = .GlobalEnv)
         },
         
   
         F22 = {
           assign("analyname","Design: Multi LOC MULTI REP SPLIT PLOT AxB", envir = .GlobalEnv)
         },
         
         F23 = {
           assign("analyname","Design:Multi LOC MULTI REP GUBD B nested in A ", envir = .GlobalEnv)
         },

         
         F24 = {
           assign("analyname","Design: Multi LOC MULTI REP BLUE STRIP AxB", envir = .GlobalEnv)
         }, 

        F25 = {
          assign("analyname","Design: Multi LOC MULTI REP RCBD B nested in A", envir = .GlobalEnv)
        },

        F26 = {
          assign("analyname","Design: Multi LOC MULTI REP CRD AxB", envir = .GlobalEnv)
        }
         
  )
}