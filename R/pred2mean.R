`.pred2mean`<-function(myperf){
		
		if(is.list(myperf)){
			tmp.myperf<-myperf
			tmp.myperf.new<-NULL
			for(j in 1:length(tmp.myperf)){
				myperf<-tmp.myperf[[j]]
				
				tmp.names<-intersect(names(myperf),names(myperf))
				myperf.new<-rep(0,length(tmp.names))
				names(myperf.new)<-tmp.names
				for(i in 1:length(tmp.names)){
					myperf.new[i]<-mean(myperf[which(names(myperf)==tmp.names[i])])
				}
				tmp.myperf.new<-rbind(tmp.myperf.new,(myperf.new))
				
			}
			myperf.new<-tmp.myperf.new
			rownames(myperf.new)<-paste("fold", 1:length(tmp.myperf), sep=".")
		}else{
			tmp.names<-intersect(names(myperf),names(myperf))
			myperf.new<-rep(0,length(tmp.names))
			names(myperf.new)<-tmp.names
			for(i in 1:length(tmp.names)){
				myperf.new[i]<-mean(myperf[which(names(myperf)==tmp.names[i])])
			}
		}
		return(myperf.new)
	}
