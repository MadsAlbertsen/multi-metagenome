### Metagenome workflow R script v.1.0 ############################################################

# Install the required packages
install.packages("vegan")           #Used for PCA/CA analysis
install.packages("plyr")            #Easy sorting and general handling of data
install.packages("RColorBrewer")    #Easy color selection
install.packages("alphahull")       #Generation of convex spaces for extraction of points

# Load the required packages
library("vegan")
library("plyr")
library(RColorBrewer)
library(alphahull)

# Load the pre-defined functions
# Calculation of basic stats of a data subset
genome.stats <- matrix(NA, nrow=0, ncol=10) 
calc.genome.stats <- function(x) c(
  "",
  sum(g.out$length),
  nrow(g.out),
  round(mean(g.out$length),1),
  max(g.out$length),
  round(sum((g.out$gc*g.out$length))/sum(g.out$length),1),
  round(sum((g.out$HPminus*g.out$length))/sum(g.out$length),1),
  round(sum((g.out$HPplus*g.out$length))/sum(g.out$length),1),
  nrow(e.out),
  length(unique(e.out$hmm.id))
)
colnames(genome.stats) <- c("name",
                            "total.length",
                            "#scaffolds", 
                            "mean.length", 
                            "max.length", 
                            "gc", "HPminus", 
                            "HPplus", 
                            "tot.ess",
                            "uni.ess"
                            )

### Read and prepare data ##########################################################################
# Read Data - all data have been generated using previous steps - R is used to combine them
HPminus <- read.csv("HPminus.scaffold.coverage.csv", header = T)               
HPplus <- read.csv("HPplus.scaffold.coverage.csv", header = T)
gc <- read.delim("assembly.gc.tab", header = T)
kmer <- read.delim("assembly.kmer.tab", header = T)
colnames(kmer)[1] = "name"
ess <- read.table("assembly.orfs.hmm.id.txt", header = F)
colnames(ess) = c("name","orf","hmm.id")
ess.tax <- read.delim("assembly.orfs.hmm.blast.tax.tab", header = F) 
colnames(ess.tax) = c("name","orf","phylum")
cons.tax <- read.delim("assembly.tax.consensus.txt", header = T)
colnames(cons.tax) = c("name","phylum","tax.color","all.assignments")

# Combine all data to one data matrix "d"
d <- as.data.frame(cbind(HPminus$Name,
                         HPplus$Reference.length,
                         gc$gc,
                         HPminus$Average.coverage,
                         HPplus$Average.coverage
                         )
                   ,row.names=F
                   )

colnames(d) = c("name",
                "length",
                "gc",
                "HPminus",
                "HPplus"
                )

d <- merge(d,cons.tax, by = "name", all = T)

# Clean up the phyla names
d$phylum<-sub("<phylum>","",d$phylum,)
d$phylum<-sub("unclassified Bacteria","TM7",d$phylum,)
d$phylum<-sub("/Chlorobi group","",d$phylum,)
d$phylum<-sub("Chlamydiae/","",d$phylum,)
d$phylum<-sub(" group","",d$phylum,)

# Combine data on essential genes 
ess<- merge(ess, d, by = "name", all.x = T)
ess<- merge(ess, ess.tax, by = c("name","orf"), all.x = T)

# Subsetting the data by scaffold length is nice for an initial overview of the data
d <- subset(d,length > 5000)                                                   
ess <- subset(ess,length > 2000)                                                   

### General assembly statistics ###################################################################

length(d$name)            #total number of scaffolds
sum(d$length)/1000000     #total length (Mbp)
mean(d$length)            #mean scaffold length (bp)
max(d$length)             #max scaffold length (bp)


### Coverage plots - Colored by GC ################################################################
# Define a good color palette and make it transparrent
gbr<-colorRampPalette(c("green","blue","orange","red"))
gcspan <- round(max(d$gc)-min(d$gc)+1)
palette(adjustcolor(gbr(70), alpha.f = 0.1))

plot(x = d$HPminus,
     y = d$HPplus,
     log="xy",
     cex = sqrt(d$length)/100,
     pch=20,
     col=d$gc-min(d$gc),
     xlim = c(7,5000), 
     ylim = c(0.005,5000),
     xlab = "Coverage HP-",
     ylab = "Coverage HP+"     
     )


# Add scaffold size legend
t<-as.character(c(200,50,10,1))
legend(x = 7,
       y = 6000,
       t,
       bty="n",
       pch=20,
       col=rgb(0,0,0,0.2),
       pt.cex=sqrt(as.integer(t)*1000)/100,
       x.intersp=4,
       y.intersp=1,
       adj=c(0.5,0.5),
       title="Length (kbp)",
       cex=0.75
       )

# Add gc legend
for (i in 1:5){
  par(new=T,plt=c(0.27,0.31,0.80,0.87))
  barplot(rep(1,gcspan), col=1:gcspan, axes=FALSE, space=0,border=NA,  horiz=TRUE)
}

axis(2,c(0,gcspan/2,gcspan), c(round(min(d$gc)),round(min(d$gc)+(max(d$gc)-min(d$gc))/2),round(max(d$gc))), las=2, cex.axis=0.75, hadj=0.25, tck=-0.15)

mtext("% GC",side=3, line = 0,cex=0.75)


### Coverage plots - Colored by phylum affiliation ################################################
dev.off()
plot(x = d$HPminus,
     y = d$HPplus,
     log="xy",
     cex = sqrt(d$length)/100,
     pch=20,
     col=rgb(0,0,0,0.05),
     xlim = c(7,5000), 
     ylim = c(0.005,5000),
     xlab = "Coverage HP-",
     ylab = "Coverage HP+"     
)

# Color by phylum level classificiation - only the 8 most abundant phyla are added
palette(brewer.pal(8,"Set1"))
points(x = d$HPminus[d$tax.color<9],
       y = d$HPplus[d$tax.color<9],
       cex = sqrt(d$length[d$tax.color<9])/100*0.7,
       col=d$tax.color[d$tax.color<9]      
)

# Add scaffold size legend
t<-as.character(c(200,50,10,1))
legend(x = 7,
       y = 6000,
       t,
       bty="n",
       pch=20,
       col=rgb(0,0,0,0.2),
       pt.cex=sqrt(as.integer(t)*1000)/100,
       x.intersp=4,
       y.intersp=1,
       adj=c(0.5,0.5),
       title="Length (kbp)"
)

# Add phylum affiliation legend
legend(x = 5000,
       y = 1,
       arrange(d[!duplicated(factor(d$tax.color)),],tax.color)[1:8,6],
       bty="n",
       pch=20,
       col=1:8,
       pt.cex=4,
       xjust=1,
       y.intersp=1.2
)

### Genome extraction ######################################################################################
# Plot the full dataset
palette(brewer.pal(8,"Set1"))
plot(x = d$HPminus, 
     y = d$HPplus, 
     log="xy", cex = sqrt(d$length)/100, 
     pch=20, 
     col=rgb(0,0,0,0.05), 
     xlim = c(7,5000),  
     ylim = c(0.005,5000), 
     xlab = "Coverage HP-", 
     ylab = "Coverage HP+"
     )

# Add phylum level classification
points(x = d$HPminus, 
       y = d$HPplus,  
       cex = sqrt(d$length)/100*0.7, 
       col=d$tax.color
       )

# Zoom on target genome
plot(x = d$HPminus, 
     y = d$HPplus, 
     log="xy", 
     cex = sqrt(d$length)/100, 
     pch=20, col=rgb(0,0,0,0.1), 
     xlim = c(55,110),  
     ylim = c(0.5,10), 
     xlab = "Coverage HP-", 
     ylab = "Coverage HP+"
     )

points(x = d$HPminus,
       y = d$HPplus,
       cex = sqrt(d$length)/100*0.7,
       col=d$tax.color,
       lwd=2
       )

# Interactively chose 6 points on the plot that include the subset of scaffolds targeted
z.def <- ahull(locator(6, type="p", pch=20), alpha=100000)                 

# Mark the defined subset on the plot
plot(z.def,add=T, col="black")

# Extract the scaffolds (g.out) and essential genes (e.out) in the defined subset
g.out <- {}
for (i in 1:nrow(d)) { if (inahull(z.def, c(d$HPminus[i],d$HPplus[i]))) g.out <- rbind(g.out,d[i,])}
e.out<-{}
for (i in 1:nrow(ess)) { if (inahull(z.def, c(ess$HPminus[i],ess$HPplus[i]))) e.out <- rbind(e.out,ess[i,])}

# Find the duplicated single copy essential genes and color them on the graph - note TIGR00436, PF01795 and PF00750 is not single copy
d.out<-e.out[which(duplicated(e.out$hmm.id) | duplicated(e.out$hmm.id, fromLast=TRUE)),] 
d.out[,c(1,3,8)]

points(d.out$HPminus,
       d.out$HPplus, 
       col=d.out$hmm.id,
       cex=3, 
       pch = 22, 
       lwd = 3
       )

# Calculate basic statistics of the scaffolds in the bin
cbind(colnames(genome.stats),calc.genome.stats())

# Do a Correspondence Analysis on the subset of scaffolds and add the results to the extracted data
ca<-scores(cca(kmer[g.out$name,2:ncol(kmer)]), choices=1:5)$sites
g.out<-cbind(g.out,ca)
e.out<-merge(e.out,g.out[,c(1,9:13)],all.x=T,by="name")
d.out<-merge(d.out,g.out[,c(1,9:13)],all.x=T,by="name")

# Store the initial results
g2g.a <- g.out
g2d.a <- z.def
g2e.a <- e.out

# Plot the 5 most imortant components of the CA analysis
palette(adjustcolor(gbr(gcspan), alpha.f = 0.5))
pairs(g.out[,9:13],
      cex = sqrt(g.out$length)/100, 
      pch=20, 
      col=g.out$gc-min(d$gc)
      )

# Plot the combination of components that gives the highest resolution
plot(x = g.out$CA1, 
     y = g.out$CA2, 
     cex = sqrt(g.out$length)/100,
     pch = 20,
     col = g.out$gc--min(d$gc),
     xlab = "CA1",
     ylab = "CA2"
     )

# Same as above but colored by phylum level classification of essential genes
plot(x = g.out$CA1, 
     y = g.out$CA2, 
     cex = sqrt(g.out$length)/100,
     pch = 20,
     col = rgb(0,0,0,0.2),
     xlab = "CA1",
     ylab = "CA2"
)

palette(brewer.pal(8,"Set1"))
points(x = e.out$CA1,
       y = e.out$CA2, 
       col = e.out$tax.color,
       cex = sqrt(e.out$length)/100*0.7, 
       lwd = 2
      )

# Interactively define a new subspace that includes the genome of interest 
z.def <- ahull(locator(6, type = "p", pch=20), alpha = 10000)   
plot(z.def,add=T, col = "black")
g.out <- {}
for (i in 1:nrow(g2g.a)) { if (inahull(z.def, c(g2g.a$CA1[i],g2g.a$CA2[i]))) g.out <- rbind(g.out,g2g.a[i,])}
e.out<-{}
for (i in 1:nrow(g2e.a)) { if (inahull(z.def, c(g2e.a$CA1[i],g2e.a$CA2[i]))) e.out <- rbind(e.out,g2e.a[i,])}

# Store the selection
g2g.b <- g.out
g2d.b <- z.def
g2e.b <- e.out

# Show the basic stats of the selection
cbind(colnames(genome.stats),calc.genome.stats())

# Find duplicated single copy genes 
d.out<-e.out[which(duplicated(e.out$hmm.id) | duplicated(e.out$hmm.id, fromLast=TRUE)),] 

# Plot the final extraction and overlay the duplicated single copy genes - note that TIGR00436, PF01795 and PF00750 are not always single copy
plot(x = g.out$CA1, 
     y = g.out$CA2, 
     cex = sqrt(g.out$length)/100,
     pch = 20,
     col = g.out$gc--min(d$gc),
     xlab = "CA1",
     ylab = "CA2"
)

points(x = d.out$CA1,
       y = d.out$CA2,
       col = d.out$hmm.id,
       cex = 3, 
       pch = 22, 
       lwd = 3
)

# Look which single copy genes were duplicated ()
d.out[,c(1,3,8)]

# Add the final bin stats to the variable "genome.stats"
genome.stats<-rbind(genome.stats,calc.genome.stats())

# Give the genome a name
genome.stats[nrow(genome.stats),1] <- "g2"

# Print the names of the scaffolds to a file
write.table(g2g.b$name,file="g2.txt",quote=F,row.names=F,col.names=F)

### Overview of extracted genomes####################################################################
# Plot the full dataset
palette(brewer.pal(8,"Set1"))
plot(x = d$HPminus, 
     y = d$HPplus, 
     log="xy", cex = sqrt(d$length)/100, 
     pch=20, 
     col=rgb(0,0,0,0.05), 
     xlim = c(7,5000),  
     ylim = c(0.005,5000), 
     xlab = "Coverage HP-", 
     ylab = "Coverage HP+"
)

# Add phylum level classification
points(x = d$HPminus, 
       y = d$HPplus,  
       cex = sqrt(d$length)/100*0.7, 
       col=d$tax.color
)

# Add a point in the center of each extracted bin
points(genome.stats[,7], 
       genome.stats[,8], 
       pch = 20, cex=4, 
       col="white"
)

points(genome.stats[,7], 
       genome.stats[,8], 
       cex=4*0.7
)

text(as.numeric(genome.stats[,7]),
     as.numeric(genome.stats[,8]),
     cex=0.5,
     font=2
)

