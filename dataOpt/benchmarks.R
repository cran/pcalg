# This script creates a graph with the timings of idaFast
timings <-  read.csv("timings.csv", stringsAsFactors=FALSE)


fast <- timings[which(timings$optimized==TRUE & timings$Myc == "Low"), ]
slow <- timings[which(timings$optimized==FALSE & timings$Myc == "Low"), ]



png("timings.png", width=400, height=400, units ="px")

r_slow <- lm(log(slow[,2]) ~ log(slow[,1]) )
r_fast <- lm(log(fast[,2]) ~ log(fast[,1]) )


plot(slow[,1:2],col="red", xlab = "Anzahl von Genen", 
   ylab = " Laufzeit (Sekunden)", pch=16, cex=1.4, xaxt= "n")
points(fast[,1:2],col="blue",pch= 8, lwd=2)

axis(1,slow[,1],slow[,1])
legend("topleft",legend = c("urspruengliche Funktion",
       "optimierte Funktion"), pch = c(16,8), 
       pt.cex = c(1.4,1), pt.lwd=c(1,2), col = c("red","blue"))
       #lwd= c(1,2), cex=c(1.4,1),col = c("red","blue"))

mtext("200", at =c(250), line = -20.55 )

dev.off()

png("timings_log.png", width=400, height=400, units ="px")

x <- c(1:5000)
y <- x^(r_slow$coefficients[2])*exp(r_slow$coefficients[1])
plot(slow[,1:2],col="red", xlab = "Number of genes, n", 
   ylab = " Time, t (sec)", ylim = c(0.01,50000), pch=16, cex=1.4, yaxt= "n", log = "xy")
lines(x,y,lty=2, col="lightgray",lwd=2)
x <- c(1:5000)
y <- x^(r_fast$coefficients[2])*exp(r_fast$coefficients[1])
lines(x,y,lty=2, col="lightgray",lwd=2)
points(fast[,1:2],col="blue",pch= 8, lwd=2)
axis(2,c(0.01,0.1,1,10,100,1000,10000),c("0.01","0.1","1","10","100","1000","10000"), las = 2)
legend("topleft",legend = c("idaFast",
       "idaFastOpt"), pch = c(16,8), 
       pt.cex = c(1.4,1), pt.lwd=c(1,2), col = c("red","blue"))
       #lwd= c(1,2), cex=c(1.4,1),col = c("red","blue"))
slope <- sprintf("%.02f", r_slow$coefficients[2])
text(200,10,bquote( 't ~ '~ O(n^.(slope)) ),srt=38)
slope <- sprintf("%.02f", r_fast$coefficients[2])
text(200,0.17,bquote( 't ~ '~ O(n^.(slope)) ), srt=25)

dev.off()
