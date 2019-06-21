# 06/21/19 Update and test in Win10 
if(!require(pksensi)) install.packages("pksensi"); library(pksensi)

# Set and check the GCC compiler 
PATH = "c:/Rtools/mingw_32/bin"
Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";"))
Sys.which("gcc")
system('gcc -v')

mcsim_install(mxstep = 5000)

mName <- "pbtk1cpt_v2"
compile_model(mName, application = "mcsim",　version = "6.1.0")
parms <- c(vdist = 0.5, ke = 0.2, km = 0.5, kgutabs = 2.0)
params <- names(parms)

# Generate parameter matrix
LL <- 0.5 
UL <- 1.5
q <- "qunif"
q.arg <- list(list(min = parms["vdist"] * LL, max = parms["vdist"] * UL),
              list(min = parms["ke"] * LL, max = parms["ke"] * UL),
              list(min = parms["km"] * LL, max = parms["km"] * UL),
              list(min = parms["kgutabs"] * LL, max = parms["kgutabs"] * UL)) 
set.seed(1234)
x <- rfast99(params, n = 800, q = q, q.arg = q.arg, replicate = 20)

cex <- 0.2
png(file="png/parameter_matrix.png",width=3300,height=2100,res=300)
par(mfrow=c(4,4),mar=c(0,0,0,0),oma=c(3,3,2,1));
plot(x$a[,1,"vdist"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,2,"vdist"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,3,"vdist"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,4,"vdist"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,2,"ke"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,3,"ke"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,4,"ke"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,2,"km"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,3,"km"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,4,"km"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], ylab = "", cex = cex)
plot(x$a[,2,"kgutabs"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,3,"kgutabs"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,4,"kgutabs"], ylab = "", yaxt="n", cex = cex)
dev.off()

# PK modeling (decoupling simulation)
Outputs <- c("Agutlument", "Aelimination", "Acompartment", "Ccompartment", "AUC", "Ametabolized")
times <- seq(from = 0.01, to = 24.01, by = 1)
conditions <- c("Agutlument = 10") # Set the initial state of Agutlument = 10 
out <- solve_mcsim(x, mName = mName, 
                   params = params,
                   vars = Outputs,
                   time = times,
                   condition = conditions)


# Uncertainty analysis
png(file="png/uncertainty.png",width=3300,height=1800,res=300)
par(mfrow = c(2,3), mar = c(2,2,2,1), oma = c(2,2,0,0))
pksim(out$y, vars = "Agutlument", main = "Agutlument")
pksim(out$y, vars = "Aelimination", legend = F, main = "Aelimination")
pksim(out$y, vars = "Acompartment", legend = F, main = "Acompartment")
pksim(out$y, vars = "Ccompartment", legend = F, main = "Ccompartment")
pksim(out$y, vars = "Ametabolized", legend = F, main = "Ametabolized")
pksim(out$y, vars = "AUC", legend = F, main = "AUC")
mtext("Time", SOUTH<-1, line=0.4, outer=TRUE)
mtext("Quantity", WEST<-2, line=0.4, outer=TRUE)
dev.off()

# Heatmap
png(file="png/heatmap_SI.png",width=3300,height=2100,res=300)
heat_check(out)
dev.off()

png(file="png/heatmap_CI.png",width=3300,height=2100,res=300)
heat_check(out, index = "CI")
dev.off()

# Time-course sensitivity and convergene indices
png(file="png/Ccompartment.png",width=3300,height=1800,res=400)
plot(out, var = 4)
dev.off()

png(file="png/Ametabolized.png",width=3300,height=1800,res=400)
plot(out, var = 6)
dev.off()

# Scatter plot of parameter influence on model ouptuts
r <- 1 # specific replication
var <- "Ccompartment"
png(file="png/scatter_Ccompartment.png",width=3300,height=2100,res=300)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(3,3,2,1));
plot(x$a[,r,"vdist"], out$y[,r,"0.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"0.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"0.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"0.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"1.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"1.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"1.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"1.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"2.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"2.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"2.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"2.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"4.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"4.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"4.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"4.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"12.01",var], ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"12.01",var], ylab = "", yaxt="n", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"12.01",var], ylab = "", yaxt="n", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"12.01",var], ylab = "", yaxt="n", cex = cex)
mtext(var, NORTH<-3, line=0.4, adj=0, cex=1.5, outer=TRUE)
dev.off()

var <- "Ametabolized"
png(file="png/scatter_Ametabolized.png",width=3300,height=2100,res=300)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(3,3,2,1));
plot(x$a[,r,"vdist"], out$y[,r,"0.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"0.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"0.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"0.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"1.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"1.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"1.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"1.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"2.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"2.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"2.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"2.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"4.01",var], xaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"4.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"4.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"4.01",var], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,r,"vdist"], out$y[,r,"12.01",var], ylab = "", cex = cex)
plot(x$a[,r,"ke"], out$y[,r,"12.01",var], ylab = "", yaxt="n", cex = cex)
plot(x$a[,r,"km"], out$y[,r,"12.01",var], ylab = "", yaxt="n", cex = cex)
plot(x$a[,r,"kgutabs"], out$y[,r,"12.01",var], ylab = "", yaxt="n", cex = cex)
mtext(var, NORTH<-3, line=0.4, adj=0, cex=1.5, outer=TRUE)
dev.off()

# Decision making
check(out, SI.cutoff = 0.05, vars = "Ccompartment")
check(out, SI.cutoff = 0.05, vars = "Ametabolized")
check(out, SI.cutoff = 0.05)

#
file.remove(c("mcsim.pbtk1cpt_v2.model.exe", "output.csv", "pbtk1cpt_v2.model.c", "setpts.out", "input.in"))

