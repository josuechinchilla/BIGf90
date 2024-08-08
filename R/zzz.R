.onLoad <- function(libname, pkgname){
  msg <- paste0("\nBIGf90 is a wrapper package for BLUPf90 family of programs, please use the following citations:\n\n",
  "BLUPf90 Reference: \n",
  "Misztal, I., S. Tsuruta, D.A.L. Lourenco, I. Aguilar, A. Legarra, and Z. Vitezica. 2014.\n",
  "Manual for BLUPF90 family of programs: http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf")
  
  packageStartupMessage(msg)
}