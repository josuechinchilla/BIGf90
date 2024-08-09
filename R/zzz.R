.onLoad <- function(libname, pkgname){
  msg <- paste0("\nBIGf90 is a wrapper package for BLUPf90 family of programs, please use the following citations:\n\n",
 "BIGf90 Reference: \n",
  "Chinchilla-Vargas J, Taniguti C, Sandercock A, Breeding Insight Team (2024). BIGf90:
  Breeding Insight Genomics R front face to blupf90 modules. R package version 0.3.2,
  https://github.com/josuechinchilla/BIGf90 \n\n",
  "BLUPF90 Reference: \n",
  "Misztal, I., S. Tsuruta, D.A.L. Lourenco, I. Aguilar, A. Legarra, and Z. Vitezica. 2014.\n",
  "Manual for BLUPF90 family of programs: http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf \n\n",
  "The BLUPF90 programs are free for research, but their use should be acknowledged in publications. For commercial use,
please contact Ignacy Misztal (ignacy@uga.edu) or Daniela Lourenco (danilino@uga.edu). Provide your name, company, and purpose of use in the email."              
               )
  
  packageStartupMessage(msg)
}
