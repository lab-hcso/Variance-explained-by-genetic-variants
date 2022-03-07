#*********************************************************************************
# Calculates the variance in liability explained by a single biallelic loci.
# PA: allele frequency of the risk allele (denoted A)
# RR1: relative risk of Aa (one risk allele) compared to aa (no risk allele)
# RR2: relative risk of AA (two risk alleles) compared to aa (no risk allele)
# K:  overall probability of disease in population
# Returns the variance explained (Vg) and the mean liability for each genotype (the overall liability is normalized to mean 0 and variance 1)
#*********************************************************************************

func.Vg <- function (PA,RR1,RR2,K) {
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0
faa= K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 
T = qnorm(1-faa) 
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)
mean.all= PAa*muAa+ PAA*muAA
Vg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
actual.Vg =  Vg/(1+Vg) 
VR = 1-actual.Vg 
actual.T = Paa*sqrt(VR)*qnorm(1-faa) + PAa*sqrt(VR)*qnorm(1-fAa) + PAA*sqrt(VR)*qnorm(1-fAA)
actual.muaa = actual.T - sqrt(VR) * qnorm(1-faa)
actual.muAa = actual.T - sqrt(VR) * qnorm(1-fAa)
actual.muAA = actual.T - sqrt(VR) * qnorm(1-fAA)

res <- list(Vg=actual.Vg,muaa=actual.muaa, muAa = actual.muAa, muAA=actual.muAA)
res

} 

##example
func.Vg(PA=0.3,RR1=1.3,RR2=1.3^2,K=0.01) 

