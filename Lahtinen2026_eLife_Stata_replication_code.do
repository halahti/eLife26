*This is the replication code for the data preparation and analysis of the article 
*Lahtinen H., Ganna A., Kaprio J., Korhonen K., Lombardi S., Silventoinen K., Martikainen P. Heterogenous associations of polygenic indices of 35 traits with mortality. Published in eLife.
*To be published in  eLife
*Software used: Stata 16 & 18
*see PLINK commands for genetic data preparation in another code file of the github repository https://github.com/halahti/eLife26/
*Hannu Lahtinen 12 January 2026


******************
*Data preparation*
******************

frame create repopgi
frame create folk
frame create ft
frame create twsb
frame create survey
frame create vl
frame create sai
frame create pop
frame create toldtwsurv
 frame create hos
frame create med


frame sai { 
     use "D:\e45\custom-made\u1543_tk_korjattu_06_2025\U1543_PYSYVAT_TIDOT_S.dta"  	
}



 frame ft {
	use "D:\e45\custom-made\u1543_finterv_17_tk_s1.dta"
	destring *, replace
	bysort thl_id: egen koul=max(floor(ututku/1e5))
	replace koul=1 if koul==.

	keep thl_id	koul 
	gen bbID=thl_id	
	duplicates drop
	
 merge 1:1 bbID using "D:\e45\external\U1543_biopankki\U1543_BB2020_8_Original_data\tutkpalv_u1543_bb2020_8_fh17_phenotypes_09072020.dta"
 destring FT17_OTOS_SUKUP, gen(sukup) ig(NA)
 
	 destring IKA_POIMINTA, gen(ika) ig(NA)
*ika_poiminta is the age 16.11.2016
	gen ikav=floor(ika)
	gen ikapvd=ika-ikav
		gen ikapvd2=ikapvd+45/365
		replace ikapvd2=ikapvd2-1 if ikapvd2>1
egen syntykk=cut(ikapvd2), gr(12)
replace syntykk=syntykk+1

gen ikav16=floor(ika+45/365)
gen syntyv=2016-ikav16

 }
 

	 
	 
	 
 
 
 frame repopgi {
  use "D:\e45\external\tutkpalv_u1543_pgirepo_frfhh2000_1.dta"
	 append using  "D:\e45\external\U1543_pgirepo_tw_kor_2025_06", gen(tw)
	 
	 
 }
 
 
 
 
    frame survey {
    append using "W:\Hannu\numericdata\fr92.dta" "W:\Hannu\numericdata\fr97.dta" "W:\Hannu\numericdata\fr02.dta" "W:\Hannu\numericdata\fr07.dta" "W:\Hannu\numericdata\fr12.dta" "W:\Hannu\numericdata\t2000.dta" "W:\Hannu\numericdata\t2011.dta" "W:\Hannu\numericdata\ft17.dta",  gen(aineisto) 
 
 
 
 
rename vuosi aineistovuosi 
 
replace aineistovuosi=2000 if aineisto==6
replace aineistovuosi=2011 if aineisto==7
replace aineistovuosi=2017 if aineisto==8

gen thl_id=bbID

duplicates tag thl_id, gen(d)

*for health2000/2011, the baseline year will be 2000 (except for the additional sample for second round)
drop if d==1 & aineistovuosi==2011

*some variables from "an extension block"
merge 1:1 bbID using "W:\Hannu\numericdata\ft17_edt.dta",
rename Q40 Q40o
merge 1:1 bbID using "W:\Hannu\numericdata\fr12_edt.dta", gen(_merge2)

merge 1:1 bbID using "W:\Hannu\numericdata\t2000_edt.dta", gen(_merge3)
merge 1:1 bbID using "W:\Hannu\numericdata\t2011_edt.dta", gen(_merge4) update 



*******risk factor phenotypes

gen alkot2000=  KYS1_K39_M2
 gen alkot2011=  KYS1_K39_M4
replace alkot2011=alkot2011/52

 gen alkofr9297=ALKI2
 gen alkofr02=ALKI2_FR02 
gen alkofr0712=ALKI2_FR07




egen alkog=rowtotal(alkofr9297 alkofr0712 alkofr02 alkot2000 alkot2011)
replace alkog=. if alkofr9297==. & alkofr0712==. & alkofr02==. & alkot2000==. & alkot2011==. 


egen alkg17=rowtotal(FT17_L1_82_1 FT17_L1_82_2 FT17_L1_82_3 FT17_L1_82_4)
replace alkg17=. if aineistovuosi!=2017


gen alcfreq17=. 
replace alcfreq17=0 if FT17_L1_79==0
replace alcfreq17=1/4 if FT17_L1_79==1
replace alcfreq17=1 if FT17_L1_79==2
replace alcfreq17=3 if FT17_L1_79==3
replace alcfreq17=7 if FT17_L1_79==4

gen alcdose17=.
replace alcdose17=2*12 if FT17_L1_80==1
replace alcdose17=4 *12 if FT17_L1_80==2
replace alcdose17=6 *12 if FT17_L1_80==3
replace alcdose17=9 *12 if FT17_L1_80==4
replace alcdose17=15 *12 if FT17_L1_80==5

gen alcg17=alcfreq17*alcdose17
replace alcg17=0 if  FT17_L1_79==0

replace alkog=alcg17 if aineistovuosi==2017

sum alkog,d
replace alkog=r(p99)  if alkog>r(p99) & alkog!=.

gen smoke=.
replace smoke=1 if TUPI3==1
replace smoke=2 if TUPI3==2
replace smoke=3 if TUPI3==3 | TUPI3==4


replace smoke=1 if FB01==0 | FB02==0 | FB03==0
replace smoke=2 if (FB02==1 | FB03==1)  & (FB06==4 |FB06==5 |FB06==6 |FB06==7 )
replace smoke=3 if (FB02==1 | FB03==1)  & (FB06==1 |FB06==2 |FB06==3 |FB06==.)


lab def smoke 1 "Never" 2 "Quitted >6 months  ago" 3 "Yes", replace
lab val smoke smoke

replace bmi=BMII_BMI if bmi==.
replace  bmi=FT17_TT2_12_BMI if bmi==.

gen selfhealth=Q40o
replace selfhealth=Q40 if aineistovuosi==2012
replace selfhealth=FT17_L1_1 if aineistovuosi==2017

replace selfhealth=BA01 if aineistovuosi==2000 |aineistovuosi==2011



*depression questions: fr(pl1992): KY62 KY60_6 ky63 H2000/11: KYS1_K77 KYS1_K1304  KYS1_K0105 KYS3_MASE Kys1_K1035 FR92: Q45J. K124_3; ft17: M_BDI6_SUM M_BDI6_2LK,  FT17_L1_5_7
*FR02_26M
gen survdep=0 if aineistovuosi!=.
replace survdep=1 if KY62==2

replace survdep=1 if Q45J==3 | K124_3<5
replace survdep=1 if M_BDI6_2LK==1 |  FT17_L1_5_7==2
replace survdep=1 if KYS1_K77==3 |MDD12==3 | CIDI_MDD12==1 | KYS1_K0105==2 | KYS1_K0105==3 | KYS3_MASE==3 |KYS3_MASE==4 | KYS3_MASE==5
replace survdep=. if KY62==. &  Q45J==. & K124_3==. & M_BDI6_2LK==. &  FT17_L1_5_7==. &  KYS1_K77==. & MDD12==. & CIDI_MDD12==. & KYS1_K0105==.  & KYS1_K0105==. & KYS3_MASE==.
   }
   
   
   

  
  
 
 frame folk {
 		use "D:\e45\custom-made\u1543_tk_korjattu_06_2025\U1543_FOLK_PERUS_LAAJA_S.dta", clear
		append using "D:\e45\custom-made\u1543_tk_korjattu_06_2025\U1543_FOLK_VL_LAAJA_S.dta", force 
	keep sukup syntyv kieli kuolv thl_id ututku_aste  
		keep if syntyv!=.

	destring *, replace
	
	replace ututku_aste=1 if  ututku_aste==.
	bysort thl_id: egen koul=max(ututku_aste)
		
	bysort thl_id: egen kuolv2=min(kuolv)	
	drop kuolv
	rename kuolv2 kuolv
	
	bysort thl_id: egen kieli=min(kieli_k)	
	
	drop ututku_aste kieli_k
	
	duplicates drop
compress


*two individuals in 80s data blocks had inconsitent sex to newer block  (both are women in the older, probably less reliable, block, and men in the newer)
duplicates tag thl_id, gen(d)
 drop if d==1 & sukup==2
 drop d
 }
 
 
frame  toldtwsurv {
	use "D:\e45\external\tutkpalv_U1543_twins\U1543_twin1975_survey_labelled_kor_2025_06", clear
	append using "D:\e45\external\tutkpalv_U1543_twins\U1543_twin1975_survey_labelled_kor_2025_06.dta"
gen thl_id=bbid
keep thl_id sex
duplicates drop
} 
 
 
 
 
 ***population data
 
 
 cwf pop
 use "D:\e45\external\u1543_fr_ft_h2000_tw_pca_kor_2025_06.dta"
 
 frlink 1:1 thl_id, frame(survey)
 frget aineistovuosi bmi smoke alkog selfhealth survdep, from(survey)
 merge 1:1 thl_id using "D:\e45\external\external_korjaus_2025_06\U1543_KSYYT_4725_FINAL_S.dta", gen(_merge) keep(master match)
  
frlink 1:1 thl_id, frame(folk)
frget  sukup kieli  koul  , from(folk)
 
 frlink 1:1 thl_id, frame(toldtwsurv)
frget  sex, from(toldtwsurv)
 
frlink 1:1 thl_id, frame(sai)
frget  syntyv syntykk  maastamuuttoeka maastamuuttovika, from(sai) 

 
 frlink 1:1 thl_id, frame(repopgi)
frget  pgs* , from(repopgi)

  frlink 1:1 thl_id, frame(ft)
 frget koul sukup syntyv syntykk , from(ft) suffix("_ft")
 replace koul=koul_ft if koul==.
   replace sukup=sukup_ft if sukup==.
   replace sukup=sex if sukup==.
   
   replace syntyv=syntyv_ft if syntyv==.
   replace syntykk=syntykk_ft if syntykk==.

   gen svkk=syntyv+syntykk/12-(1/24)

*school years after basic
gen koulv=.
replace koulv=0 if koul==1
replace koulv=3 if koul==3
replace koulv=4 if koul==4
replace koulv=5 if koul==5
replace koulv=6 if koul==6
replace koulv=8 if koul==7
replace koulv=12 if koul==8

	  gen female=sukup-1
  
	gen koul4=koul if koul==1
replace koul4=2 if koul==3 |koul==4
replace koul4=3 if koul==5 |koul==6
replace koul4=4 if koul==7 |koul==8

gen koul3=koul4
replace koul3=3 if koul4==4

gen ika95jan1=1995-svkk
  
  
 gen svkk2=dofm(svkk)
gen day25=svkk2+(25*365.25)

gen svkk3=ym(syntyv,syntykk)
gen svkk4=dofm(svkk3)

gen kdum=1 if kvuosi!=.
replace kdum=0 if kvuosi==.

gen twcohort=1     if   aineistovuosi==. & syntyv < 1962
replace twcohort=2 if   aineistovuosi==. & syntyv >1962 & syntyv<1980
replace twcohort=3 if   aineistovuosi==. & syntyv >1980 & syntyv<1990

gen startyr=1995
replace startyr=aineistovuosi+0.5 if aineistovuosi>1995 & aineistovuosi!=.

replace startyr=2004 if twcohort==2
replace startyr=2013 if twcohort==3
replace startyr=svkk+25 if ((aineistovuosi+0.5)-svkk<25) & (aineistovuosi>1970 ) & (svkk+25>1995 )
replace startyr=svkk+25 if twcohort==2 & syntyv==1979

replace startyr=. if startyr>=2020
gen sm1=round((startyr-floor(startyr))*12+0.6)
gen startmonth=ym(floor(startyr),sm1)

gen startday=dofm(startmonth)

replace startday=. if startday>kpaiv


gen risktime=kpaiv-startday if kpaiv>startday
replace risktime = td(31dec2019)-startday   if kpaiv==.

replace risktime=. if risktime<0

gen startage=round((startday-svkk4)/365.25,.01)  
gen startagedays=(startday-svkk4) if startyr!=.

gen endday=startday+risktime
gen endage=round((startagedays+risktime)/365.25,.01) if startyr!=.
gen endagedays=startagedays+risktime
  
replace endday=. if endday<startday 
  replace endage=. if endday<startday 
    
forvalues i=1/20{
	sum pc`i'
	replace pc`i'=(pc`i'-r(mean))/r(sd)
	
}
  
*hospital admissions and medicines, needed for depression
 frame hos {
use "D:\e45\external\external_korjaus_2025_06\U1543_THL_4725_96_20_S.dta", clear

gen tupv=dofc(tupva)
gen lapv=dofc(lpvm)
drop tupva lpvm

append using "D:\e45\external\external_korjaus_2025_06\U1543_THL_4725_94_95_S.dta", gen(k94)  force
append using "D:\e45\external\external_korjaus_2025_06\U1543_THL_4725_87_93_S.dta", gen(k93) 

replace  tupv=tupva if k94==1 |k93==1
replace  lapv=lpvm if k94==1 |k93==1

gen tuyr=year(tupv)
gen layr=year(lapv)

keep if layr>=1993 & tuyr<=2019

frlink m:1 thl_id, frame(pop)
frget startyr startday, from(pop)

gen kee=1 if (startday-tupv<731 & startday-tupv>0) 
replace kee=1 if  (startday-lapv<731 & startday-lapv>0)
replace kee=1 if  (startday-lapv<0 & startday-tupv>731)

keep if kee==1

foreach i in  pdgo pdge sdg1o sdg1e sdg2o sdg2e {
di "`i'	"
gen `i'3=	substr(`i',1,3) if k94==0 & k93==0
gen `i'dep=0
replace `i'dep=1 if `i'3 == "F32" | `i'3 == "F33" |  `i'3 == "F34" |  `i'3 == "F35" |  `i'3 == "F36" |  `i'3 == "F37" |  `i'3 == "F38" |  `i'3 == "F39" 
	replace `i'dep=0 if (`i'=="F323" | `i'=="F333") & (k94==0 & k93==0)
}
*old data (1993-1994) had ICD9 classification
foreach i in  pdg  sdg1  sdg2 sdg3 {
gen `i'3=	substr(`i',1,3) if k94==1 | k93==1
gen `i'4=	substr(`i',1,4) if k94==1 | k93==1
gen `i'dep=0
replace `i'dep=1 if  `i'3=="296" |  `i'4=="3090" |  `i'4=="3091"  |  `i'3=="311"
}

 egen hosdep0=rowmax(pdgodep pdgedep sdg1odep sdg1edep sdg2odep sdg2edep pdgodep  sdg1dep  sdg2dep sdg3dep)
egen hosdep=max(hosdep0), by(thl_id)

keep thl_id hosdep  
duplicates drop
 }

 frame med {
 use   "D:\e45\external\external_korjaus_2025_06\U1543_KELA_RES_S.dta", clear
 
 
frlink m:1 thl_id, frame(pop)
frget startyr startday, from(pop)

gen kee=1 if startday-otpvm<731 &  startday-otpvm>0
replace kee=1 if year(otpvm)<1997 & startyr<1997

  keep if kee==1
  
  gen atc4= substr(atc,1,4)
    gen atc5= substr(atc,1,5)
	keep if atc4=="N06A" | atc5=="N06CA"
	gen meddep=1 
	keep thl_id meddep
	duplicates drop
 }  
 
 
  frlink 1:1 thl_id, frame(hos)
frget hosdep, from(hos)
 
 replace hosdep=0 if hosdep==.
 
  frlink 1:1 thl_id, frame(med)
frget meddep, from(med)
 
  replace meddep=0 if meddep==.

  
  
  **sibling link
  
  frame twsb: use "W:\Hannu\families_ids\sibsz.dta"
  
  frlink 1:1 thl_id, frame(twsb)
  
frget fam_id dz, from(twsb)

merge 1:1 thl_id using "W:\Hannu\families_ids\mztws.dta", keep(master) keepusing(thl_id) nogenerate

drop if syntyv==. | risktime==.
drop if sukup==.

bysort fam_id: drop if _N==1
  
 gen pairid=fam_id
 
 
 drop if dz==. & aineisto==.

 **************************

 gen deathage=round(startage+(risktime/365.25),.001) if kdum==1
 
 
 gen kdum25=0
replace kdum25=1 if deathage>=25 & deathage<65

gen kdum65=0
replace kdum65=1 if deathage>=65 & deathage<80

gen kdum80=0
replace kdum80=1 if deathage>=80 & deathage<110

 egen startage3=cut(startage), at(25 64 79 110)
 
 
 
 gen startday25=startday if startage>=25 & startage<65
 gen startday65=startday if startage>=65 & startage<80 & kdum25==0
  gen startday80=startday if startage>=80  & kdum25==0 & kdum65==0
  
  replace startday65=svkk4+65*365.25 if startday65==. & svkk4+65*365.25<td(31dec2019)
    replace startday80=svkk4+80*365.25 if startday80==. & svkk4+80*365.25<td(31dec2019)
	

gen endday25=.
replace endday25=kpaiv if startday25!=. & kdum25==1
replace endday25=svkk4+65*365.25 if startday25!=. & kdum25==0
replace endday25=td(31dec2019) if startday25!=. & endday25>td(31dec2019)

gen endday65=.
replace endday65=kpaiv if startday65!=. & kdum65==1
replace endday65=svkk4+80*365.25 if startday65!=. & kdum65==0
replace endday65=td(31dec2019) if startday65!=. & endday65>td(31dec2019)


gen endday80=.
replace endday80=kpaiv if startday80!=. & kdum80==1
replace endday80=td(31dec2019) if startday80!=. & kdum80==0 & endday80>td(31dec2019)

replace endday25=. if startday25>endday25
replace endday65=. if startday65>endday65
replace endday80=. if startday80>endday80


gen startage25=round((startday25-svkk4)/365.25,.01)  
gen startage65=round((startday65-svkk4)/365.25,.01)  
gen startage80=round((startday80-svkk4)/365.25,.01)  


gen endage25=round((endday25-svkk4)/365.25,.01)  
gen endage65=round((endday65-svkk4)/365.25,.01)  
gen endage80=round((endday80-svkk4)/365.25,.01)  


*external and natural causes of death
gen external=0
replace external=1 if tpksaika>40 & tpksaika<54 & risktime!=. & kdum==1

gen nonexternal=0
replace nonexternal=1 if tpksaika<=40  & kdum==1

gen depsum=hosdep +survdep+meddep


*missing phenotype data
sum  bmi smoke selfhealth alkog koulv  depsum hosdep survdep meddep if  dz==. &pc1!=.

tab aineistovuosi if alkog==. & dz==. & pc1!=.
tab aineistovuosi if survdep==. & dz==. & pc1!=.
tab aineistovuosi if bmi==. & dz==. & pc1!=.
tab aineistovuosi if smoke==. & dz==. & pc1!=.
tab aineistovuosi if selfhealth==. & dz==. & pc1!=.

tab aineistovuosi if survdep==.  & alkog!=. & dz==. & pc1!=.



gen aineistovuosi2=aineistovuosi
replace aineistovuosi2=1995 if twcohort==1
replace aineistovuosi2=2004 if twcohort==2
replace aineistovuosi2=2013 if twcohort==3

qui reg kdum bmi smoke selfhealth alkog koulv depsum pgs_bmi  pgs_dep1 pgs_dpw2 pgs_ea6 pgs_eversmoke2 pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 c.startage##c.startage if  dz==.
gen phsample=e(sample)

foreach i in smoke bmi selfhealth alkog koulv depsum {
   egen `i'_std=std(`i') if phsample==1   
}

replace  selfhealth_std=selfhealth_std*-1

gen selfhealth_rev=6-selfhealth


gen bmi_cat=1 if bmi<18.5
replace bmi_cat=2 if bmi>=18.5 & bmi<25
replace bmi_cat=3 if bmi>=25 & bmi<30
replace bmi_cat=4 if bmi>=30 & bmi<35
replace bmi_cat=5 if bmi>=35 & bmi<99

gen alkog_cat=1 if alkog==0
replace alkog_cat=2 if alkog>0 & alkog<=24
replace alkog_cat=3 if alkog>24 & alkog<=72
replace alkog_cat=4 if alkog>72 & alkog<=192
replace alkog_cat=5 if alkog>192 & alkog<=999


tostring(syntyv), gen(sv_str)
tostring(syntykk), gen(sk_str)

gen syntystr=sv_str+sk_str+"15"
replace syntystr=sv_str+"0"+sk_str+"15" if syntykk<10


gen dobirth=date(syntystr,"YMD")

*
exit
*

******
******
******

*analysis code

cd "W:\Hannu\repomort\rep\"

log using "replication.smcl"

 stset endage, failure(kdum==1) enter(startage)  origin(time 0)
 
 putexcel set "mortality_pgi_jan25.xlsx", modify sheet("PGI")

****************************************************************************************
*Main analysis + within-sibship + Sex + education interactions						   *
*Figure 1 & Supplemetary table S5 + Panels A & B in figure 2 & Supplementary Table S6  *
****************************************************************************************


local j=4
foreach i in pgs_activity1 pgs_adhd1 pgs_adventure1  pgs_afb2 pgs_asteczrhi1 pgs_asthma1 pgs_audit1 pgs_bmi2 pgs_cannabis2 pgs_cp2 pgs_cpd2 pgs_dep1 pgs_dpw2 pgs_ea6 pgs_eversmoke2 pgs_extra2 pgs_famsat1 pgs_friendsat1 pgs_hayfever1 pgs_height3 pgs_highmath1 pgs_leftout1 pgs_menarche1 pgs_migraine1 pgs_morning1 pgs_narcis1 pgs_nearsighted1 pgs_nebwomen2 pgs_neuro2 pgs_open2 pgs_reading1 pgs_religatt1 pgs_risk2 pgs_selfhealth1 pgs_selfmath1 {

di "`i':"

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==., 
local popb=_b[`i']
local popse = _se[`i']
local popn= e(N)
local popd= e(N_fail)
local popr=e(risk) 
di "pop: " exp(`popb') " 95%ci: " exp(`popb' +invnormal(0.025)*`popse') " " exp(`popb' +invnormal(0.975)*`popse')  " z: " abs(`popb'/`popse')

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz!=.,  strata(fam_id)
local sibb=_b[`i']
local sibse = _se[`i']
local sibn= e(N)
local sibd= e(N_fail)
local sibr=e(risk) 
di "sib: " exp(`sibb') " 95%ci: " exp(`sibb' +invnormal(0.025)*`sibse') " " exp(`sibb' +invnormal(0.975)*`sibse')  " z: " abs(`sibb'/`sibse')

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==. & female==0 
local menb=_b[`i']
local mense = _se[`i']
local menn= e(N)
local mend= e(N_fail)
local menr=e(risk) 
di "men: " exp(`menb') " 95%ci: " exp(`menb' +invnormal(0.025)*`mense') " " exp(`menb' +invnormal(0.975)*`mense')  " z: " abs(`menb'/`mense')

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==. & female==1 
local womb=_b[`i']
local womse = _se[`i']
local womn= e(N)
local womd= e(N_fail)
local womr=e(risk) 
di "women: " exp(`womb') " 95%ci: " exp(`womb' +invnormal(0.025)*`womse') " " exp(`womb' +invnormal(0.975)*`womse')  " z: " abs(`womb'/`womse')

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==. & koul3==1 
local basb=_b[`i']
local basse = _se[`i']
local basn= e(N)
local basd= e(N_fail)
local basr=e(risk) 
di "bas: " exp(`basb') " 95%ci: " exp(`basb' +invnormal(0.025)*`basse') " " exp(`basb' +invnormal(0.975)*`basse')  " z: " abs(`basb'/`basse')

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==. & koul3==2
local secb=_b[`i']
local secse = _se[`i']
local secn= e(N)
local secd= e(N_fail)
local secr=e(risk)
di "sec: " exp(`secb') " 95%ci: " exp(`secb' +invnormal(0.025)*`secse') " " exp(`secb' +invnormal(0.975)*`secse')  " z: " abs(`secb'/`secse')

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==. & koul3==3
local ltrb=_b[`i']
local ltrse = _se[`i']
local ltrn= e(N)
local ltrd= e(N_fail)
local ltrr=e(risk)
di "ltr: " exp(`ltrb') " 95%ci: " exp(`ltrb' +invnormal(0.025)*`ltrse') " " exp(`ltrb' +invnormal(0.975)*`ltrse')  " z: " abs(`ltrb'/`ltrse')

putexcel b`j'=`popb'   c`j'=`popse'   d`j'=`sibb'   e`j'=`sibse'   f`j'=`menb'   g`j'=`mense'   h`j'=`womb'  i`j'=`womse'  j`j'=`basb'   k`j'=`basse'   l`j'=`secb'  m`j'=`secse'  n`j'=`ltrb'   o`j'=`ltrse'   
local j=`j'+1
}

putexcel b39=`popn' b40=`popd' b41=`popr' d39=`sibn' d40=`sibd' d41=`sibr' f39=`menn' f40=`mend' f41=`menr' h39=`womn' h40=`womd' h41=`womr' j39=`basn' j40=`basd' j41=`basr' l39=`secn' l40=`secd' l41=`secr' n39=`ltrn' n40=`ltrd'   n41=`ltrr' 


*number of strata
xtlogit kdum pgs_bmi2 female pc1-pc10 i.aineistovuosi2 if  dz!=. ,  i(fam_id) or
putexcel d42= (e(N_g))


***********************************************************************
*Age-specific follow-up (Panel C of Figure 2/ supplementary table S6) *
***********************************************************************

local j=4  
foreach i in pgs_activity1 pgs_adhd1 pgs_adventure1  pgs_afb2 pgs_asteczrhi1 pgs_asthma1 pgs_audit1 pgs_bmi2 pgs_cannabis2 pgs_cp2 pgs_cpd2 pgs_dep1 pgs_dpw2 pgs_ea6 pgs_eversmoke2 pgs_extra2 pgs_famsat1 pgs_friendsat1 pgs_hayfever1 pgs_height3 pgs_highmath1 pgs_leftout1 pgs_menarche1 pgs_migraine1 pgs_morning1 pgs_narcis1 pgs_nearsighted1 pgs_nebwomen2 pgs_neuro2 pgs_open2 pgs_reading1 pgs_religatt1 pgs_risk2 pgs_selfhealth1 pgs_selfmath1  {

 qui  stset endage25, failure(kdum25==1)  enter(startage25)   
	
di "`i':"

qui stcox `i' female pc1-pc10 i.aineistovuosi2  if  dz==.
local age25b=_b[`i']
local age25se = _se[`i']
local age25n= e(N)
local age25d= e(N_fail)
local age25r=e(risk)
di "age2564: " exp(`age25b') " 95%ci: " exp(`age25b' +invnormal(0.025)*`age25se') " " exp(`age25b' +invnormal(0.975)*`age25se')  " z: " abs(`age25b'/`age25se')


 qui     stset endage65, failure(kdum65==1)  enter(startage65)  

qui stcox `i' female pc1-pc10 i.aineistovuosi2   if  dz==.
local age65b=_b[`i']
local age65se = _se[`i']
local age65n= e(N)
local age65d= e(N_fail)
local age65r=e(risk)
di "age6579: " exp(`age65b') " 95%ci: " exp(`age65b' +invnormal(0.025)*`age65se') " " exp(`age65b' +invnormal(0.975)*`age65se')  " z: " abs(`age65b'/`age65se')

qui     stset endage80, failure(kdum80==1) enter(startage80)   
 

qui stcox `i' female pc1-pc10 i.aineistovuosi2   if  dz==.
local age80b=_b[`i']
local age80se = _se[`i']
local age80n= e(N)
local age80d= e(N_fail)
local age80r=e(risk)
di "age80: " exp(`age80b') " 95%ci: " exp(`age80b' +invnormal(0.025)*`age80se') " " exp(`age80b' +invnormal(0.975)*`age80se')  " z: " abs(`age80b'/`age80se')



putexcel p`j'=`age25b'   q`j'=`age25se'   r`j'=`age65b'   s`j'=`age65se'   t`j'=`age80b'   u`j'=`age80se'  
local j=`j'+1
}


putexcel p39=`age25n' p40=`age25d' p41=`age25r' r39=`age65n' r40=`age65d' r41=`age65r' t39=`age80n' t40=`age80d' t41=`age80r'

**********************************************************************************
*External/natural causes of death (Panel D of FIgure 2 / Supplementary Table S6) *
**********************************************************************************

local j=4  
foreach i in pgs_activity1 pgs_adhd1 pgs_adventure1  pgs_afb2 pgs_asteczrhi1 pgs_asthma1 pgs_audit1 pgs_bmi2 pgs_cannabis2 pgs_cp2 pgs_cpd2 pgs_dep1 pgs_dpw2 pgs_ea6 pgs_eversmoke2 pgs_extra2 pgs_famsat1 pgs_friendsat1 pgs_hayfever1 pgs_height3 pgs_highmath1 pgs_leftout1 pgs_menarche1 pgs_migraine1 pgs_morning1 pgs_narcis1 pgs_nearsighted1 pgs_nebwomen2 pgs_neuro2 pgs_open2 pgs_reading1 pgs_religatt1 pgs_risk2 pgs_selfhealth1 pgs_selfmath1 {

 qui     stset endage, failure(external==1) enter(startage)   origin(time 0)

di "`i':"

qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==., 
local extb=_b[`i']
local extse = _se[`i']
local extn= e(N)
local extd= e(N_fail)
local extr= e(risk)
di "ext: " exp(`extb') " 95%ci: " exp(`extb' +invnormal(0.025)*`extse') " " exp(`extb' +invnormal(0.975)*`extse')  " z: " abs(`extb'/`extse')

 qui     stset endage, failure(nonexternal==1)  enter(startage)   origin(time 0)


qui stcox `i' female pc1-pc10 i.aineistovuosi2 if  dz==., 
local nonexb=_b[`i']
local nonexse = _se[`i']
local nonexn= e(N)
local nonexd= e(N_fail)
local nonexr=e(risk)
di "nonex: " exp(`nonexb') " 95%ci: " exp(`nonexb' +invnormal(0.025)*`nonexse') " " exp(`nonexb' +invnormal(0.975)*`nonexse')  " z: " abs(`nonexb'/`nonexse')

putexcel v`j'=`extb'   w`j'=`extse'   x`j'=`nonexb'  y`j'=`nonexse'    
local j=`j'+1
}

putexcel v39=`extn' v40=`extd' v41=`extr' x39=`nonexn'  x40=`nonexd'   x41=`nonexr' 


******************************
*Phenotype analysis (Table 1)*
******************************

      stset endage, failure(kdum==1) enter(startage)   origin(time 0)
	  
 putexcel set "mortality_pgi_jan25.xlsx", modify sheet("PGI_pheno")
  
	  
stcox bmi_std female pc1-pc10 i.aineistovuosi2 if phsample==1
mat m1=(r(table)[1,1]\r(table)[5..6,1])'
stcox pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')

 stcox alkog_std female pc1-pc10 i.aineistovuosi2 if phsample==1
 mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
  stcox pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')

   
 stcox depsum_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
    mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
	
 stcox koulv_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
   mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')

 stcox smoke_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')

 stcox selfhealth_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat m1=m1\((r(table)[1,1]\r(table)[5..6,1])')
 
 putexcel b5=matrix(m1)
 
stcox bmi_std pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
mat m2=(r(table)[1,1..2]\r(table)[5..6,1..2])'
 stcox alkog_std pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
mat m2=m2\((r(table)[1,1..2]\r(table)[5..6,1..2])')
  stcox depsum_std pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat m2=m2\((r(table)[1,1..2]\r(table)[5..6,1..2])')
 stcox koulv_std pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat m2=m2\((r(table)[1,1..2]\r(table)[5..6,1..2])')
 stcox smoke_std pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat m2=m2\((r(table)[1,1..2]\r(table)[5..6,1..2])')
 stcox selfhealth_std pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat m2=m2\((r(table)[1,1..2]\r(table)[5..6,1..2])')

 putexcel e5=matrix(m2)

 
 stcox bmi_std alkog_std depsum_std koulv_std smoke_std selfhealth_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat m3f=(r(table)[1,1..6]\r(table)[5..6,1..6])'
 stcox pgs_bmi pgs_dpw2 pgs_dep1  pgs_ea6 pgs_eversmoke2 pgs_selfhealth1  female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat m3g=(r(table)[1,1..6]\r(table)[5..6,1..6])'

 mat m3=m3f[1,1...]\m3g[1,1...]
 forvalues i=2/6 {
   mat m3=m3\(m3f[`i',1...]\m3g[`i',1...])  
 }
 
  putexcel h5=matrix(m3)

  stcox  bmi_std pgs_bmi alkog_std pgs_dpw2 depsum_std pgs_dep1 koulv_std pgs_ea6 smoke_std pgs_eversmoke2 selfhealth_std      pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 

  mat m4=(r(table)[1,1..12]\r(table)[5..6,1..12])'

    putexcel k5=matrix(m4)
	

**************************************************** 	
*Categorical phenotypes (Supplementary table S7)   *
****************************************************
  
stcox ib2.bmi_cat female pc1-pc10 i.aineistovuosi2 if phsample==1
mat mo1=(r(table)[1,1..5]\r(table)[5..6,1..5])'
stcox c.pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
mat mo1=mo1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox ib2.alkog_cat female pc1-pc10 i.aineistovuosi2 if phsample==1
 mat mo1=mo1\((r(table)[1,1..5]\r(table)[5..6,1..5])')
  stcox c.pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
  mat mo1=mo1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox i.depsum female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat mo1=mo1\((r(table)[1,1..4]\r(table)[5..6,1..4])')
 stcox c.pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
   mat mo1=mo1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox i.koul4 female pc1-pc10 i.aineistovuosi2 if phsample==1 
    mat mo1=mo1\((r(table)[1,1..4]\r(table)[5..6,1..4])')
  stcox c.pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
     mat mo1=mo1\((r(table)[1,1]\r(table)[5..6,1])')
  stcox i.smoke female pc1-pc10 i.aineistovuosi2 if phsample==1 
     mat mo1=mo1\((r(table)[1,1..3]\r(table)[5..6,1..3])')
  stcox c.pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1
     mat mo1=mo1\((r(table)[1,1]\r(table)[5..6,1])')
 stcox i.selfhealth female pc1-pc10 i.aineistovuosi2 if phsample==1 
    mat mo1=mo1\((r(table)[1,1..5]\r(table)[5..6,1..5])')
 stcox c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
    mat mo1=mo1\((r(table)[1,1]\r(table)[5..6,1])')
 
     putexcel b20=matrix(mo1)

stcox ib2.bmi_cat c.pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
mat mo2=(r(table)[1,1..6]\r(table)[5..6,1..6])'
 stcox ib2.alkog_cat c.pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
mat mo2=mo2\((r(table)[1,1..6]\r(table)[5..6,1..6])')
 stcox i.depsum c.pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat mo2=mo2\((r(table)[1,1..5]\r(table)[5..6,1..5])')
 stcox i.koul4 c.pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat mo2=mo2\((r(table)[1,1..5]\r(table)[5..6,1..5])')
 stcox i.smoke c.pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1 
 mat mo2=mo2\((r(table)[1,1..4]\r(table)[5..6,1..4])')
 stcox i.selfhealth c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
 mat mo2=mo2\((r(table)[1,1..6]\r(table)[5..6,1..6])')
  
  putexcel e20=matrix(mo2) 
	  
 stcox ib2.bmi_cat ib2.alkog_cat i.depsum i.koul4 i.smoke i.selfhealth female pc1-pc10 i.aineistovuosi2 if phsample==1 
 mat mo3f=(r(table)[1,1..26]\r(table)[5..6,1..26])'
 stcox c.pgs_bmi c.pgs_dpw2 c.pgs_dep1  c.pgs_ea6 c.pgs_eversmoke2 c.pgs_selfhealth1  female pc1-pc10 i.aineistovuosi2 if phsample==1 
mat mo3g=(r(table)[1,1..6]\r(table)[5..6,1..6])'  
mat mo3=mo3f[1..5,1..3]\mo3g[1,1..3]\mo3f[6..10,1..3]\mo3g[2,1..3]\mo3f[11..14,1..3]\mo3g[3,1..3]\mo3f[15..18,1..3]\mo3g[4,1..3]\mo3f[19..21,1..3]\mo3g[5,1..3]\mo3f[22..26,1..3]\mo3g[6,1..3]
	
	putexcel h20=matrix(mo3)

 

  stcox  ib2.bmi_cat   c.pgs_bmi ib2.alkog_cat c.pgs_dpw2  i.depsum c.pgs_dep1 i.koul4 c.pgs_ea6  i.smoke c.pgs_eversmoke2  i.selfhealth   c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
  mat mo4=(r(table)[1,1..32]\r(table)[5..6,1..32])'
 
 putexcel k20=matrix(mo4)
   
   
 *************************************************
 *Schoenfeld residuals (Supplementary Table S4)  *
 *************************************************
 
      stset endage, failure(kdum==1) enter(startage)   origin(time 0)
 
mat scoe=J(35,2,.)  
local j=1
foreach i in pgs_activity1 pgs_adhd1 pgs_adventure1  pgs_afb2 pgs_asteczrhi1 pgs_asthma1 pgs_audit1 pgs_bmi2 pgs_cannabis2 pgs_cp2 pgs_cpd2 pgs_dep1 pgs_dpw2 pgs_ea6 pgs_eversmoke2 pgs_extra2 pgs_famsat1 pgs_friendsat1 pgs_hayfever1 pgs_height3 pgs_highmath1 pgs_leftout1 pgs_menarche1 pgs_migraine1 pgs_morning1 pgs_narcis1 pgs_nearsighted1 pgs_nebwomen2 pgs_neuro2 pgs_open2 pgs_reading1 pgs_religatt1 pgs_risk2 pgs_selfhealth1 pgs_selfmath1 {



qui stcox `i'  female pc1-pc10 i.aineistovuosi2 if  dz==., 
local popb=_b[`i']
local popse = _se[`i']
 estat phtest, detail
 di "`i' est: " exp(`popb') " 95% ci: " exp(`popb' +invnormal(0.025)*`popse') " " exp(`popb' +invnormal(0.975)*`popse')  " z: " abs(`popb'/`popse') " N: " e(N) " Deaths: " e(N_fail) " Risk years: " e(risk)
di "schoenfeld corr: " r(phtest)[1,1]  "; p: " r(phtest)[1,4]
mat scoe[`j',1]=r(phtest)[1,1]
mat scoe[`j',2]=r(phtest)[1,4]
local j=`j'+1
}
   

mat scoeph=J(6,2,.)  
local j=1
foreach i in  bmi_std alkog_std depsum_std koulv_std smoke_std selfhealth_std  {



qui stcox `i'  female pc1-pc10 i.aineistovuosi2 if phsample==1, 
local popb=_b[`i']
local popse = _se[`i']
 estat phtest, detail
 di "`i' est: " exp(`popb') " 95% ci: " exp(`popb' +invnormal(0.025)*`popse') " " exp(`popb' +invnormal(0.975)*`popse')  " z: " abs(`popb'/`popse') " N: " e(N) " Deaths: " e(N_fail) " Risk years: " e(risk)
di "schoenfeld corr: " r(phtest)[1,1]  "; p: " r(phtest)[1,4]
mat scoeph[`j',1]=r(phtest)[1,1]
mat scoeph[`j',2]=r(phtest)[1,4]
local j=`j'+1
}


 putexcel set "mortality_pgi_jan25.xlsx", modify sheet("schoenfeld")
putexcel b4=matrix(scoe)	   b41=matrix(scoeph)	

**************************************
*Correaltions (Supplementary Table S3)*
**************************************

 putexcel set "mortality_pgi_jan25.xlsx", modify sheet("correlations")

cor pgs_activity1 pgs_adhd1 pgs_adventure1  pgs_afb2 pgs_asteczrhi1 pgs_asthma1 pgs_audit1 pgs_bmi2 pgs_cannabis2 pgs_cp2 pgs_cpd2 pgs_dep1 pgs_dpw2 pgs_ea6 pgs_eversmoke2 pgs_extra2 pgs_famsat1 pgs_friendsat1 pgs_hayfever1 pgs_height3 pgs_highmath1 pgs_leftout1 pgs_menarche1 pgs_migraine1 pgs_morning1 pgs_narcis1 pgs_nearsighted1 pgs_nebwomen2 pgs_neuro2 pgs_open2 pgs_reading1 pgs_religatt1 pgs_risk2 pgs_selfhealth1 pgs_selfmath1 if  dz==.
putexcel b4=matrix(r(C))
scalar nob=r(N)
putexcel a2="A: PGIs (N: `=nob' )"

 cor bmi_std alkog_std depsum_std koulv_std smoke_std selfhealth_std  pgs_bmi pgs_dpw2 pgs_dep1  pgs_ea6 pgs_eversmoke2 pgs_selfhealth1  if phsample==1 
putexcel b41=matrix(r(C))

scalar nob=r(N)

putexcel a39="B: Selected PGIs and phenotypes (N: `=nob' )"

**********************************************
*PGIs with/without risk phenotypes (Table 2) *
**********************************************

 putexcel set "mortality_pgi_jan26.xlsx", modify sheet("riskfactor")
	  
      stset endage, failure(kdum==1) enter(startage)   origin(time 0)
	  
	  
 stcox c.pgs_eversmoke2  female pc1-pc10 i.aineistovuosi2 if smoke==1  & dz==.
 mat ests=_b[pgs_eversmoke2],_se[pgs_eversmoke2], e(N),e(N_fail)
 putexcel b3=matrix(ests)
 stcox c.pgs_eversmoke2  female pc1-pc10 i.aineistovuosi2 if smoke!=1 & smoke!=. & dz==.
 mat ests=_b[pgs_eversmoke2],_se[pgs_eversmoke2], e(N),e(N_fail)
 putexcel b4=matrix(ests)
	  
 stcox c.pgs_bmi  female pc1-pc10 i.aineistovuosi2 if bmi_cat==2  & dz==.
  mat ests=_b[pgs_bmi],_se[pgs_bmi], e(N),e(N_fail)
 putexcel b5=matrix(ests)
 stcox c.pgs_bmi  female pc1-pc10 i.aineistovuosi2 if bmi_cat!=2 & bmi_cat!=. & dz==.
  mat ests=_b[pgs_bmi],_se[pgs_bmi], e(N),e(N_fail)
 putexcel b6=matrix(ests)
 
 stcox c.pgs_dep1  female pc1-pc10 i.aineistovuosi2 if depsum==0  & dz==.
  mat ests=_b[pgs_dep1],_se[pgs_dep1], e(N),e(N_fail)
 putexcel b7=matrix(ests)
 stcox c.pgs_dep1  female pc1-pc10 i.aineistovuosi2 if depsum!=0 & depsum!=. & dz==.
 mat ests=_b[pgs_dep1],_se[pgs_dep1], e(N),e(N_fail)
 putexcel b8=matrix(ests)
 
 stcox c.pgs_dpw2  female pc1-pc10 i.aineistovuosi2 if alkog_cat==1  & dz==.
  mat ests=_b[pgs_dpw2],_se[pgs_dpw2], e(N),e(N_fail)
 putexcel b9=matrix(ests)
 stcox c.pgs_dpw2  female pc1-pc10 i.aineistovuosi2 if alkog_cat!=1 & alkog_cat!=. & dz==.
  mat ests=_b[pgs_dpw2],_se[pgs_dpw2], e(N),e(N_fail)
 putexcel b10=matrix(ests)
 
 stcox c.pgs_ea6  female pc1-pc10 i.aineistovuosi2 if koul4==4  & dz==.
  mat ests=_b[pgs_ea6],_se[pgs_ea6], e(N),e(N_fail)
 putexcel b11=matrix(ests)
 stcox c.pgs_ea6  female pc1-pc10 i.aineistovuosi2 if koul4!=4 & koul4!=. & dz==.
  mat ests=_b[pgs_ea6],_se[pgs_ea6], e(N),e(N_fail)
 putexcel b12=matrix(ests)

  stcox c.pgs_selfhealth1  female pc1-pc10 i.aineistovuosi2 if selfhealth_rev==5  & dz==.
   mat ests=_b[pgs_selfhealth1],_se[pgs_selfhealth1], e(N),e(N_fail)
 putexcel b13=matrix(ests)
 stcox c.pgs_selfhealth1  female pc1-pc10 i.aineistovuosi2 if selfhealth_rev!=5 & selfhealth_rev!=. & dz==.
  mat ests=_b[pgs_selfhealth1],_se[pgs_selfhealth1], e(N),e(N_fail)
 putexcel b14=matrix(ests)
 
 
 
 
 
*****************************************
*Descriptives (Supplementary Table S2) *
*****************************************

  putexcel set "mortality_pgi_jan25.xlsx", modify sheet("descriptives")


mean bmi alkog depsum koulv smoke selfhealth_rev   if dz==.
 putexcel b4=matrix(e(b)') c4=matrix(e(sd)')
 
 
 tab bmi_cat   if dz==., matcell(freq)
 mat pct=freq/r(N)*100
  putexcel b13=matrix(pct) d4=(r(N))

   
 tab alkog_cat   if dz==., matcell(freq)
 mat pct=freq/(r(N))*100
  putexcel b19=matrix(pct) d5=(r(N))
  
  
   
 tab depsum   if dz==., matcell(freq)
 mat pct=freq/(r(N))*100
  putexcel b25=matrix(pct) d6=(r(N))
  
   tab koul4   if dz==., matcell(freq)
 mat pct=freq/(r(N))*100
  putexcel b30=matrix(pct) d7=(r(N))
  
     tab smoke   if dz==., matcell(freq)
 mat pct=freq/(r(N))*100
  putexcel b36=matrix(pct) d8=(r(N))
  
    
     tab selfhealth_rev   if dz==., matcell(freq)
 mat pct=freq/(r(N))*100
  putexcel b40=matrix(pct) d9=(r(N))


*************************************************
*Information criteria (Supplmentary table S8 )  *
*************************************************
       stset endage, failure(kdum==1) enter(startage)   origin(time 0)

	   *phenotype-pgi-analysis ics/model tests
	 
  *continuous
  
 * models 1
   qui stcox smoke_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
   est sto smop
estat ic
mat cm1a=r(S)[1,5..6]
  qui stcox c.pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1
  est sto smog
estat ic 
mat cm1b=r(S)[1,5..6]
 
qui stcox bmi_std female pc1-pc10 i.aineistovuosi2 if phsample==1
est sto bmip
estat ic 
mat cm1c=r(S)[1,5..6]
qui stcox c.pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
est sto bmig
estat ic 
mat cm1d=r(S)[1,5..6]




 qui stcox depsum_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto depp
estat ic 
mat cm1e=r(S)[1,5..6]
 qui stcox c.pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
 est sto depg
 estat ic 
mat cm1f=r(S)[1,5..6]

 qui stcox alkog_std female pc1-pc10 i.aineistovuosi2 if phsample==1
 est sto alkp
estat ic 
mat cm1g=r(S)[1,5..6]
  qui stcox c.pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
  est sto alkg
estat ic 
mat cm1h=r(S)[1,5..6]

 qui stcox koulv_std  female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto edup
 estat ic 
mat cm1i=r(S)[1,5..6]

  qui stcox c.pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto edug
  estat ic
mat cm1j=r(S)[1,5..6]



 qui stcox selfhealth_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto srhp
 estat ic 
mat cm1k=r(S)[1,5..6]


 qui stcox c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto srhg
 estat ic
mat cm1l=r(S)[1,5..6]


*models 2

 qui stcox smoke_std c.pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1 
 estat ic 
 mat cm2a=r(S)[1,5..6]

lrtest smop
scalar cm1am2a=r(p)
lrtest smog
scalar cm1bm2a=r(p)


qui stcox bmi_std c.pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
estat ic
 mat cm2b=r(S)[1,5..6]
 
lrtest bmip
scalar cm1cm2b=r(p)
lrtest bmig
scalar cm1dm2b=r(p)


 qui stcox depsum_std c.pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
estat ic 
 mat cm2c=r(S)[1,5..6]

lrtest depp
scalar cm1em2c=r(p)
lrtest depg
scalar cm1fm2c=r(p)


 qui stcox alkog_std c.pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
estat ic 
 mat cm2d=r(S)[1,5..6]

lrtest alkp
scalar cm1gm2d=r(p)
lrtest alkg
scalar cm1hm2d=r(p)



 qui stcox koulv_std  c.pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
estat ic 
 mat cm2e=r(S)[1,5..6]

lrtest edup
scalar cm1im2e=r(p)
lrtest edug
scalar cm1jm2e=r(p)


 qui stcox selfhealth_std c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
estat ic 
 mat cm2f=r(S)[1,5..6]

lrtest srhp	
scalar cm1km2f=r(p)

lrtest srhg
scalar cm1lm2f=r(p)


*model 3

 qui stcox smoke_std bmi_std depsum_std alkog_std  koulv_std   selfhealth_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto fullp
estat ic 
 mat cm3a=r(S)[1,5..6]

qui stcox c.pgs_eversmoke2  c.pgs_bmi  c.pgs_dep1 c.pgs_dpw2  c.pgs_ea6 c.pgs_selfhealth1  female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto fullg
estat ic 
 mat cm3b=r(S)[1,5..6]



*model 4
  qui stcox   smoke_std c.pgs_eversmoke2  bmi_std   c.pgs_bmi depsum_std c.pgs_dep1  alkog_std c.pgs_dpw2   koulv_std  c.pgs_ea6  selfhealth_std   c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
  estat ic 
   mat cm4=r(S)[1,5..6]

lrtest fullp
scalar cm3am4=r(p)
lrtest fullg 
scalar cm3bm4=r(p)

mat cm2p=cm1am2a\cm1bm2a\cm1cm2b\cm1dm2b\cm1em2c\cm1fm2c\cm1gm2d\cm1hm2d\cm1im2e\cm1jm2e\cm1km2f\cm1lm2f

  
 mat conic=cm1a\cm1b\cm1c\cm1d\cm1e\cm1f\cm1g\cm1h\cm1i\cm1j\cm1k\cm1l\cm2a\cm2b\cm2c\cm2d\cm2e\cm2f\cm3a\cm3b\cm4
  
  putexcel set "modeltests.xlsx", modify
  putexcel b4=matrix(conic)
    putexcel e4=matrix(cm2p) e22=cm3am4   e23=cm3bm4

  
  *categorical phenotypes*
  

* models 1
   qui stcox i.smoke female pc1-pc10 i.aineistovuosi2 if phsample==1 
   est sto smop
estat ic
mat dm1a=r(S)[1,5..6]
  qui stcox c.pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1
  est sto smog
estat ic 
mat dm1b=r(S)[1,5..6]
 
qui stcox ib2.bmi_cat  female pc1-pc10 i.aineistovuosi2 if phsample==1
est sto bmip
estat ic 
mat dm1c=r(S)[1,5..6]
qui stcox c.pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
est sto bmig
estat ic 
mat dm1d=r(S)[1,5..6]




 qui stcox i.depsum female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto depp
estat ic 
mat dm1e=r(S)[1,5..6]
 qui stcox c.pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
 est sto depg
 estat ic 
mat dm1f=r(S)[1,5..6]

 qui stcox ib2.alkog_cat female pc1-pc10 i.aineistovuosi2 if phsample==1
 est sto alkp
estat ic 
mat dm1g=r(S)[1,5..6]
  qui stcox c.pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
  est sto alkg
estat ic 
mat dm1h=r(S)[1,5..6]

 qui stcox i.koul4  female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto edup
 estat ic 
mat dm1i=r(S)[1,5..6]

  qui stcox c.pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto edug
  estat ic
mat dm1j=r(S)[1,5..6]



 qui stcox selfhealth_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto srhp
 estat ic 
mat dm1k=r(S)[1,5..6]


 qui stcox c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto srhg
 estat ic
mat dm1l=r(S)[1,5..6]


*models 2

 qui stcox i.smoke c.pgs_eversmoke2 female pc1-pc10 i.aineistovuosi2 if phsample==1 
 estat ic 
 mat dm2a=r(S)[1,5..6]

lrtest smop
scalar dm1am2a=r(p)
lrtest smog
scalar dm1bm2a=r(p)


qui stcox ib2.bmi_cat  c.pgs_bmi female pc1-pc10 i.aineistovuosi2 if phsample==1
estat ic
 mat dm2b=r(S)[1,5..6]
 
lrtest bmip
scalar dm1cm2b=r(p)
lrtest bmig
scalar dm1dm2b=r(p)


 qui stcox i.depsum c.pgs_dep1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
estat ic 
 mat dm2c=r(S)[1,5..6]

lrtest depp
scalar dm1em2c=r(p)
lrtest depg
scalar dm1fm2c=r(p)


 qui stcox ib2.alkog_cat c.pgs_dpw2 female pc1-pc10 i.aineistovuosi2 if phsample==1
estat ic 
 mat dm2d=r(S)[1,5..6]

lrtest alkp
scalar dm1gm2d=r(p)
lrtest alkg
scalar dm1hm2d=r(p)



 qui stcox i.koul4  c.pgs_ea6 female pc1-pc10 i.aineistovuosi2 if phsample==1 
estat ic 
 mat dm2e=r(S)[1,5..6]

lrtest edup
scalar dm1im2e=r(p)
lrtest edug
scalar dm1jm2e=r(p)


 qui stcox selfhealth_std c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
estat ic 
 mat dm2f=r(S)[1,5..6]

lrtest srhp	
scalar dm1km2f=r(p)

lrtest srhg
scalar dm1lm2f=r(p)


*model 3

 qui stcox i.smoke ib2.bmi_cat  i.depsum ib2.alkog_cat  i.koul4   selfhealth_std female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto fullp
estat ic 
 mat dm3a=r(S)[1,5..6]

qui stcox c.pgs_eversmoke2  c.pgs_bmi  c.pgs_dep1 c.pgs_dpw2  c.pgs_ea6 c.pgs_selfhealth1  female pc1-pc10 i.aineistovuosi2 if phsample==1 
est sto fullg
estat ic 
 mat dm3b=r(S)[1,5..6]



*model 4
  qui stcox   i.smoke c.pgs_eversmoke2  ib2.bmi_cat    c.pgs_bmi i.depsum c.pgs_dep1  ib2.alkog_cat c.pgs_dpw2   i.koul4  c.pgs_ea6  selfhealth_std   c.pgs_selfhealth1 female pc1-pc10 i.aineistovuosi2 if phsample==1 
  estat ic 
   mat dm4=r(S)[1,5..6]

lrtest fullp
scalar dm3am4=r(p)
lrtest fullg 
scalar dm3bm4=r(p)

mat dm2p=dm1am2a\dm1bm2a\dm1dm2b\dm1dm2b\dm1em2c\dm1fm2c\dm1gm2d\dm1hm2d\dm1im2e\dm1jm2e\dm1km2f\dm1lm2f
  
 mat disic=dm1a\dm1b\dm1c\dm1d\dm1e\dm1f\dm1g\dm1h\dm1i\dm1j\dm1k\dm1l\dm2a\dm2b\dm2c\dm2d\dm2e\dm2f\dm3a\dm3b\dm4
  
  putexcel set "modeltests.xlsx", modify
  putexcel b27=matrix(disic)
    putexcel e27=matrix(dm2p) e45=dm3am4   e46=dm3bm4


log close

*
exit
*


***
*To finish, some stastics based on aggregate coefficients
***


import excel "mortality_pgi_jan26.xlsx",  cellrange(A3:Y38) firstrow clear sheet("PGI")


rename (loghr logse  D E F G H I J K L M N O P Q R S T U V W X Y) (pop_loghr pop_logse sib_loghr sib_logse men_loghr men_logse women_loghr women_logse basic_loghr basic_logse sec_loghr sec_logse tert_loghr tert_logse age25_loghr age25_logse  age65_loghr age65_logse  age80_loghr age80_logse  extr_loghr extr_logse  nonex_loghr nonex_logse)

lab def labpgi	1 "Physical Activity"
lab def labpgi	2 "ADHD", add
lab def labpgi	3 "Adventurousness", add
lab def labpgi	4 "Age at First Birth", add
lab def labpgi	5 "Asthma/Eczema/Rhinitis", add
lab def labpgi	6 "Asthma", add
lab def labpgi	7 "Alcohol Misuse", add
lab def labpgi	8 "Body Mass Index", add
lab def labpgi	9 "Cannabis Use", add
lab def labpgi	10 "Cognitive Performance", add
lab def labpgi	11 "Cigarettes per Day", add
lab def labpgi	12 "Depressive Symptoms", add
lab def labpgi	13 "Drinks per Week", add
lab def labpgi	14 "Educational Attainment", add
lab def labpgi	15 "Ever Smoker", add
lab def labpgi	16 "Extraversion", add
lab def labpgi	17 "Life Satisfaction - Family", add
lab def labpgi	18 "Life Satisfaction - Friend", add
lab def labpgi	19 "Hayfever (Allergic Rhinitis)", add
lab def labpgi	20 "Height", add
lab def labpgi	21 "Highest Math", add
lab def labpgi	22 "Left out of Social Activity", add
lab def labpgi	23 "Age at Menarche", add
lab def labpgi	24 "Migraine", add
lab def labpgi	25 "Morning Person", add
lab def labpgi	26 "Narcissism", add
lab def labpgi	27 "Nearsightedness", add
lab def labpgi	28 "Number Ever Born (women)", add
lab def labpgi	29 "Neuroticism", add
lab def labpgi	30 "Openness", add
lab def labpgi	31 "Childhood Reading", add
lab def labpgi	32 "Religious Attendance", add
lab def labpgi	33 "Risk Tolerance", add
lab def labpgi	34 "Self-Rated Health ", add
lab def labpgi	35 "Self-Rated Math Ability", add


encode pgi, gen(pgi2) 
lab val pgi2 labpgi


gen b1="("
gen b2=")"
gen b3="; "



foreach i in pop sib men women basic sec tert age25 age65  age80  extr  nonex  {
gen hr_`i'=exp(`i'_loghr)
gen lci_`i'=exp(`i'_loghr+ invnormal(0.025)*`i'_logse)
gen uci_`i'=exp(`i'_loghr+ invnormal(0.975)*`i'_logse)
egen ci_`i'=concat(b1 lci_`i' b3 uci_`i' b2), format(%9.2f)
gen p_`i'=normal(-abs(`i'_loghr/`i'_logse))*2



list pgi2 hr_`i' ci_`i' p_`i'
}



* Benjamini–Hochberg p-values
gen ord_origin=_n
gen ord=.

foreach i in  p_pop p_sib p_men p_women p_basic p_sec p_tert p_age25 p_age65 p_age80 p_extr p_nonex {
sort `i'
replace ord=_n
gen  `i'_BH= `i'/(ord/35)
	forvalues j=34(-1)1 {
	replace `i'_BH=`i'_BH[_n+1] in `j' if `i'_BH>`i'_BH[_n+1]
	}
}
sort ord_origin




gen sibdiff=abs(sib_loghr)-abs(pop_loghr)
gen sibdiffse=sqrt(sib_logse^2+pop_logse^2)


gen gdiff=abs(men_loghr)-abs(women_loghr)
gen gdiffse=sqrt(men_logse^2+women_logse^2)

gen secdiff=abs(sec_loghr)-abs(basic_loghr)
gen secdiffse=sqrt(sec_logse^2+basic_logse^2)

gen tertdiff=abs(tert_loghr)-abs(basic_loghr)
gen tertdiffse=sqrt(tert_logse^2+basic_logse^2)


gen age65diff=abs(age65_loghr)-abs(age25_loghr)
gen age65diffse=sqrt(age65_logse^2+age25_logse^2)

gen age80diff=abs(age80_loghr)-abs(age25_loghr)
gen age80diffse=sqrt(age80_logse^2+age25_logse^2)

gen extrdiff=abs(extr_loghr)-abs(nonex_loghr)
gen extrdiffse=sqrt(extr_logse^2+nonex_logse^2)



meta set sibdiff sibdiffse, fixed
meta summarize

mat medif=r(theta),r(ci_lb), r(ci_ub),r(p)

meta set gdiff gdiffse, fixed
meta summarize

local gdadhd=gdiff in 2
local gdadhd_se=gdiffse in 2
di "ADHD-sex-difference p: " normal(-abs(`gdadhd'/`gdadhd_se'))*2
di "ADHD-sex-difference p, BH-Adjusted:  " normal(-abs(`gdadhd'/`gdadhd_se'))*2/(2/35)
*edit 9.2.2026: fixed p-value order (1/35->2/35)

local gdea=gdiff in 14 
local gdea_se=gdiffse in 14
di "education-sex-difference p: " normal(-abs(`gdea'/`gdea_se'))*2
di "education-sex-difference p:, BH-Adjusted: "  normal(-abs(`gdea'/`gdea_se'))*2/(1/35)
*edit 9.2.2026: larger p than for ADHD, so we use that


mat medif=medif\(r(theta),r(ci_lb), r(ci_ub),r(p))

meta set secdiff secdiffse, fixed
meta summarize

mat medif=medif\(r(theta),r(ci_lb), r(ci_ub),r(p))

meta set tertdiff tertdiffse, fixed
meta summarize

mat medif=medif\(r(theta),r(ci_lb), r(ci_ub),r(p))

meta set age65diff age65diffse, fixed
meta summarize

local adadhd=age65diff in 2
local adadhd_se=age65diffse in 2
di "ADHD-age-difference 25-64 vs 65-79 p: " normal(-abs(`adadhd'/`adadhd_se'))*2 ", BH-Adjusted: "  normal(-abs(`adadhd'/`adadhd_se'))*2/(1/35)


gen age80diffp=normal(-abs(age80diff/age80diffse))*2
sort age80diffp
replace ord=_n
gen  age80diffp_BH= age80diffp/(ord/35)
	forvalues j=34(-1)1 {
	replace age80diffp_BH=age80diffp_BH[_n+1] in `j' if age80diffp_BH>age80diffp_BH[_n+1]
	}
sort pgi
list pgi age80diff age80diffse age80diffp age80diffp_BH
*edit 9.2.2026: added calculations for age80diffp_BH

mat medif=medif\(r(theta),r(ci_lb), r(ci_ub),r(p))

meta set extrdiff extrdiffse, fixed
meta summarize
mat medif=medif\(r(theta),r(ci_lb), r(ci_ub),r(p))

*as some estimates have opposte signs, let's calculate once more (without absolute value in the first stage)
gen extrdiff2=extr_loghr-nonex_loghr
gen extrdiffp=normal(-abs((extrdiff2)/extrdiffse))*2
list pgi  hr_extr hr_nonex  extr_loghr nonex_loghr extrdiff  extrdiff2 extrdiffse extrdiffp

