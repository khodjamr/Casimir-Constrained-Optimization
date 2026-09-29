(* ::Package:: *)

$HistoryLength=0;

(*----------output directory----------*)
baseDirectory=If[StringQ[$InputFileName]&&StringLength[$InputFileName]>0,DirectoryName[$InputFileName],If[$FrontEnd=!=Null,NotebookDirectory[],Directory[]]];
dataDirectory=FileNameJoin[{baseDirectory,"data_mathematica9_final"}];
figureDirectory=FileNameJoin[{baseDirectory,"figures_mathematica9_final"}];
If[!DirectoryQ[dataDirectory],CreateDirectory[dataDirectory]];
If[!DirectoryQ[figureDirectory],CreateDirectory[figureDirectory]];

(*----------constants----------*)
hbar=1.054571817*10^-34;
c0=299792458.;
kB=1.380649*10^-23;
eVJ=1.602176634*10^-19;
eps0vac=8.8541878128*10^-12;
T0=300.;
EgGaAs=1.424;(*eV*)etaPhysicalMax=0.950;
etaNominalMax=0.945;
deltaA=10.*10^-9;
deltaEta=0.005;
aMin=0.25*10^-6;
aMax=1.50*10^-6;
cavityContrastMin=0.10*10^-3;(*Pa*)(*GaAs dielectric limits used for the imaginary-axis completion*)epsStatic=12.9;
epsHigh=10.89;
Ephon=0.035;
Euv=10.0;

(*----------Adachi 1989 GaAs optical model----------*)
E0=1.42;
Delta0=1.77-E0;
E1=2.90;
Delta1=3.13-E1;
E2=4.7;
EgD=1.73;
AA=3.45;
B1=6.37;
B11=13.08;
Gamma1=0.10;
CC=2.39;
gammaC=0.146;
DD=24.2;
epsInfAdachi=1.6;

heav[x_?NumericQ]:=If[x>=0.,1.,0.];

epsilonAdachi[E_?NumericQ]:=Module[{chi0,chiso,f0,fso,eps1A,eps2A,chi1r,rawB,eps2B,chi1,eps1Bcomplex,epsB,chi2,denC,epsC,eps2D},chi0=(E+0. I)/E0;
chiso=(E+0. I)/(E0+Delta0);
f0=chi0^-2 (2.-Sqrt[1.+chi0]-Sqrt[(1.-E/E0) heav[1.-E/E0]+0. I]);
fso=chiso^-2 (2.-Sqrt[1.+chiso]-Sqrt[(1.-E/(E0+Delta0)) heav[1.-E/(E0+Delta0)]+0. I]);
eps1A=AA E0^-1.5 (f0+0.5 (E0/(E0+Delta0))^1.5 fso);
eps2A=AA/E^2 (Sqrt[Max[E-E0,0.]] heav[E/E0-1.]+0.5 Sqrt[Max[E-E0-Delta0,0.]] heav[E/(E0+Delta0)-1.]);
chi1r=E/E1;
rawB=If[chi1r==0.,0.,Pi chi1r^-2 (B1-B11 Sqrt[Max[E1-E,0.]] heav[1.-chi1r])];
eps2B=Max[rawB,0.];
chi1=(E+I Gamma1)/E1;
eps1Bcomplex=-B1 chi1^-2 Log[1.-chi1^2];
epsB=Re[eps1Bcomplex]+I eps2B;
chi2=(E+0. I)/E2;
denC=(1.-chi2^2)^2+(chi2 gammaC)^2;
epsC=CC (1.-chi2^2)/denC+I CC chi2 gammaC/denC;
eps2D=DD/E^2 (E-EgD)^2 heav[1.-EgD/E] heav[1.-E/E1];
N[eps1A+I eps2A+epsB+epsC+I eps2D+epsInfAdachi]];

(*----------Gauss-Legendre rule,implemented without post-v9 functions----------*)
ClearAll[glStandard,glRule];
glStandard[n_Integer]:=glStandard[n]=Module[{beta,mat,vals,vecs,ord,nodes,weights},beta=Table[N[j/Sqrt[4. j^2-1.]],{j,1,n-1}];
mat=DiagonalMatrix[ConstantArray[0.,n]]+DiagonalMatrix[beta,1]+DiagonalMatrix[beta,-1];
{vals,vecs}=Eigensystem[mat];
ord=Ordering[vals];
nodes=vals[[ord]];
vecs=vecs[[ord]];
weights=2. (vecs[[All,1]]^2);
{nodes,weights}];

glRule[n_Integer,a_?NumericQ,b_?NumericQ]:=Module[{x,w},{x,w}=glStandard[n];
Transpose[{0.5 (b-a) x+0.5 (a+b),0.5 (b-a) w}]];

(*----------causal imaginary-axis dielectric response----------*)
kkRule=glRule[320,0.001,6.0];
kkE=kkRule[[All,1]];
kkW=kkRule[[All,2]];
kkEps2=Map[Max[Im[epsilonAdachi[#]],0.]&,kkE];
eps0Interband=1.+(2./Pi) Total[kkW kkEps2/kkE];
deltaUV=epsHigh-eps0Interband;
deltaPh=epsStatic-epsHigh;

epsilonImag[xi_?NumericQ]:=Module[{X,inter},X=hbar xi/eVJ;
inter=1.+(2./Pi) Total[kkW kkE kkEps2/(kkE^2+X^2)];
N[inter+deltaUV Euv^2/(Euv^2+X^2)+deltaPh Ephon^2/(Ephon^2+X^2)]];

(*----------equilibrium Lifshitz pressure:GaAs half-space vs PEC----------*)
peq[a_?NumericQ,nU_Integer:180]:=Module[{sum=0.,n=0,xi,u0,upper,rule,eps,term,u,us,rte,rtm,qte,qtm,integrand},While[True,xi=2. Pi n kB T0/hbar;
u0=xi a/c0;
If[n>0&&u0>30.,Break[]];
upper=Max[35.,u0+8.];
rule=glRule[nU,u0,upper];
eps=If[n==0,epsStatic,epsilonImag[xi]];
term=Total[Map[Function[pair,u=pair[[1]];
us=Sqrt[u^2+(eps-1.) u0^2];
rte=(u-us)/(u+us);
rtm=(eps u-us)/(eps u+us);
qte=(-rte) Exp[-2. u];
qtm=rtm Exp[-2. u];
integrand=u^2 (qte/(1.-qte)+qtm/(1.-qtm));
pair[[2]] integrand],rule]];
If[n==0,term=0.5 term];
sum=sum+term;
n++;
If[n>5000,Print["ERROR: Matsubara sum failed to terminate."];Abort[]];];
N[-(kB T0)/(Pi a^3) sum]];

(*----------real-frequency Fresnel amplitudes----------*)
fresnelReal[E_?NumericQ,q_?NumericQ]:=Module[{eps,kz0,kzs,rte,rtm},eps=epsilonAdachi[E];
kz0=Sqrt[1.-q^2+0. I];
kzs=Sqrt[eps-q^2+0. I];
If[Im[kzs]<0.,kzs=-kzs];
rte=(kz0-kzs)/(kz0+kzs);
rtm=(eps kz0-kzs)/(eps kz0+kzs);
{rte,rtm,kzs}];

occupation[E_?NumericQ,eta_?NumericQ,eg_?NumericQ]:=1./(Exp[(E-eta eg) eVJ/(kB T0)]-1.);
occupation[E_?NumericQ,eta_?NumericQ]:=occupation[E,eta,EgGaAs];

excessOccupation[E_?NumericQ,eta_?NumericQ,eg_?NumericQ]:=occupation[E,eta,eg]-occupation[E,0.,eg];
excessOccupation[E_?NumericQ,eta_?NumericQ]:=excessOccupation[E,eta,EgGaAs];

(*Rows returned:{E[eV],weight[eV],Kpw[Pa/eV/occupation],Kew[...]}.*)
ClearAll[noneqKernels];
noneqKernels[a_?NumericQ,nE_Integer,nPW_Integer,nEW_Integer,eWindow_?NumericQ,alphaMax_?NumericQ]:=Module[{eRule,tRule,alphaRule,rows,E,wE,k0,ip,ie,t,wt,q,fres,rte,rtm,phase,denTE,denTM,absTE,absTM,alpha,wa,phaseE,kPW,kEW},eRule=glRule[nE,EgGaAs,EgGaAs+eWindow];
tRule=glRule[nPW,0.,1.];
alphaRule=glRule[nEW,0.,alphaMax];
rows=Map[Function[epair,E=epair[[1]];wE=epair[[2]];
k0=E eVJ/(hbar c0);
ip=Total[Map[Function[tpair,t=tpair[[1]];wt=tpair[[2]];
q=Sqrt[Max[1.-t^2,0.]];
fres=fresnelReal[E,q];rte=fres[[1]];rtm=fres[[2]];
phase=Exp[2. I k0 a t];
denTE=Abs[1.-rte (-1.) phase]^2;
denTM=Abs[1.-rtm phase]^2;
absTE=Max[1.-Abs[rte]^2,0.];
absTM=Max[1.-Abs[rtm]^2,0.];
wt t^2 2. (absTE/denTE+absTM/denTM)],tRule]];
kPW=eVJ/(4. Pi^2) k0^3 ip;
ie=Total[Map[Function[apair,alpha=apair[[1]];wa=apair[[2]];
q=Sqrt[1.+alpha^2];
fres=fresnelReal[E,q];rte=fres[[1]];rtm=fres[[2]];
phaseE=Exp[-2. k0 a alpha];
denTE=Abs[1.-rte (-1.) phaseE]^2;
denTM=Abs[1.-rtm phaseE]^2;
wa alpha^2 phaseE ((-Im[rte])/denTE+Im[rtm]/denTM)],alphaRule]];
kEW=-eVJ/Pi^2 k0^3 ie;
{E,wE,Re[kPW],Re[kEW]}],eRule];
rows];
noneqKernels[a_?NumericQ,nE_Integer,nPW_Integer,nEW_Integer,eWindow_?NumericQ]:=noneqKernels[a,nE,nPW,nEW,eWindow,10.];
noneqKernels[a_?NumericQ,nE_Integer,nPW_Integer,nEW_Integer]:=noneqKernels[a,nE,nPW,nEW,0.60,10.];
noneqKernels[a_?NumericQ,nE_Integer,nPW_Integer]:=noneqKernels[a,nE,nPW,180,0.60,10.];
noneqKernels[a_?NumericQ,nE_Integer]:=noneqKernels[a,nE,200,180,0.60,10.];
noneqKernels[a_?NumericQ]:=noneqKernels[a,180,200,180,0.60,10.];

noneqFromKernels[kernels_List,eta_?NumericQ]:=Module[{pw,ew},pw=Total[Map[(#[[2]] #[[3]] excessOccupation[#[[1]],eta])&,kernels]];
ew=Total[Map[(#[[2]] #[[4]] excessOccupation[#[[1]],eta])&,kernels]];
{N[pw],N[ew]}];

pressureComponents[a_?NumericQ,eta_?NumericQ,nU_Integer:180,nE_Integer:180,nPW_Integer:200,nEW_Integer:180]:=Module[{peqv,kernels,ne},peqv=peq[a,nU];
If[eta<=0.,Return[{peqv,0.,0.,peqv}]];
kernels=noneqKernels[a,nE,nPW,nEW];
ne=noneqFromKernels[kernels,eta];
{peqv,ne[[1]],ne[[2]],peqv+ne[[1]]+ne[[2]]}];

(*Used only for energy-window and evanescent-cutoff convergence checks.*)
pressureComponentsCutoffs[a_?NumericQ,eta_?NumericQ,nU_Integer,nE_Integer,nPW_Integer,nEW_Integer,eWindow_?NumericQ,alphaMax_?NumericQ]:=Module[{peqv,kernels,ne},peqv=peq[a,nU];
If[eta<=0.,Return[{peqv,0.,0.,peqv}]];
kernels=noneqKernels[a,nE,nPW,nEW,eWindow,alphaMax];
ne=noneqFromKernels[kernels,eta];
{peqv,ne[[1]],ne[[2]],peqv+ne[[1]]+ne[[2]]}];

ClearAll[finiteSlabRT,finiteCavityU,finiteSlabEqSide,
  finiteSlabEquilibrium,finiteSlabNoneq,
  finiteSlabPressureComponents,RunFiniteSlabValidation];

finiteSlabRT[r_?NumberQ,kzs_?NumberQ,k0_?NumericQ,b_?NumericQ]:=
 Module[{z,den},
  z=Exp[I k0 kzs b];
  den=1.-r^2 z^2;
  {r (1.-z^2)/den,(1.-r^2) z/den}
 ];

finiteCavityU[rho_?NumberQ,tau_?NumberQ,rL_?NumberQ,rR_?NumberQ]:=
 Module[{det},
  det=(1.-rho rL) (1.-rho rR)-tau^2 rL rR;
  {{1.-rho rR,tau rR},{tau rL,1.-rho rL}}/det
 ];

finiteSlabEqSide[aSide_?NumericQ,aOther_?NumericQ,b_?NumericQ,
  nU_Integer]:=
 Module[{sum=0.,n=0,xi,u0,upper,rule,eps,term,u,us,z,
   rList,rcList={-1.,1.},p,r,rc,den,rho,tau,rOther,rEff,q,
   integrand},
  While[True,
   xi=2. Pi n kB T0/hbar;
   u0=xi aSide/c0;
   If[n>0&&u0>30.,Break[]];
   upper=Max[35.,u0+8.];
   rule=glRule[nU,u0,upper];
   eps=If[n==0,epsStatic,epsilonImag[xi]];
   term=Total[Map[Function[pair,
      u=pair[[1]];
      us=Sqrt[u^2+(eps-1.) u0^2];
      z=Exp[-us b/aSide];
      rList={(u-us)/(u+us),(eps u-us)/(eps u+us)};
      integrand=0.;
      Do[
       r=rList[[p]];
       rc=rcList[[p]];
       den=1.-r^2 z^2;
       rho=r (1.-z^2)/den;
       tau=(1.-r^2) z/den;
       rOther=rc Exp[-2. u aOther/aSide];
       rEff=rho+tau^2 rOther/(1.-rho rOther);
       q=rc rEff Exp[-2. u];
       integrand=integrand+q/(1.-q),
       {p,1,2}];
      pair[[2]] u^2 integrand],rule]];
   If[n==0,term=0.5 term];
   sum=sum+term;
   n++;
   If[n>5000,
    Print["ERROR: finite-slab Matsubara sum failed to terminate."];
    Abort[]];
  ];
  N[-(kB T0)/(Pi aSide^3) sum]
 ];

finiteSlabEqSide[aSide_?NumericQ,aOther_?NumericQ,b_?NumericQ]:=
 finiteSlabEqSide[aSide,aOther,b,300];

finiteSlabEquilibrium[aL_?NumericQ,aR_?NumericQ,b_?NumericQ,
  nU_Integer]:=
 {finiteSlabEqSide[aL,aR,b,nU],
  finiteSlabEqSide[aR,aL,b,nU]};

finiteSlabEquilibrium[aL_?NumericQ,aR_?NumericQ,b_?NumericQ]:=
 finiteSlabEquilibrium[aL,aR,b,300];

finiteSlabNoneq[aL_?NumericQ,aR_?NumericQ,b_?NumericQ,
  eta_?NumericQ,nE_Integer,nPW_Integer,nEW_Integer,
  eWindow_?NumericQ,alphaMax_?NumericQ]:=
 Module[{eRule,tRule,alphaRule,pw={0.,0.},ew={0.,0.},E,wE,k0,dn,
   ip,ie,t,wt,q,fres,rList,kzs,p,r,rc,rt,rho,tau,phaseVec,
   rL,rR,uMat,src,cov,absorp,cross,alpha,wa},
  eRule=glRule[nE,EgGaAs,EgGaAs+eWindow];
  tRule=glRule[nPW,0.,1.];
  alphaRule=glRule[nEW,0.,alphaMax];

  Do[
   E=eRule[[i,1]];
   wE=eRule[[i,2]];
   k0=E eVJ/(hbar c0);
   dn=excessOccupation[E,eta];

   ip={0.,0.};
   Do[
    t=tRule[[j,1]];
    wt=tRule[[j,2]];
    q=Sqrt[Max[1.-t^2,0.]];
    fres=fresnelReal[E,q];
    rList={fres[[1]],fres[[2]]};
    kzs=fres[[3]];

    Do[
     r=rList[[p]];
     rc=If[p==1,-1.,1.];
     rt=finiteSlabRT[r,kzs,k0,b];
     rho=rt[[1]];
     tau=rt[[2]];
     rL=rc Exp[2. I k0 t aL];
     rR=rc Exp[2. I k0 t aR];
     uMat=finiteCavityU[rho,tau,rL,rR];
     absorp=Max[Re[1.-Abs[rho]^2-Abs[tau]^2],0.];
     cross=-(rho Conjugate[tau]+tau Conjugate[rho]);
     src={{absorp,cross},{Conjugate[cross],absorp}};
     cov=uMat.src.ConjugateTranspose[uMat];
     ip=ip+wt t^2 2. Re[Diagonal[cov]],
     {p,1,2}],
    {j,1,Length[tRule]}];

   ie={0.,0.};
   Do[
    alpha=alphaRule[[j,1]];
    wa=alphaRule[[j,2]];
    q=Sqrt[1.+alpha^2];
    fres=fresnelReal[E,q];
    rList={fres[[1]],fres[[2]]};
    kzs=fres[[3]];
    phaseVec={Exp[-2. k0 alpha aL],Exp[-2. k0 alpha aR]};

    Do[
     r=rList[[p]];
     rc=If[p==1,-1.,1.];
     rt=finiteSlabRT[r,kzs,k0,b];
     rho=rt[[1]];
     tau=rt[[2]];
     rL=rc phaseVec[[1]];
     rR=rc phaseVec[[2]];
     uMat=finiteCavityU[rho,tau,rL,rR];
     src={{Im[rho],Im[tau]},{Im[tau],Im[rho]}};
     cov=uMat.src.ConjugateTranspose[uMat];
     ie=ie+wa alpha^2 rc phaseVec Re[Diagonal[cov]],
     {p,1,2}],
    {j,1,Length[alphaRule]}];

   pw=pw+wE dn eVJ/(4. Pi^2) k0^3 ip;
   ew=ew-wE dn eVJ/Pi^2 k0^3 ie,
   {i,1,Length[eRule]}];

  {N[pw],N[ew]}
 ];

finiteSlabNoneq[aL_?NumericQ,aR_?NumericQ,b_?NumericQ,
  eta_?NumericQ,nE_Integer,nPW_Integer,nEW_Integer]:=
 finiteSlabNoneq[aL,aR,b,eta,nE,nPW,nEW,0.60,10.];

finiteSlabNoneq[aL_?NumericQ,aR_?NumericQ,b_?NumericQ,
  eta_?NumericQ]:=
 finiteSlabNoneq[aL,aR,b,eta,260,420,320,0.60,10.];

finiteSlabPressureComponents[aL_?NumericQ,aR_?NumericQ,b_?NumericQ,
  eta_?NumericQ,nU_Integer,nE_Integer,nPW_Integer,nEW_Integer]:=
 Module[{eq,ne},
  eq=finiteSlabEquilibrium[aL,aR,b,nU];
  ne=If[eta<=0.,{{0.,0.},{0.,0.}},
    finiteSlabNoneq[aL,aR,b,eta,nE,nPW,nEW]];
  Table[{eq[[i]],ne[[1,i]],ne[[2,i]],
    eq[[i]]+ne[[1,i]]+ne[[2,i]]},{i,1,2}]
 ];

finiteSlabPressureComponents[aL_?NumericQ,aR_?NumericQ,b_?NumericQ,
  eta_?NumericQ]:=
 finiteSlabPressureComponents[aL,aR,b,eta,300,260,420,320];

RunFiniteSlabValidation[]:=
 Module[{aStar=1.3512217912113391*10^-6,bStar=13.5*10^-6,
   calc,nominalRow,activeSpecs,activeRows,activeRow,rows,aL,aR,eta,
   reducedL,reducedR,reduced,full,minFull,delta,file},

  calc[label_String,aL_?NumericQ,aR_?NumericQ,eta_?NumericQ]:=
   Module[{},
    reducedL=pressureComponents[aL,eta,300,260,420,320][[4]];
    reducedR=pressureComponents[aR,eta,300,260,420,320][[4]];
    reduced=Min[reducedL,reducedR];
    full=finiteSlabPressureComponents[aL,aR,bStar,eta,
      300,260,420,320];
    minFull=Min[full[[All,4]]];
    delta=minFull-reduced;
    {label,10.^6 aL,10.^6 aR,eta,1000. reduced,
     1000. full[[1,4]],1000. full[[2,4]],1000. minFull,
     1000. delta,100. Abs[delta/reduced]}
   ];

  nominalRow=calc["nominal",aStar,aStar,0.945];

  activeSpecs={{aStar-deltaA,aStar-deltaA},
    {aStar-deltaA,aStar+deltaA},
    {aStar+deltaA,aStar+deltaA}};

  activeRows=Table[
    calc["active worst case",activeSpecs[[i,1]],
     activeSpecs[[i,2]],0.945-deltaEta],
    {i,1,Length[activeSpecs]}];

  activeRow=First[Sort[activeRows,#1[[8]]<#2[[8]]&]];
  rows={nominalRow,activeRow};

  file=exportCSV["finite_slab_validation.csv",
    {"point","aL_um","aR_um","eta","P_reduced_mPa",
     "P_finite_L_mPa","P_finite_R_mPa","P_finite_min_mPa",
     "Delta_min_mPa","Abs_relative_difference_percent"},rows];

  Print[Grid[Prepend[rows,
    {"point","aL_um","aR_um","eta","P_reduced_mPa",
     "P_finite_L_mPa","P_finite_R_mPa","P_finite_min_mPa",
     "Delta_min_mPa","Abs_rel_percent"}],Frame->All]];

  Print["Finite-slab validation written to ",file];
  rows
 ];


(*----------analytical far-field plateau/material screen----------*)
farFieldPressure[eg_?NumericQ,eta_?NumericQ,nE_Integer:360]:=Module[{rule,total,E,w,k0,dn},rule=glRule[nE,eg,eg+1.50];
total=Total[Map[Function[pair,E=pair[[1]];w=pair[[2]];
k0=E eVJ/(hbar c0);
dn=excessOccupation[E,eta,eg];
w k0^3 dn],rule]];
N[eVJ/(3. Pi^2) total]];

(*Continuous Band-Gap Optimization*)

(* The legacy name bandGapScreen is retained so RunStudy[] needs no other change.
   A coarse scan identifies the global basin, followed by a continuous
   golden-section refinement. *)
ClearAll[bandGapScreen];
bandGapScreen[eta_?NumericQ]:=
 Module[{grid,vals,k,seed,lo,hi,phi,x1,x2,f1,f2,i,
   egBest,pBest},

  grid=Range[0.08,4.00,0.02];
  vals=Map[farFieldPressure[#,eta,240]&,grid];
  k=First[Ordering[vals,-1]];
  seed=grid[[k]];

  lo=Max[0.08,seed-0.08];
  hi=Min[4.00,seed+0.08];
  phi=(Sqrt[5.]-1.)/2.;

  x1=hi-phi (hi-lo);
  x2=lo+phi (hi-lo);
  f1=farFieldPressure[x1,eta,420];
  f2=farFieldPressure[x2,eta,420];

  For[i=1,i<=60,i++,
   If[f1<f2,
    lo=x1;
    x1=x2;
    f1=f2;
    x2=lo+phi (hi-lo);
    f2=farFieldPressure[x2,eta,420],
    hi=x2;
    x2=x1;
    f2=f1;
    x1=hi-phi (hi-lo);
    f1=farFieldPressure[x1,eta,420]
   ];
  ];

  egBest=0.5 (lo+hi);
  pBest=farFieldPressure[egBest,eta,520];
  {N[egBest],N[pBest]}
 ];

RegenerateBandGapOptima[]:=
 Module[{screenEtas,screenRows,eta,egBest,pBest,file},

  screenEtas={0.94,0.945,0.95};

  screenRows=Table[
    {egBest,pBest}=bandGapScreen[eta];
    {eta,egBest,1000. pBest},
    {eta,screenEtas}];

  file=exportCSV["bandgap_optima.csv",
    {"eta","EgOpt_eV","PinfOpt_mPa"},screenRows];

  Print[Grid[Prepend[screenRows,
    {"eta","EgOpt_eV","PinfOpt_mPa"}],Frame->All]];

  Print["Continuous band-gap optima written to ",file];
  screenRows
 ];

(*----------finite-thickness reduction indicators----------*)
opacityLeakage[b_?NumericQ,eta_?NumericQ,nE_Integer,nQ_Integer]:=Module[{eRule,qRule,den=0.,num=0.,E,wE,q,wq,k0,dn,fres,rte,rtm,kzs,wang,atten,base,leak},eRule=glRule[nE,EgGaAs,EgGaAs+0.60];
qRule=glRule[nQ,0.,1.];
Do[E=eRule[[i,1]];wE=eRule[[i,2]];
k0=E eVJ/(hbar c0);dn=excessOccupation[E,eta];
base=0.;leak=0.;
Do[q=qRule[[j,1]];wq=qRule[[j,2]];
fres=fresnelReal[E,q];rte=fres[[1]];rtm=fres[[2]];kzs=fres[[3]];
wang=q Sqrt[1.-q^2] (Max[1.-Abs[rte]^2,0.]+Max[1.-Abs[rtm]^2,0.]);
atten=Exp[-2. k0 Im[kzs] b];
base=base+wq wang;
leak=leak+wq wang atten;,{j,1,Length[qRule]}];
den=den+wE dn k0^3 base;
num=num+wE dn k0^3 leak;,{i,1,Length[eRule]}];
N[num/den]];
opacityLeakage[b_?NumericQ,eta_?NumericQ,nE_Integer]:=opacityLeakage[b,eta,nE,100];
opacityLeakage[b_?NumericQ,eta_?NumericQ]:=opacityLeakage[b,eta,160,100];
opacityLeakage[b_?NumericQ]:=opacityLeakage[b,etaPhysicalMax,160,100];

equilibriumCouplingIndicator[b_?NumericQ,a_?NumericQ]:=(1.+b/a)^-3;
equilibriumCouplingIndicator[b_?NumericQ]:=equilibriumCouplingIndicator[b,aMax];

findEmissionThickness[target_?NumericQ]:=Module[{lo=0.,hi=30.*10^-6,mid,i},For[i=1,i<=45,i++,mid=0.5 (lo+hi);
If[opacityLeakage[mid,etaPhysicalMax,120,80]>target,lo=mid,hi=mid];];
hi];
findEmissionThickness[]:=findEmissionThickness[10^-3];

adoptedThickness[]:=Max[findEmissionThickness[10^-3],9. aMax];

(*----------robust objective----------*)
ClearAll[pressureKernelCache,peqCache];
pressureKernelCache[a_?NumericQ,nE_Integer,nPW_Integer,nEW_Integer]:=pressureKernelCache[a,nE,nPW,nEW]=noneqKernels[a,nE,nPW,nEW];
peqCache[a_?NumericQ,nU_Integer]:=peqCache[a,nU]=peq[a,nU];

pressureCached[a_?NumericQ,eta_?NumericQ,nU_Integer:120,nE_Integer:120,nPW_Integer:160,nEW_Integer:120]:=Module[{peqv,ne,kk},peqv=peqCache[a,nU];
If[eta<=0.,Return[peqv]];
kk=pressureKernelCache[a,nE,nPW,nEW];
ne=noneqFromKernels[kk,eta];
N[peqv+ne[[1]]+ne[[2]]]];

robustOneGap[a0_?NumericQ,eta0_?NumericQ,nU_Integer:120,nE_Integer:120,nPW_Integer:160,nEW_Integer:120]:=Module[{vals},vals=Flatten[Table[pressureCached[a0+sa deltaA,eta0+se deltaEta,nU,nE,nPW,nEW],{sa,{-1.,1.}},{se,{-1.,1.}}]];
Min[vals]];

robustCornerRows[a0_?NumericQ,eta0_?NumericQ,nU_Integer,nE_Integer,nPW_Integer,nEW_Integer]:=Flatten[Table[{sa deltaA,se deltaEta,pressureCached[a0+sa deltaA,eta0+se deltaEta,nU,nE,nPW,nEW]},{sa,{-1.,1.}},{se,{-1.,1.}}],1];

uncertaintyGridRows[a0_?NumericQ,eta0_?NumericQ,nA_Integer,nEta_Integer,nU_Integer,nE_Integer,nPW_Integer,nEW_Integer]:=Module[{daVals,deVals},
If[nA<2||nEta<2,Return[$Failed]];
daVals=Table[-deltaA+2. deltaA i/(nA-1),{i,0,nA-1}];
deVals=Table[-deltaEta+2. deltaEta j/(nEta-1),{j,0,nEta-1}];
Flatten[Table[{da,de,pressureCached[a0+da,eta0+de,nU,nE,nPW,nEW]},{da,daVals},{de,deVals}],1]];

cavityContrastRobust[a0_?NumericQ,eta0_?NumericQ,nU_Integer:120,nE_Integer:120,nPW_Integer:160,nEW_Integer:120]:=robustOneGap[a0,eta0,nU,nE,nPW,nEW]-farFieldPressure[EgGaAs,eta0-deltaEta,300];

feasibleQ[a0_?NumericQ,eta0_?NumericQ,b_?NumericQ]:=And[aMin<=a0<=aMax,0.85<=eta0<=etaNominalMax,EgGaAs (1.-(eta0+deltaEta))>=2.5 kB T0/eVJ,opacityLeakage[b,etaPhysicalMax,100,64]<=10^-3,equilibriumCouplingIndicator[b,aMax]<=10^-3];

(*Coarse-to-fine search.The symmetry theorem in the manuscript reduces the bilateral problem exactly to this one-gap robust search.*)
coarseSearch[]:=Module[{aGrid,etaGrid,rows,a,eta,j,cc,b},b=adoptedThickness[];
(*The endpoints are constructed algebraically so the exported map covers the complete design box exactly.*)
aGrid=Table[aMin+(aMax-aMin) i/50.,{i,0,50}];
etaGrid=Table[0.85+(etaNominalMax-0.85) j/19.,{j,0,19}];
rows=Flatten[Table[j=robustOneGap[a,eta,90,80,100,80];
cc=j-farFieldPressure[EgGaAs,eta-deltaEta,220];
{a,eta,j,cc,If[cc>=cavityContrastMin,1.,0.]},{eta,etaGrid},{a,aGrid}],1];
rows];

refineBestGap[aSeed_?NumericQ,eta0_?NumericQ]:=Module[{f,sol,a,lo,hi},(*At the selected lobe,max-min robustness is obtained by balancing the two gap endpoints. The seed now comes from the full-box coarse search.*)
lo=Max[aMin,aSeed-0.03*10^-6];
hi=Min[aMax,aSeed+0.03*10^-6];
f[x_?NumericQ]:=pressureCached[x-deltaA,eta0-deltaEta,220,300,400,300]-pressureCached[x+deltaA,eta0-deltaEta,220,300,400,300];
sol=Quiet[Check[FindRoot[f[a]==0.,{a,aSeed,lo,hi},AccuracyGoal->8,PrecisionGoal->8,MaxIterations->30],$Failed]];
If[sol===$Failed,Print["FATAL: endpoint-balancing refinement failed in the best coarse lobe."];Return[$Failed]];
a/.sol];
refineBestGap[eta0_?NumericQ]:=refineBestGap[1.351*10^-6,eta0];
refineBestGap[]:=refineBestGap[1.351*10^-6,etaNominalMax];

(*----------exports----------*)
exportCSV[name_String,header_List,rows_List]:=Module[{file},file=FileNameJoin[{dataDirectory,name}];
Export[file,Prepend[rows,header],"CSV"];
file];

numericVectorQ[v_,n_Integer]:=ListQ[v]&&Length[v]==n&&And@@Map[NumberQ,v];
numericTableQ[t_,n_Integer]:=ListQ[t]&&Length[t]>0&&And@@Map[Function[row,numericVectorQ[row,n]],t];

showAndExportPDF[name_String,g_]:=Module[{file},
Print[g];
file=FileNameJoin[{figureDirectory,name}];
Export[file,g,"PDF"];
file];

RunStudy[]:=Module[{screenEtas,screenRows,candidates,materialRows,bEmission,bOpt,thicknessTargets,thicknessRows,coarse,feasibleRows,bestCoarse,aOpt,etaOpt,comp,cornerRows,cornerMin,validationSpecs,validationTables,validationActive,validationSummary,validationRows,activeValidation,jRob,pInfNom,pInfRob,cNom,cRob,convRows,cutoffRows,pressureRows,aVals,etaVals,p,kk,ne,decompRows,opacityRows,bVals,asymRows,offsets,aL,aR,pL,pR,spectrumRows,nominalKernels,mats,egs,i,eta,egBest,pBest,fig2Plot,fig3Plot,fig4Plot,fig5Plot,fig6Plot,fig7Plot,mapPressureRows,mapPressureDisplayRows,mapContrastRows,mapDensity,mapZeroContour,mapCavityContour,mapMarker,jMin,jMax,jDisplayMin},Print["--- Robust nonequilibrium Casimir optimization (Mathematica 9) ---"];
Print["Adachi KK completion: eps0Interband = ",N[eps0Interband,10],", deltaUV = ",N[deltaUV,10]];

(*Reduced band-gap optimization and material screen*)screenEtas={0.94,0.945,0.95};
screenRows=Table[{egBest,pBest}=bandGapScreen[eta];
{eta,egBest,1000. pBest},{eta,screenEtas}];
exportCSV["bandgap_optima.csv",{"eta","EgOpt_eV","PinfOpt_mPa"},screenRows];
mats={"InSb","InAs","InP","GaAs","ZnS"};
egs={0.170,0.354,1.344,1.424,3.680};
materialRows=Table[{mats[[i]],egs[[i]],1000. farFieldPressure[egs[[i]],0.94,300],1000. farFieldPressure[egs[[i]],0.945,300],1000. farFieldPressure[egs[[i]],0.95,300]},{i,1,Length[mats]}];
exportCSV["material_screen.csv",{"material","Eg_eV","Pinf_eta094_mPa","Pinf_eta0945_mPa","Pinf_eta095_mPa"},materialRows];

(*Thickness constraints*)
bEmission=findEmissionThickness[10^-3];
bOpt=Max[bEmission,9. aMax];
Print["Emission-only thickness at 10^-3 = ",10^6 bEmission," um"];
Print["Adopted thickness = ",10^6 bOpt," um"];
Print["Leakage at adopted b = ",opacityLeakage[bOpt,etaPhysicalMax,180,120]];
Print["Equilibrium coupling indicator = ",equilibriumCouplingIndicator[bOpt,aMax]];
thicknessTargets={10^-2,10^-3,10^-4};
thicknessRows=Table[bEmission=findEmissionThickness[target];{target,10^6 bEmission,10^6 aMax (target^(-1./3.)-1.),10^6 Max[bEmission,aMax (target^(-1./3.)-1.)]},{target,thicknessTargets}];
exportCSV["thickness_tolerance.csv",{"indicator_target","b_emission_um","b_equilibrium_um","b_required_um"},thicknessRows];

(*Coarse search is exported to document global lobe selection.*)coarse=coarseSearch[];
If[!numericTableQ[coarse,5],Print["FATAL: coarseSearch[] did not return a numeric N x 5 table."];Return[$Failed]];
exportCSV["fig4_coarse_robust_map.csv",{"a_m","eta","J_Pa","Ccav_Pa","cavityConstraintPass"},coarse];
feasibleRows=Select[coarse,#[[5]]==1.&];
If[Length[feasibleRows]==0,Print["FATAL: no feasible point was found in the coarse design map."];Return[$Failed]];
bestCoarse=First[Sort[feasibleRows,#1[[3]]>#2[[3]]&]];
Print["Best feasible coarse point: a = ",10^6 bestCoarse[[1]]," um, eta = ",bestCoarse[[2]],", J = ",1000. bestCoarse[[3]]," mPa"];

(*Refinement of the best lobe found in the complete declared design box.*)etaOpt=bestCoarse[[2]];
aOpt=refineBestGap[bestCoarse[[1]],etaOpt];
If[!NumberQ[aOpt],Print["FATAL: refineBestGap[] did not return a numeric gap."];Return[$Failed]];
comp=pressureComponents[aOpt,etaOpt,360,480,600,480];
If[!numericVectorQ[comp,4],Print["FATAL: pressureComponents[] did not return four numeric values."];Return[$Failed]];
cornerRows=robustCornerRows[aOpt,etaOpt,300,420,600,420];
cornerMin=Min[cornerRows[[All,3]]];
exportCSV["robust_corners.csv",{"delta_a_nm","delta_eta","P_mPa"},Map[{10^9 #[[1]],#[[2]],1000. #[[3]]}&,cornerRows]];

(*Resolve the continuous uncertainty rectangle numerically. The nested grids share nodes, so cached gap kernels are reused. The finest grid has 0.25 nm and 0.00025 resolution.*)
validationSpecs={{21,11},{41,21},{81,41}};
validationTables=Table[uncertaintyGridRows[aOpt,etaOpt,validationSpecs[[i,1]],validationSpecs[[i,2]],300,420,600,420],{i,1,Length[validationSpecs]}];
If[!And@@Map[numericTableQ[#,3]&,validationTables],Print["FATAL: continuous-uncertainty validation returned nonnumeric data."];Return[$Failed]];
validationActive=Map[Function[t,First[Sort[t,#1[[3]]<#2[[3]]&]]],validationTables];
validationSummary=Table[{validationSpecs[[i,1]],validationSpecs[[i,2]],10^9 validationActive[[i,1]],validationActive[[i,2]],1000. validationActive[[i,3]]},{i,1,Length[validationSpecs]}];
exportCSV["uncertainty_grid_convergence.csv",{"N_delta_a","N_delta_eta","active_delta_a_nm","active_delta_eta","Jgrid_mPa"},validationSummary];
validationRows=Last[validationTables];
activeValidation=Last[validationActive];
exportCSV["uncertainty_validation.csv",{"delta_a_nm","delta_eta","P_mPa"},Map[{10^9 #[[1]],#[[2]],1000. #[[3]]}&,validationRows]];
jRob=activeValidation[[3]];
If[jRob<cornerMin-10^-9,Print["WARNING: an interior uncertainty point is worse than the four corners; the dense-grid minimum has been used."]];
pInfNom=farFieldPressure[EgGaAs,etaOpt,420];
pInfRob=farFieldPressure[EgGaAs,etaOpt-deltaEta,420];
cNom=comp[[4]]-pInfNom;
cRob=jRob-pInfRob;
If[cRob<cavityContrastMin,Print["FATAL: the dense uncertainty minimum fails the robust cavity-contrast constraint."];Return[$Failed]];
If[Abs[validationActive[[-1,3]]-validationActive[[-2,3]]]>10^-8,Print["WARNING: uncertainty-grid refinement has not converged to 10^-8 Pa."]];
Print["OPTIMUM: aL=aR = ",10^6 aOpt," um; eta = ",etaOpt,"; b = ",10^6 bOpt," um"];
Print["Robust bilateral pressure = ",1000. jRob," mPa"];
Print["Active uncertainty point {delta a [nm], delta eta} = ",{10^9 activeValidation[[1]],activeValidation[[2]]}];
Print["Nominal {Peq, PPW, PEW, Ptotal} mPa = ",1000. comp];
Print["Robust cavity excess = ",1000. cRob," mPa"];
exportCSV["headline_results.csv",{"a_um","eta","b_um","Jrob_mPa","PinfRob_mPa","CcavRob_mPa","PinfNom_mPa","CcavNom_mPa","Peq_mPa","PPW_mPa","PEW_mPa","Ptotal_mPa"},{{10^6 aOpt,etaOpt,10^6 bOpt,1000. jRob,1000. pInfRob,1000. cRob,1000. pInfNom,1000. cNom,1000. comp[[1]],1000. comp[[2]],1000. comp[[3]],1000. comp[[4]]}}];


(****************************************)
(*Figure 2: Thickness Indicators*)
(****************************************)
(*Include 0.10 um explicitly because the manuscript contrasts it with the opaque-slab regime.*)
bVals=Join[{0.10*10^-6},Range[2.*10^-6,16.*10^-6,0.25*10^-6]];
opacityRows=Table[{10^6 b,opacityLeakage[b,etaPhysicalMax,100,64],equilibriumCouplingIndicator[b,aMax]},{b,bVals}];
exportCSV["fig2_opacity.csv",{"b_um","emissionLeakage","equilibriumIndicator"},opacityRows];
If[!numericTableQ[opacityRows,3],Print["FATAL: opacity data are not fully numeric."];Return[$Failed]];

(*Clean ticks for the ListLogPlot below*)
yMajorTicks = Table[{10^k, Superscript[10, k], {0.01, 0}}, {k, -4, 0}];
yIntermediateTicks = Flatten[Table[{n*10^k, "", {0.005, 0}}, {k, -4, -1}, {n, 2, 9}], 1];
allYTicks = Join[yMajorTicks, yIntermediateTicks];

(* Ticks for the right frame without labels *)
yMajorTicksRight = Table[{10^k, "", {0.01, 0}}, {k, -4, 0}];
yIntermediateTicksRight = Flatten[Table[{n*10^k, "", {0.005, 0}}, {k, -4, -1}, {n, 2, 9}], 1];
allYRightTicks = Join[yMajorTicksRight, yIntermediateTicksRight];

fig2Plot=ListLogPlot[{opacityRows[[All,{1,2}]],opacityRows[[All,{1,3}]]},Joined->True,Frame->True,Axes->False,
FrameLabel -> {
 Style[
  "Slab thickness \!\(\*StyleBox[\"b\",FontSlant->\"Italic\"]\) (\[Mu]m)",
  FontFamily -> "Times", FontSize -> 10
 ],
 Style[
  "Reduction-validity indicator",
  FontFamily -> "Times", FontSize -> 10
 ]
},
FrameTicksStyle -> Directive[FontFamily -> "Times", FontSize -> 10],
FrameTicks -> {{allYTicks, allYRightTicks}, {Automatic, Automatic}},
PlotLegends->Placed[LineLegend[{"Emission leakage","Equilibrium coupling indicator"}],{0.2,0.1}],PlotRange->All,ImageSize->500,
BaseStyle->{FontFamily->"Times",FontSize->10},PlotRangePadding->None];

(****************************************)
(*Figure 3: Pressure curves*)
(****************************************)
aVals=Range[0.60*10^-6,1.55*10^-6,0.02*10^-6];
etaVals={0.,0.90,0.925,0.94,0.945};
pressureRows=Flatten[Table[kk=pressureKernelCache[a,90,120,90];
p=If[eta==0.,peqCache[a,90],ne=noneqFromKernels[kk,eta];peqCache[a,90]+ne[[1]]+ne[[2]]];
{10^6 a,eta,1000. p},{a,aVals},{eta,etaVals}],1];
exportCSV["fig3_pressure_curves.csv",{"gap_um","eta","P_mPa"},pressureRows];
If[!numericTableQ[pressureRows,3],Print["FATAL: pressure curve data are not fully numeric."];Return[$Failed]];
fig3Plot=ListLinePlot[Table[Select[pressureRows,#[[2]]==eta&][[All,{1,3}]],{eta,etaVals}],
Frame->True,GridLines->{None,{0}},Axes->False,FrameLabel -> {
 Style[
  "Gap \!\(\*StyleBox[\"a\",FontSlant->\"Italic\"]\) (\[Mu]m)",
  FontFamily -> "Times", FontSize -> 10
 ],
 Style[
  "Outward pressure (mPa)",
  FontFamily -> "Times", FontSize -> 10
 ]
},
FrameTicksStyle -> Directive[FontFamily -> "Times", FontSize -> 10],
PlotLegends->Placed[Map[Function[x,"\[Eta] = "<>ToString[x]],etaVals],{0.1,0.8}],PlotRange->{-5,3.2},ImageSize->500,
BaseStyle->{FontFamily->"Times",FontSize->10},PlotRangePadding->None];

(****************************************)
(*Figure 4: robust design plane. The solid contour is J=0; the dashed
  contour is the imposed robust cavity-excess threshold Ccav=0.10 mPa.*)
(****************************************)
mapPressureRows=Map[{10^6 #[[1]],#[[2]],1000. #[[3]]}&,coarse];
mapContrastRows=Map[{10^6 #[[1]],#[[2]],1000. #[[4]]}&,coarse];
jMin=Min[mapPressureRows[[All,3]]];
jMax=Max[mapPressureRows[[All,3]]];
jDisplayMin=Max[jMin,-5.];
mapPressureDisplayRows=Map[{#[[1]],#[[2]],Min[jMax,Max[jDisplayMin,#[[3]]]]}&,mapPressureRows];
mapDensity=ListDensityPlot[mapPressureDisplayRows,InterpolationOrder->0,Frame->True,Axes->False,
FrameLabel -> {
 Style[
  "Nominal gap \!\(\*SubscriptBox[\(a\), \(0\)]\) (\[Mu]m)",
  FontFamily -> "Times", FontSize -> 10
 ],
 Style[
  "Nominal drive \!\(\*SubscriptBox[\(\[Eta]\), \(0\)]\)",
  FontFamily -> "Times", FontSize -> 10
 ]
},
FrameTicksStyle -> Directive[FontFamily -> "Times", FontSize -> 10],PlotRange->All,ColorFunction->"TemperatureMap",
ColorFunctionScaling->True,ImageSize->520];
mapZeroContour=ListContourPlot[mapPressureRows,Contours->{0.},ContourShading->False,
ContourStyle->{{Black,Thick}},InterpolationOrder->1];
mapCavityContour=ListContourPlot[mapContrastRows,Contours->{1000. cavityContrastMin},ContourShading->False,
ContourStyle->{{Black,Dotted,Thick}},InterpolationOrder->1];
mapMarker=Graphics[{Black,PointSize[0.018],Point[{10^6 aOpt,etaOpt}],White,PointSize[0.007],Point[{10^6 aOpt,etaOpt}]}];
fig4Plot=Style[Legended[Show[mapDensity,mapZeroContour,mapCavityContour,mapMarker,
PlotRange->All],Placed[BarLegend[{"TemperatureMap",{jDisplayMin,jMax}},
LegendLabel->Placed["Four-corner robust-pressure upper bound (mPa)",Right,Rotate[#,90 Degree]&],
LegendMarkerSize->{15,500}],Right]],
FontFamily->"Times",FontSize->10];

(*One-factor-at-a-time quadrature convergence. All unvaried orders are held at their final values.*)
convRows=Join[
Table[comp=pressureComponents[aOpt,etaOpt,n,480,600,480];Join[{"Nu",n},1000. comp],{n,{180,260,360}}],
Table[comp=pressureComponents[aOpt,etaOpt,360,n,600,480];Join[{"NE",n},1000. comp],{n,{180,300,480}}],
Table[comp=pressureComponents[aOpt,etaOpt,360,480,n,480];Join[{"NPW",n},1000. comp],{n,{200,400,600}}],
Table[comp=pressureComponents[aOpt,etaOpt,360,480,600,n];Join[{"NEW",n},1000. comp],{n,{180,300,480}}]
];
exportCSV["convergence.csv",{"varied_parameter","order","Peq_mPa","PPW_mPa","PEW_mPa","Ptotal_mPa"},convRows];

(*Energy-window and evanescent-cutoff convergence; no additional figure is needed.*)
cutoffRows=Join[
Table[comp=pressureComponentsCutoffs[aOpt,etaOpt,360,480,600,480,w,10.];Join[{"energy_window_eV",w},1000. comp],{w,{0.40,0.60,0.80}}],
Table[comp=pressureComponentsCutoffs[aOpt,etaOpt,360,480,600,480,0.60,al];Join[{"alpha_max",al},1000. comp],{al,{8.,10.,12.}}]
];
exportCSV["cutoff_convergence.csv",{"varied_parameter","value","Peq_mPa","PPW_mPa","PEW_mPa","Ptotal_mPa"},cutoffRows];

(****************************************)
(*Figure 5: Decomposition*)
(****************************************)
decompRows=Table[comp=pressureComponents[a,etaOpt,100,110,140,110];
{10^6 a,1000. comp[[1]],1000. comp[[2]],1000. comp[[3]],1000. comp[[4]]},{a,Range[0.70*10^-6,1.55*10^-6,0.02*10^-6]}];
exportCSV["fig5_decomposition.csv",{"gap_um","Peq_mPa","PPW_mPa","PEW_mPa","Ptotal_mPa"},decompRows];
If[!numericTableQ[decompRows,5],Print["FATAL: decomposition data are not fully numeric."];Return[$Failed]];
fig5Plot=ListLinePlot[Table[decompRows[[All,{1,j}]],{j,2,5}],Frame->True,Axes->False,
FrameLabel -> {
 Style[
  "Gap \!\(\*StyleBox[\"a\",FontSlant->\"Italic\"]\) (\[Mu]m)",
  FontFamily -> "Times", FontSize -> 10
 ],
 Style[
  "Pressure contribution (mPa)",
  FontFamily -> "Times", FontSize -> 10
 ]
},
FrameTicksStyle -> Directive[FontFamily -> "Times", FontSize -> 10],
PlotLegends->Placed[{"Equilibrium","Propagating excess","Evanescent excess","Total"},{0.85,0.2}],
PlotRange->{-3,4},ImageSize->500,BaseStyle->{FontFamily->"Times",FontSize->10},PlotRangePadding->None];

(****************************************)
(*Figure 6: Spectrum at optimum*)
(****************************************)
nominalKernels=noneqKernels[aOpt,300,400,300];
spectrumRows=Map[{#[[1]],#[[2]],1000. #[[3]] excessOccupation[#[[1]],etaOpt],1000. #[[4]] excessOccupation[#[[1]],etaOpt],1000. (#[[3]]+#[[4]]) excessOccupation[#[[1]],etaOpt]}&,nominalKernels];
exportCSV["fig6_spectrum_optimum.csv",{"E_eV","quadrature_weight_eV","dPPW_dE_mPaPerEV","dPEW_dE_mPaPerEV","dPNE_dE_mPaPerEV"},spectrumRows];
If[!numericTableQ[spectrumRows,5],Print["FATAL: optimum spectrum data are not fully numeric."];Return[$Failed]];
fig6Plot=ListLinePlot[{spectrumRows[[All,{1,3}]],spectrumRows[[All,{1,4}]],spectrumRows[[All,{1,5}]]},
Frame->True,Axes->False,FrameLabel -> {
 Style[
  "photon energy (eV)",
  FontFamily -> "Times", FontSize -> 10
 ],
 Style[
  "dP/dE (mPa/eV)",
  FontFamily -> "Times", FontSize -> 10
 ]
},
FrameTicksStyle -> Directive[FontFamily -> "Times", FontSize -> 10],
PlotLegends->Placed[{"Propagating","Evanescent","Sum"},{0.85,0.85}],PlotRange->{{1.42,2.02},{-5,122}},ImageSize->500,
BaseStyle->{FontFamily->"Times",FontSize->10},PlotRangePadding->None];

(****************************************)
(*Figure 7: Asymmetry diagnostic*)
(****************************************)
offsets=Range[-120.*10^-9,120.*10^-9,10.*10^-9];
asymRows=Table[aL=aOpt+off;aR=aOpt-off;
pL=pressureCached[aL,etaOpt,300,420,600,420];
pR=pressureCached[aR,etaOpt,300,420,600,420];
{10^9 off,1000. pL,1000. pR,1000. Min[pL,pR]},{off,offsets}];
exportCSV["fig7_asymmetry.csv",{"offset_nm","PL_mPa","PR_mPa","min_mPa"},asymRows];
If[!numericTableQ[asymRows,4],Print["FATAL: asymmetry data are not fully numeric."];Return[$Failed]];
fig7Plot=ListLinePlot[{asymRows[[All,{1,2}]],asymRows[[All,{1,3}]],asymRows[[All,{1,4}]]},
Frame->True,Axes->False,FrameLabel -> {
 Style[
  "Slab offset from symmetric design (nm)",
  FontFamily -> "Times", FontSize -> 10
 ],
 Style[
  "Outward pressure (mPa)",
  FontFamily -> "Times", FontSize -> 10
 ]
},
FrameTicksStyle -> Directive[FontFamily -> "Times", FontSize -> 10],
PlotLegends->Placed[{"\!\(\*SubscriptBox[\(P\), \(L\)]\)","\!\(\*SubscriptBox[\(P\), \(R\)]\)","min(\!\(\*SubscriptBox[\(P\), \(L\)]\),\!\(\*SubscriptBox[\(P\), \(R\)]\))"},{0.9,0.86}],
PlotRange->{2,3.1},ImageSize->500,BaseStyle->{FontFamily->"Times",FontSize->10},PlotRangePadding->None];

(*Display and export all six computational manuscript figures in order.*)
showAndExportPDF["fig2_opacity.pdf",fig2Plot];
showAndExportPDF["fig3_pressure_curves.pdf",fig3Plot];
showAndExportPDF["fig4_robust_design_map.pdf",fig4Plot];
showAndExportPDF["fig5_decomposition.pdf",fig5Plot];
showAndExportPDF["fig6_spectrum_optimum.pdf",fig6Plot];
showAndExportPDF["fig7_asymmetry.pdf",fig7Plot];
Print["Study complete. CSV data: ",dataDirectory];
Print["Figures: ",figureDirectory];
{aOpt,etaOpt,bOpt,jRob,comp}];

(*----------numerical self-tests for a Mathematica 9 run----------*)
SelfTest[]:=Module[{pinf,leak,eqi,pass},pinf=farFieldPressure[EgGaAs,0.95,360];
leak=opacityLeakage[13.5*10^-6,0.95,140,90];
eqi=equilibriumCouplingIndicator[13.5*10^-6,aMax];
Print["Self-test Pinf(GaAs, eta=.95) = ",1000. pinf," mPa (reference ~3.65028 mPa)"];
Print["Self-test leakage(b=13.5 um) = ",leak," (reference ~1.77e-4)"];
Print["Self-test equilibrium indicator = ",eqi," (reference 1e-3)"];
pass=Abs[1000. pinf-3.65028]<=0.02&&leak<=2.5*10^-4&&Abs[eqi-10^-3]<=10^-10;
If[pass,Print["SELF-TEST PASSED. It is safe to execute RunStudy[]."],Print["SELF-TEST FAILED. Do not execute RunStudy[] until the warnings are resolved."]];
pass];

Print["Definitions loaded. Execute SelfTest[] for a quick check and RunStudy[] for the full optimization."];
