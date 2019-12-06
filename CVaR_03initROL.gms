$TITLE Value at Risk and Conditional Value at Risk models
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
option optcr=0, reslim=120;

SET
         scen
         Asset
         p
;
ALIAS(Asset,i);
ALIAS(scen, s);


PARAMETER
         RetScen(i, s)
         AllRetScen(p,i,s)
         
;


$GDXIN RScen
$LOAD p, scen, Asset, AllRetScen
$GDXIN


Parameter
        LastHold(i) 'last hold of assets'
        InitHold(i)  'be';
        
$GDXIN Inithold
$LOAD LastHold=x.l
$GDXIN

SCALARS

        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit'
;


alpha  = 0.95;


PARAMETERS
        pr(s)       'Scenario probability'
        FP(i,s)      'Final values'
        EP(i)       'Expected final values'
       
        
;






POSITIVE VARIABLES
         x(i)            'Holdings of assets in monetary units (not proportions)'
         VaRDev(s)       'Measures of the deviation from the VaR'
         buy(i)          'Amount to buy'
         sell(i)         'Amount to be sold'
         x_0(i)          'Holdings for the flat cost regime'
         x_1(i)          'Holdings for the linear cost regime'
         
BINARY VARIABLE
    Y(i)          'Indicator variable for assets included in the portfolio';

PARAMETER
    xlow(i)    'lower bound for active variables' ;

// In case short sales are allowed these bounds must be set properly.
xlow(i) = 0.0;
x.UP(i) = 100000.0;
;

SCALARS
  FlatCost 'Normalized fixed transaction cost' / 20 /
  PropCost 'Normalized proportional transaction cost' / 0.001 /;
  

VARIABLES
         losses(s)       'The scenario loss function'
         VaR             'The alpha Value-at-Risk'
         CVaR            'The alpha Conditional Value-at-Risk'
         ExpectedReturn  'Expected return of the portfolio'
         obj             'objective function value'
         PortReturn
;

x_0.UP(i) = 20000;

Parameters
    MuTarget,
    CVaRTarget,
    BestCase,
    WorstCase
    
;



EQUATIONS
         BalanceCon(i)    'Same sell and buy amount'
       
         ReturnCon       'Equation defining the portfolio expected return'
         LossDefCon(s)   'Equation defining the losses'
         VaRDevCon(s)    'Equation defining the VaRDev variable'
         CVaRDefCon      'Equation defining the CVaR'
         ObjectivFunc    'lambda formulation of the MeanCVaR model'
         CVaRLimCon      'Constraint limiting the CVaR'
         ReturnLimCon    'Constraint on a minimum expected return'
         
         
         HoldingCon(i)        'Constraint defining the holdings'
         ReturnDefWithCost    'Equation defining the portfolio return with cost'
         FlatCostBounds(i)    'Upper bounds for flat transaction fee'
         LinCostBounds(i)     'Upper bonds for linear transaction fee'
;



*--Objective------

*--s.t.-----------
BalanceCon(i)..          x(i) - buy(i) + sell(i) =e=  LastHold(i);


ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, FP(i, s)*x(i) );

VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

ReturnDefWithCost..      PortReturn =e= SUM(i, ( EP(i)*x_0(i) - FlatCost*Y(i) ) ) +
                          SUM(i, (EP(i) - PropCost)*x_1(i));

ObjectivFunc ..          Obj =E= (1-lambda)*PortReturn - lambda*CVaR;

CVaRLimCon ..            CVaR =L= CVaRLim;

ReturnLimCon ..          ExpectedReturn =G= ExpRetLim;


HoldingCon(i)..           x(i) =e= x_0(i) + x_1(i); 

FlatCostBounds(i)..       x_0(i) =l= x_0.UP(i) * Y(i);

LinCostBounds(i)..        x_1(i) =l= Y(i);



 

*--Models-----------
MODEL RevisionCVaRModel 'The Conditional Value-at-Risk Model' /BalanceCon,ReturnCon,LossDefCon,VaRDevCon,CVaRDefCon,
ObjectivFunc,CVaRLimCon,ReturnLimCon,HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds /;


loop(p$(ord(p)>1 and ord(p)<4),
    pr(s) = 1.0 / CARD(s);

    Retscen(i,s) = AllRetScen(p,i,s);
    InitHold(i) = x.l*SUM(s, pr(s) * (AllRetScen(p-1,i,s)+1));
    ExpRetLim = -SUM(i,InitHold(i));
    CVaRLim = SUM(i,InitHold(i))*10;
    
   
    FP(i,s) = 1 + RetScen( i, s );
    EP(i) = SUM(s, pr(s) * FP(i,s));

     DISPLAY  EP;
    
       );
$exit 
    X.fx(i) = CVaRLim/180;
   

    lambda = 0.999;
    SOLVE RevisionCVaRModel Maximizing obj Using MIP;
    MuTarget = ExpectedReturn.l;
    CVaRTarget = CVaR.L;

*-------------------------------------
    X.lo(i) = 0;
    X.up(i) = CVaRLim;

//Let's minimize CVaR with the target return from the equal weight portfolio


Lambda = 1;
CVaRLim = CVaRLim+100000;
SOLVE CVaRModel Maximizing OBJ Using MIP;
display X.l, ExpectedReturn.l,sell.l, buy.l,CVaR.l;


//Let's maximize expected return with the target CVaR from the equal weight portfolio
CVaRLim = CVaRtarget;


Lambda = 0;
ExpRetLim = -100000;
SOLVE CVaRModel Maximizing OBJ Using MIP;
display X.l, ExpectedReturn.l, CVaR.l;



















