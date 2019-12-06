$TITLE sccenario generation
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
option optcr=0, reslim=120;


SET Asset Underlying asset;

SET Date Time periods;

ALIAS(Asset,i);
ALIAS(Date,t);
Parameter
    IndexReturns(i,t) 'Index return';

$GDXIN Data
$LOAD Asset, Date, IndexReturns
$GDXIN

display i,t, IndexReturns;


set
       
*       t /3-8-1998, 31-8-1998/   t /3-8-1998, 31-8-1998, 28-9-1998/
         scen /s1*s250/
         w    /w1*w4/
         period(t) /3-1-2005/
;


// outcomment next line, if you want to reseed the random number generator at each new run of the file
*execseed=gmillisec(jnow)

SCALAR
        randnum
        start
        end;

Parameter RetScenWeeks(i,scen,w);
Parameter RetScen(t, i, scen);
Parameter RollingScen(t);

loop(period,
    start = 1+4*(ord(period)-1);
    end = 336+4*(ord(period)-1);
    loop(scen,
             loop(w, 
                     randnum = uniformint(start, end);
                     RetScenWeeks(i, scen,  w) = SUM(t$(ord(t)=randnum), IndexReturns(i, t)) ;
             );
    );
    RetScen(period, i,scen) = prod(w, (1 + RetScenWeeks(i, scen,w ))) - 1;
 );

display RetScen;
EXECUTE_UNLOAD 'Scenarios01.gdx', period, scen, Asset, RetScen;
