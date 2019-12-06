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



set
         
         scen /s1*s250/
         w    /w1*w4/
;

Alias(Time,t);
// outcomment next line, if you want to reseed the random number generator at each new run of the file
*execseed=gmillisec(jnow)

SCALAR
        randnum
        ;

Parameter RetScenWeeks(i,scen,w);
Parameter RetScen(i, scen);


loop(scen,
    loop(w, 
        randnum = uniformint(5, 340);
        RetScenWeeks(i, scen,  w) = SUM(t$(ord(t)=randnum), IndexReturns(i, t)) ;
        );
        RetScen(i,scen) = prod(w, (1 + RetScenWeeks(i, scen,w ))) - 1;
    );


display RetScen;
EXECUTE_UNLOAD 'Scenarios02.gdx', scen, Asset, RetScen;
