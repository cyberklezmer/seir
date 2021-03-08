Seir Filter

preliminary and very incomplete documentation. 

1. Preparing input dataset - seirwork/src.ods
(under construction)

Uses files from https://onemocneni-aktualne.mzcr.cz/api/v2/covid-19 and file 
modely_05_datumy.csv from restricted dataset (apply here https://docs.google.com/forms/d/e/1FAIpQLSfncCCbPngFtfdHV8XMBcU5lGiPMMRGj-6BYVLw2Nj6PvXOFA/viewform)

- decide on estimation and evaluation horizon and put the dates into sheet $dates
- run uzis2uzis(todays date) from src/new.cpp, ane copy result to $uzis (adjust path accordingly in src/rates.cpp)
- run mzcr2mzcr(todays date), ane copy result to $mzcr (osoby.csv and umrti.csv needed, adjust path accordingly)
- supply hospital detection probabilities (from PES dataset) to $pdet
- copy testy-pcr-antigenni.csv to $testy.csv (be sure to preserve positions of dates)
- copy content nakazeni-vyleceni-umrti-testy.csv to $nvut.csv (be sure to preserve dates positions, i.e. copy only the rows starting from 24.2.2020)
- copy hospitalizace.csv to $hospitalizace.csv (be sure to keep dates' position)
- copy paq_data.csv to $paq_data.csv (chagnes bi-weekliy, on demand at bisop)
- possibly adjust column Scenario R0 in $work_paq_week by values of R0 (if not filled, computed from data
- if, errors are shown $epi.csv then decrease the estimation horizon accordingly 
- save sheet epi.csv as your input file (as csv) to folder input.

2. Running the model
(under construction)

Run either run4 (cohort model) or run1 (single model). What the procedures do is determined by the #prg argument. It is possible to 
- determine the initial value of final estimation
- estimate using wls (with initial params from prg)
- reestimate using wls with regularization for getting unit standard residuals
- output the result (using init parameters if not estiamted). Output is done to files X_output1.csv (one day predictions) and X_output7.csv (7 day predistions). 
Remark: in the initial parameters list of prg, the logical variables determine whether the parameter is estimated (false) or is fixed (true)

3. Displaying the results
(under construction)

use seirwork/result1.ods (single model) or  seirwork/result4.ods (cohort model). COpy the X_output_7.csv to $output7 and see the graphs 





