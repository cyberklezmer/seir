Replication of computations from "SEIF Filter - Stochasti Epidemic Model"

To get materials for graphs:
1. Build project described by crc/cmakelist.cpp
2. Run it
3. See results in filed output*.csv (in the build directory)

To run or estimated with different data
1. Change seirwork/src.ods accordingly (files mzcr.csv and uzis.csv are created by procedures mzcr2mzcr and uzis2uzis, the latter requiring no-public data, in src/rates.cpp, the rest of inputs is publicly available)
2. Save the sheet "epi.csv" into input folder
3. Change g4 (or perhaps g) accordingly in new.cpp
4. Run procedure nw (new.cpp)
 
