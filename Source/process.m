// Make LMFDB Artin database, TD Apr 2013
// Modified JJ Nov 2015


load "nfdb.m";

Process(x^4-2);
Process(x^8+2);
Process(x^8-2);
Process(x^8-3);
Process(x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 10*x - 9);
Process(x^8 - x^7 + 7*x^6 + 7*x^5 + 7*x^4 + 7*x^3 + 7*x^2 + 5*x + 2);

ProcessNFDBdatFile(LMFDBdir*"nfdb_bordeaux.dat"); // Bordeaux database
//ProcessDirectory(OldArtdir);

//IndexDatabase();
//ExportToJSON(: test:=false);

/*
Statistics();
*/

