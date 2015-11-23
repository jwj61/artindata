// Make LMFDB Artin database, TD Apr 2013
// Modified JJ Nov 2015


load "nfdb.m";

//Process(x);
//Process(x^4-2);
//Process(x^4+2);
//Process(x^8+2);
//Process(x^8-2);
//Process(x^8-3);
//Process(x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 10*x - 9);
//Process(x^8 - x^7 + 7*x^6 + 7*x^5 + 7*x^4 + 7*x^3 + 7*x^2 + 5*x + 2);
//Process(x^6 - 3*x^5 + 5*x^4 - 5*x^3 + 5*x^2 - 3*x + 1);
//Process(x^17 - 2*x^16 + 8*x^13 + 16*x^12 - 16*x^11 + 64*x^9 - 32*x^8 - 80*x^7 + 32*x^6 + 40*x^5 + 80*x^4 + 16*x^3 - 128*x^2 - 2*x + 68);

//ProcessNFDBdatFile(LMFDBdir*"nfdb_bordeaux.dat"); // Bordeaux database
ProcessDirectory(OldArtdir);

//IndexDatabase();
//ExportToJSON(: test:=false);

/*
Statistics();
*/

