// Make LMFDB Artin database, TD Apr 2013
// Modified JJ Nov 2015


load "nfdb.m";

IndexDatabase();


//  Export the database to json, write statistics to stats.inc

ExportToJSON(: test:=false);

/*
Statistics();
*/

