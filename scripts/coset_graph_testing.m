SetDebugOnError(true);
load "config.m";

// user input
d := 4;
group := 5;
pass := 2;

// code
tdb := TransitiveGroups(d);
passports := PassportRepresentatives(tdb[group]);
passport := passports[pass];
sigma := passport[2][1];
name := Name(sigma);

G, sidepairing := CosetGraph(sigma);
str := CosetGraphToDotString(G);

WriteCosetGraphToFile(sigma);
