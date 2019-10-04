SetDebugOnError(true);
load "config.m";
SetVerbose("Passport", true);

// user input
d := 47;
group := 5;
pass := 1;

// code
tdb := TransitiveGroups(d);
passports := PassportRepresentatives(tdb[group]);
passport := passports[pass];
sigma := passport[2][1];
name := Name(sigma);

G, sidepairing := CosetGraph(sigma);
str := CosetGraphToDotString(G);

WriteCosetGraphToFile(sigma);
