SetDebugOnError(true);
load "config.m";

// user input
d := 2;
group := 1;
pass := 1;

// code
db := TransitiveGroups(d);
passports := PassportRepresentatives(db[group]);
passport := passports[pass];
sigma := passport[2][1];
name := Name(sigma);

G, sidepairing := CosetGraph(sigma);
str := CosetGraphToDotString(G);

WriteCosetGraphToFile(sigma);
