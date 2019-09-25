SetDebugOnError(true);
load "config.m";

db := TransitiveGroups(5);
passports := PassportRepresentatives(db[4]);
passport := passports[3];
sigma := passport[2][1];

G, sidepairing := CosetGraph(sigma);
