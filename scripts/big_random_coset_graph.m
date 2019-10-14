SetDebugOnError(true);
load "config.m";
SetVerbose("Passport", true);

// user input
d := 100;

// code
sigma0 := Random(Sym(d));
sigma1 := Random(Sym(d));
sigmaoo := (sigma1*sigma0)^-1;
assert sigmaoo*sigma1*sigma0 eq Id(Sym(d));
sigma := [sigma0, sigma1, sigmaoo];
time assert IsTransitive(sub<Sym(d)|sigma>);
/* time passports := PassportRepresentatives(sigma); */
/* passport := passports[1]; */
/* name := Name(sigma); */

/* time G, sidepairing := CosetGraph(sigma); */
/* str := CosetGraphToDotString(G); */

/* WriteCosetGraphToFile(sigma); */
