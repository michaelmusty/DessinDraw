SetDebugOnError(true);
load "config.m";
SetVerbose("Passport", true);

// user input
sigma := [Sym(3)|(1,2),(2,3),(1,3,2)];
name := Name(sigma);

G, sidepairing := CosetGraph(sigma);
str := CosetGraphToDotString(G);

WriteCosetGraphToFile(sigma);
