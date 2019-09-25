load "config.m";
SetDebugOnError(true);

// input triple
sigma0 := Sym(6)!(1,3,5,4)(2,6);
sigma1 := Sym(6)!(1,5,2,3,6,4);
sigmaoo := (sigma1*sigma0)^-1;
sigma := [sigma0, sigma1, sigmaoo];

G, sidepairing := CosetGraph(sigma);
