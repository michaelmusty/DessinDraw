SetDebugOnError(true);
load "config_without_triangles.m";

db := TransitiveGroups(5);
passports := PassportRepresentatives(db[4]);
passport := passports[3];
sigma := passport[2][1];

sigma0, sigma1, sigmaoo := Explode(sigma);
a := Order(sigma0);
b := Order(sigma1);
c := Order(sigmaoo);
d := Degree(Parent(sigma0));
F3<da,db,dc> := FreeGroup(3);
Delta<da,db,dc> := quo< F3 | [da^a, db^b, dc^c, dc*db*da] >;
mon := sub< Sym(d) | sigma>;
pi := hom< Delta -> mon | sigma>;

G := MultiDigraph<d | >;
AssignLabel(~G, Vertices(G)[1], Delta!1);
sidepairing := [];

/* CosetGraph(sigma); */
