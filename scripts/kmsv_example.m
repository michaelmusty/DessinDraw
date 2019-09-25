SetDebugOnError(true);

// input triple
sigma0 := Sym(6)!(1,3,5,4)(2,6);
sigma1 := Sym(6)!(1,5,2,3,6,4);
sigmaoo := (sigma1*sigma0)^-1;
sigma := [sigma0, sigma1, sigmaoo];

// initialize graph
a := Order(sigma0);
b := Order(sigma1);
c := Order(sigmaoo);
d := Degree(Parent(sigma0));
F3<da,db,dc> := FreeGroup(3);
Delta<da,db,dc> := quo< F3 | [da^a, db^b, dc^c, dc*db*da] >;
G := MultiDigraph<d | >;
AssignLabel(~G, Vertices(G)[1], Delta!1);
sidepairing := [];
mon := sub< Sym(d) | sigma>;
pi := hom< Delta -> mon | sigma>;

// first iteration
j := 1;
vj := Vertices(G)[j];
alphaj := Label(vj);
epses := [da, da^-1, db, db^-1];
for eps in epses do
  printf "\n";
  printf "eps = %o\n", eps;
  printf "alphaj = %o\n", alphaj;
  printf "alphajeps = %o\n", alphaj*eps;
  i := 1^pi(alphaj*eps);
  printf "i=1^pi(alphaj*eps)=%o\n", i;
  vi := Vertices(G)[i];
  AddEdge(~G, vj, vi, eps);
  printf "edge vertex %o -> vertex %o labelled %o added\n", j, i, eps;
  if IsLabelled(vi) then
    alphai := Label(vi);
    printf "vertex %o is labelled %o\n", i, alphai;
    gamma := alphaj*eps*alphai^-1;
    printf "alphaj=%o\n", alphaj;
    printf "eps=%o\n", eps;
    printf "alphai=%o\n", alphai;
    printf "alphai^-1=%o\n", alphai^-1;
    printf "gamma=alphaj*eps*alphai^-1=%o\n", gamma;
    if not (gamma eq Identity(Parent(gamma))) then
      printf "IsIdentity(%o) = false\n", gamma;
      printf "new side pairing element:\n";
      printf "gamma=%o\n", gamma;
      printf "(j,eps)=(%o,%o)\n", j, eps;
      printf "(i,eps^-1)=(%o,%o)\n", i, eps^-1;
      Append(~sidepairing, [* gamma, <j, eps>, <i, eps^-1> *]);
    else
      printf "IsIdentity(%o) = true\n", gamma;
      printf "no side pairing added\n";
    end if;
  else
    alphai := alphaj*eps;
    printf "vertex %o is not labelled\n", i;
    printf "new label for vi = alphai = %o\n", alphai;
    AssignLabel(~G, vi, alphai);
  end if;
end for;
