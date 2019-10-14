load "config.m";
Attach("belyi_code/nonhyperbolic.m");

d := 12;
abc := [2,3,6];
/* abc := [2,4,4]; */
time passports := PassportRepresentatives(d : abc := abc);

triples := [];
for passport in passports do
  l := passport[2];
  for i := 1 to #l do
    triples cat:= l[i][2];
  end for;
end for;

good := [];
for sigma in triples do
  nonuniform := [];
  for perm in sigma do
    if #CycleStructure(perm) gt 1 then
      Append(~nonuniform, perm);
    end if;
  end for;
  if #nonuniform eq 1 then
    Append(~good, sigma);
  end if;
end for;

#good;

curves := [* *];
maps := [* *];
for i := 1 to #triples do
  X,phi := EuclideanBelyiMap(triples[i]);
  Append(~curves, X);
  Append(~maps, phi);
end for;
