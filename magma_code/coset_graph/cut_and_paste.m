intrinsic IsAdjacentSidepairing(side::SeqEnum) -> BoolElt
  {True if the sidepairing defined by side identifies sides that are adjacent.}
  gamma, initial, terminal := Explode(side);
  G := Parent(initial);
end intrinsic;
