intrinsic IsEdgeSidepairing(e::GrphEdge) -> BoolElt
  {True if e is a side and false if e is an interior edge.}
  return Label(e)[2] eq "side";
end intrinsic;

intrinsic CosetGraphToDotString(G::GrphMultDir) -> MonStgElt
  {Given a coset graph G, output a string to produce a graphviz file in the dot language.}
  str := "digraph G {\n";
  for i := 1 to #Vertices(G) do
    v := Vertices(G)[i];
    if i eq 1 then
      str *:= Sprintf("v_%o [ label = \"id\" ];\n", i);
    else
      str *:= Sprintf("v_%o [ label = \"%o\" ];\n", i, Label(v));
    end if;
  end for;
  for e in Edges(G) do
    if IsEdgeSidepairing(e) then
      str *:= Sprintf("v_%o -> v_%o [ color=lightgray, label = \"%o\" ];\n", InitialVertex(e), TerminalVertex(e), Label(e)[1]);
    else
      str *:= Sprintf("v_%o -> v_%o [ label = \"%o\" ];\n", InitialVertex(e), TerminalVertex(e), Label(e)[1]);
    end if;
  end for;
  str *:= "}";
  return str;
end intrinsic;

intrinsic WriteCosetGraphToFile(sigma::SeqEnum[GrpPermElt]) -> Any
  {}
  G, S := CosetGraph(sigma);
  str := CosetGraphToDotString(G);
  dir := GetCurrentDirectory();
  name := Name(sigma);
  filename := dir cat "/graphviz/" cat name cat ".gv";
  Write(filename, str : Overwrite := true);
  returnText := Sprintf("%o.gv written to file\n%o", name, filename);
  return returnText;
end intrinsic;
