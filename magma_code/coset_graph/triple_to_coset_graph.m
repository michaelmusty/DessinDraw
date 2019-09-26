intrinsic CosetGraph(sigma::SeqEnum[GrpPermElt]: Al := "Full") -> Any
  {Returns a coset graph for the triangle subgroup Gamma corresponding to the permutation triple sigma. Gamma is an index d subgroup of the triangle group Delta(a,b,c) where a,b,c are the orders of sigma0,sigma1,sigmaoo respectively and d is the degree of the ambient symmetric group. If Al eq "Petal", then prefer 'a' moves; otherwise, "Full" gives 'a' and 'b' moves.}
  sigma0, sigma1, sigmaoo := Explode(sigma);
  a := Order(sigma0);
  b := Order(sigma1);
  c := Order(sigmaoo);
  d := Degree(Parent(sigma0));
  F3<da,db,dc> := FreeGroup(3);
  Delta<da,db,dc> := quo< F3 | [da^a, db^b, dc^c, dc*db*da] >;
  G := MultiDigraph<d | >;
  AssignLabel(~G, Vertices(G)[1], Delta!1);
  mon := sub< Sym(d) | sigma>;
  pi := hom< Delta -> mon | sigma>;
  sidepairing := [];
  if Al eq "Full" then
    frontier := [1];
    epses := [da, da^-1, db, db^-1];
    while not IsEmpty(frontier) do
      j := frontier[1];
      Remove(~frontier, 1);
      vj := Vertices(G)[j];
      alphaj := Label(vj);
      for eps in epses do
        i := 1^(pi(alphaj*eps));
        vi := Vertices(G)[i];
        /* AddEdge(~G, vj, vi, eps); */
        if IsLabelled(vi) then
          alphai := Label(vi);
          gamma := alphaj*eps*alphai^-1;
          if not (gamma eq Identity(Parent(gamma))) then
            AddEdge(~G, vj, vi, [* eps, "side" *]);
            /* Append(~sidepairing, [* gamma, <j, eps>, <i, eps^-1> *]); */
            Append(~sidepairing, [* gamma, <vj, eps>, <vi, eps^-1> *]);
          else
            AddEdge(~G, vj, vi, [* eps, "interior" *]);
          end if;
        else
          AddEdge(~G, vj, vi, [* eps, "interior" *]);
          AssignLabel(~G, vi, alphaj*eps);
          Append(~frontier, i);
        end if;
      end for;
    end while;
  end if;
  /* else // Al eq "Petal" */
    /* frontierA := [1]; */
    /* frontierB := [1]; */

    /* finishAllAs := function(G, frontierA, frontierB, sidepairing); */
    /*   while not IsEmpty(frontierA) do */
    /*     dists := [Distance(D0, Label(Vertices(G)[i])*D0) : i in frontierA]; */
    /*     _, jind := Min(dists); */
    /*     j := frontierA[jind]; */
    /*     vprintf Shimura : "A: Taking j = %o, with distance %o and label %o\n", */
    /*         j, RealField(6)!dists[jind], Label(Vertices(G)[j]); */
    /*     Remove(~frontierA, jind); */

    /*     alphaj := Label(Vertices(G)[j]); */
    /*     donepos := false; */
    /*     doneneg := false; */
    /*     ind := 0; */
    /*     // rotate around the chosen vertex */
    /*     while not (donepos and doneneg) do */
    /*       ind +:= 1; */
    /*       epspows := []; */
    /*       if not donepos then */
    /*         Append(~epspows, ind); */
    /*       end if; */
    /*       if not doneneg then */
    /*         Append(~epspows, -ind); */
    /*       end if; */
    /*       for epspow in epspows do */
    /*         eps := Delta.1^Sign(epspow); */
    /*         jp := j^pi(Delta.1^(epspow-Sign(epspow))); */
    /*         alphajp := Label(Vertices(G)[jp]); */
    /*         i := 1^(pi(alphajp*eps)); */
    /*         AddEdge(~G, Vertices(G)[jp], Vertices(G)[i], eps); */

    /*         vprintf Shimura : "A: Rotating by delta_a^%o, ", epspow; */

    /*         if IsLabelled(Vertices(G)[i]) then */
    /*           alphai := Label(Vertices(G)[i]); */

    /*           gamma := alphajp*eps*alphai^-1; */
    /*           if not IsScalar(Quaternion(gamma)) then */
    /*             Append(~sidepairing, [* gamma, <jp, eps>, <i, eps^-1> *]); */
    /*           end if; */
    /*           if epspow gt 0 then */
    /*             donepos := true; */
    /*           end if; */
    /*           if epspow lt 0 then */
    /*             doneneg := true; */
    /*           end if; */

    /*           vprintf Shimura : "already identified %o (label %o)", i, Label(Vertices(G)[i]); */

    /*           if not IsScalar(Quaternion(gamma)) then */
    /*             vprintf Shimura : ", sidepairing %o\n", gamma; */
    /*           else */
    /*             vprintf Shimura : "\n"; */
    /*           end if; */
    /*         else */
    /*           AssignLabel(~G, Vertices(G)[i], alphajp*eps); */
    /*           Append(~frontierB, i); */

    /*           vprintf Shimura : "new coset %o\n", i; */
    /*         end if; */
    /*       end for; */
    /*     end while; */
    /*   end while; */
    /*   return G, frontierA, frontierB, sidepairing; */
    /* end function; */

    /* // Build basic graph */
    /* while not IsEmpty(frontierA cat frontierB) do */
    /*   G, frontierA, frontierB, sidepairing := finishAllAs(G, frontierA, frontierB, sidepairing); */

    /*   // Now try a "B" */
    /*   dists := [Distance(D0, Label(Vertices(G)[i])*D0) : i in frontierB]; */
    /*   _, jind := Min(dists); */
    /*   j := frontierB[jind]; */
    /*   vprintf Shimura : "B: Taking j = %o, with distance %o and label %o\n", */
    /*         j, RealField(6)!dists[jind], Label(Vertices(G)[j]); */
    /*   Remove(~frontierB, jind); */

    /*   alphaj := Label(Vertices(G)[j]); */
    /*   for epspow in [1,-1] do */
    /*     eps := Delta.2^epspow; */
    /*     vprintf Shimura : "B: Rotating by delta_b^%o, ", epspow; */

    /*     i := 1^(pi(alphaj*eps)); */
    /*     AddEdge(~G, Vertices(G)[j], Vertices(G)[i], eps); */
    /*     if IsLabelled(Vertices(G)[i]) then */
    /*       alphai := Label(Vertices(G)[i]); */
    /*       gamma := alphaj*eps*alphai^-1; */
    /*       if not IsScalar(Quaternion(gamma)) then */
    /*         Append(~sidepairing, [* gamma, <j, eps>, <i, eps^-1> *]); */
    /*       end if; */

    /*       vprintf Shimura : "already identified %o (label %o)", i, Label(Vertices(G)[i]); */
    /*       if not IsScalar(Quaternion(gamma)) then */
    /*         vprintf Shimura : ", sidepairing %o\n", gamma; */
    /*       else */
    /*         vprintf Shimura : "\n"; */
    /*       end if; */
    /*     else */
    /*       AssignLabel(~G, Vertices(G)[i], alphaj*eps); */
    /*       Append(~frontierA, i); */
    /*       Append(~frontierB, i); */

    /*       vprintf Shimura : "new coset %o\n", i; */
    /*     end if; */

    /*     G, frontierA, frontierB, sidepairing := finishAllAs(G, frontierA, frontierB, sidepairing); */
    /*   end for; */
    /* end while; */
  /* end if; */
  return G, sidepairing;
end intrinsic;
