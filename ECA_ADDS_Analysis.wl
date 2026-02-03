(* ::Package:: *)

(* ECA Block-Sequential Analysis - Supplementary Code *)
(* D. Gimigliano, P.P. Balbi - Universidade Presbiteriana Mackenzie *)
(* Based on CAMaT (Cellular Automata Mathematica Toolbox) *)

(* ============================================================ *)
(* PARAMETERS *)
(* ============================================================ *)

latticeSize = 5;
radius = 1;

(* ============================================================ *)
(* RULE REPRESENTATION *)
(* ============================================================ *)

RuleSpaceSize[k_Integer:2,r_:1] := k^k^\[LeftFloor]2r+1\[RightFloor];

AllNeighbourhoods[k_Integer:2,r_:1] := Tuples[Range[k-1,0,-1],\[LeftFloor]2r+1\[RightFloor]];

RuleTableFromKAry[kAryRuleTable_,k_Integer:2,r_:1] :=
  MapThread[List[#1,#2]&, {Tuples[Range[k-1,0,-1],\[LeftFloor]2r+1\[RightFloor]], kAryRuleTable}];

RuleTable[rnum_Integer,k_Integer:2,r_:1] :=
  RuleTableFromKAry[PadLeft[IntegerDigits[rnum,k],k^(\[LeftCeiling]2r\[RightCeiling]+1)],k,r];
RuleTable[{rnum_Integer,k_Integer:2,r_:1}] := RuleTable[rnum,k,r];

KAryFromRuleTable[ruleTable_] := #[[2]]& /@ ruleTable;

RnumFromRuleTable[ruletable_,k_Integer:2] := FromDigits[KAryFromRuleTable[ruletable],k];

(* ============================================================ *)
(* DYNAMICAL EQUIVALENCE *)
(* ============================================================ *)

AllStatePermutations[k_Integer:2] :=
  Module[{permuts=Permutations[Range[0,k-1]]},
    MapThread[Thread[#1->#2]&, {Table[First@permuts,{Length[permuts]-1}], Rest@permuts}]];

ConjugateTransform[statetransition_,k_Integer:2,r_:1] :=
  (statetransition/.#)& /@ AllStatePermutations[k];

ReflectedTransform[statetransition_] :=
  {Reverse[statetransition[[1]]], statetransition[[2]]};

ConjugateReflectedTransform[statetransition_,k_Integer:2,r_:1] :=
  ConjugateTransform[ReflectedTransform[statetransition],k,r];

ConjugateRules[rnum_Integer,k_Integer:2,r_:1] :=
  {#,k,r}& /@ Union[RnumFromRuleTable[SortBy[#,k-#[[1]]&],k]& /@
    Transpose[ConjugateTransform[#,k,r]& /@ RuleTable[rnum,k,r]]];
ConjugateRules[{rnum_Integer,k_Integer:2,r_:1}] := ConjugateRules[rnum,k,r];

ReflectedRule[rnum_Integer,k_Integer:2,r_:1] :=
  If[IntegerQ[r]===True,
    {RnumFromRuleTable[SortBy[ReflectedTransform /@ RuleTable[rnum,k,r], k-#[[1]]&], k], k, r},
    {}];
ReflectedRule[{rnum_Integer,k_Integer:2,r_:1}] := ReflectedRule[rnum,k,r];

ConjugateReflectedRules[rnum_Integer,k_Integer:2,r_:1] :=
  If[IntegerQ[r]===True,
    Union[{RnumFromRuleTable[SortBy[#,k-#[[1]]&],k],k,r}& /@
      Transpose[ConjugateReflectedTransform[#,k,r]& /@ RuleTable[rnum,k,r]]],
    {}];
ConjugateReflectedRules[{rnum_Integer,k_Integer:2,r_:1}] := ConjugateReflectedRules[rnum,k,r];

DynamicalEquivalenceClass[rnum_Integer,k_Integer:2,r_:1] :=
  DeleteCases[Union[Join[{{rnum,k,r}}, ConjugateRules[rnum,k,r],
    {ReflectedRule[rnum,k,r]}, ConjugateReflectedRules[rnum,k,r]]], {}];
DynamicalEquivalenceClass[{rnum_Integer,k_Integer:2,r_:1}] := DynamicalEquivalenceClass[rnum,k,r];

(* ============================================================ *)
(* INTERACTION GRAPH *)
(* ============================================================ *)

CAInteractionGraph[r_:1,configLen_Integer] :=
  Graph[Flatten[Table[(First[x]\[DirectedEdge]#)& /@ Drop[x,1],
    {x, RotateLeft[#,\[LeftCeiling]r\[RightCeiling]]& /@
      Partition[Range[1,configLen], IntegerPart[2*r+1], 1, 1]}]]];

(* ============================================================ *)
(* UPDATE SCHEDULE EQUIVALENCE CLASSES *)
(* ============================================================ *)

Subdigraph[intG_Graph,C_List,D_List] :=
  If[C=={}||D=={},
    If[C=={}&&D=={}, Graph[{}],
      With[{elist=Sort[Intersection[EdgeList[intG],
        Flatten[Table[DirectedEdge[x,y],{x,Union[C,D]},{y,Union[C,D]}],1]]]},
        Graph[elist, EdgeWeight->Table["+",{Length[elist]}],
          EdgeLabels->"+", EdgeStyle->Red, VertexLabels->"Name"]]],
    With[{m=Table[Flatten[{x,If[MemberQ[C,First[x]]&&MemberQ[D,Last[x]],
        {"-",Darker[Green,0.8]},{"+",Red}]}],
        {x,Sort[Intersection[EdgeList[intG],
          Flatten[Table[DirectedEdge[x,y],{x,Union[C,D]},{y,Union[C,D]}],1]]]}]},
      If[m!={},
        Graph[Transpose[m][[1]], EdgeWeight->Transpose[m][[2]],
          EdgeLabels->MapThread[Rule,{Transpose[m][[1]],Transpose[m][[2]]}],
          EdgeStyle->MapThread[Rule,{Transpose[m][[1]],Transpose[m][[3]]}],
          VertexLabels->"Name"],
        Graph[{}]]]];

SameDigraphsQ[g1_Graph,g2_Graph] :=
  With[{order1=Flatten[Ordering[EdgeList[g1]]], order2=Flatten[Ordering[EdgeList[g2]]]},
    If[order1=={}||order2=={},
      If[order1=={}&&order2=={}, True, False],
      If[(EdgeList[g1][[order1]]==EdgeList[g2][[order2]]) &&
         (Last[Last[AbsoluteOptions[g1,EdgeWeight]]][[order1]] ==
          Last[Last[AbsoluteOptions[g2,EdgeWeight]]][[order2]]),
        True, False]]];

MoveVerticesQ[G_Graph,C_List,D_List] :=
  Module[{flag=If[D=={},True,False],
          subsets=SortBy[Cases[Subsets[D],x_/;x!={}],-Length[#]&], i=1},
    If[C=={}, False,
      While[!flag && i<=Length[subsets],
        If[SameDigraphsQ[Subdigraph[G,C,D],
            Subdigraph[G, Union[C,subsets[[i]]], Complement[D,subsets[[i]]]]],
          flag=True];
        i++];
      flag]];

DigraphUD[G_Graph,UDS_List,A_List,B_List] :=
  Module[{UD={}, U=If[UDS!={}, Last[UDS], {}]},
    If[U=={}&&A=={},
      UD=Join[UD,{{B}}];
      Table[UD=Union[UD, DigraphUD[G,{},x,Complement[B,x]]],
        {x, SortBy[Cases[Subsets[B],x_/;x!={}&&x!=B], -Length[#]&]}],
      If[!MoveVerticesQ[G,U,A],
        If[!MoveVerticesQ[G,A,B], UD=Join[UD,{Join[UDS,{A},{B}]}]];
        If[Length[B]>1,
          Table[UD=Union[UD, DigraphUD[G, Join[UDS,{A}], x, Complement[B,x]]],
            {x, SortBy[Cases[Subsets[B],x_/;x!={}&&x!=B], -Length[#]&]}]]]];
    UD];

UDSReps[intG_Graph] :=
  Module[{graphWeights=Flatten[Normal@WeightedAdjacencyMatrix[intG]]},
    If[MemberQ[graphWeights,"+"] || MemberQ[graphWeights,"-"],
      "Input has to be the interaction graph with no labels.",
      DigraphUD[intG, {}, {}, Sort[VertexList[intG]]]]];

(* ============================================================ *)
(* UPDATE SCHEDULE OPERATIONS *)
(* ============================================================ *)

ReflectedUDS[uds_List] :=
  With[{L=Max[Flatten[uds]]},
    Sort[#]& /@ (uds /. (MapThread[Rule,{Range[1,L], Range[L,1,-1]}]))];

RotatedUDS[uds_List,j_] :=
  With[{L=Max[Flatten[uds]]},
    Sort[#]& /@ (uds /. (MapThread[Rule,{Range[1,L], (1+Mod[#+j-1,L])& /@ Range[1,L]}]))];

RotatedUDS[uds_List] := Sort@Table[RotatedUDS[uds,j], {j,1,Max[Flatten[uds]]}];

CanonicalRotationUDSForm[uds_List] := First[Sort[RotatedUDS[uds]]];

RotationEquivUDSReps[intG_] := First /@ GatherBy[UDSReps[intG], CanonicalRotationUDSForm[#]&];

(* ============================================================ *)
(* ASYNCHRONISM DEGREE *)
(* ============================================================ *)

GetNeighbours[config_List,pos_Integer,r_] :=
  RotateLeft[config[[#]]& /@ (1+Mod[#,Length[config]]& /@
    Range[pos-\[LeftCeiling]r\[RightCeiling]-1, pos+\[LeftFloor]r\[RightFloor]-1]), \[LeftCeiling]r\[RightCeiling]];

VertexAsynchronismPattern[r_,UDS_List] :=
  Module[{depthlist=Table[If[MemberQ[First[UDS],i],1,0], {i,1,Length[Flatten[UDS]]}],
          len=Length[Flatten[UDS]], uds=UDS, aux},
    uds=Drop[uds,1];
    While[uds!={},
      aux=depthlist;
      Table[aux[[i]]=1+Total[depthlist[[GetNeighbours[Range[1,len],i,r]]]], {i,First[uds]}];
      depthlist=aux;
      uds=Drop[uds,1]];
    depthlist];

TotalAsynchronismDegree[r_,UDS_List,normalised_Symbol:False] :=
  With[{len=Length[Flatten[UDS]],
        max=Total[VertexAsynchronismPattern[r, Partition[Range[1,Length[Flatten[UDS]]],1]]-1]},
    If[\[Not]normalised,
      Total[VertexAsynchronismPattern[r,UDS]-1],
      Total[VertexAsynchronismPattern[r,UDS]-1]/max]];

(* ============================================================ *)
(* ASYNCHRONOUS CA AND BASINS OF ATTRACTION *)
(* ============================================================ *)

GetNeighbourPos[configLen_Integer,pos_Integer,r_] :=
  (1+Mod[#1,configLen]&) /@ Range[pos-\[LeftCeiling]r\[RightCeiling]-1, pos+\[LeftFloor]r\[RightFloor]-1];

fastAsynchronousCA[{rnum_Integer,k_Integer,r_},ic_List,tsteps_Integer,UDS_List] :=
  If[Mean[Length[#]& /@ UDS] < 0.03*Max[Flatten[UDS]],
    Block[{rule=Last[#]& /@ Reverse[RuleTable[rnum,k,r]],
           evo=Table[0,{i,1,tsteps+1},{j,1,Length[ic]}],
           neighbourslist=GetNeighbourPos[Length[ic],#,r]& /@ Range[1,Length[ic]]},
      evo[[1]]=ic;
      Do[evo[[i]]=evo[[i-1]];
        Do[evo[[i,u]]=Table[rule[[FromDigits[evo[[i,neighbourslist[[p]]]],k]+1]], {p,u}],
          {u,UDS}],
        {i,2,tsteps+1}];
      evo],
    Block[{evo=Table[0,{i,1,tsteps+1},{j,1,Length[ic]}],
           neighbourslist=GetNeighbourPos[Length[ic],#,r]& /@ Range[1,Length[ic]]},
      evo[[1]]=ic;
      Do[evo[[i]]=evo[[i-1]];
        Do[evo[[i,u]]=CellularAutomaton[{rnum,k,r},evo[[i]]][[u]], {u,UDS}],
        {i,2,tsteps+1}];
      evo]];

AsynchronousUdsBAField[{rnum_Integer,k_Integer:2,r_:1},uds_List] :=
  With[{n=Length[Flatten[uds]]},
    Graph[Table[DirectedEdge[ic, Last@fastAsynchronousCA[{rnum,k,r},ic,1,uds]],
      {ic,Tuples[Range[0,k-1],n]}]]];

(* ============================================================ *)
(* ANALYSIS FUNCTIONS *)
(* ============================================================ *)

getSortedBoA[regra_Integer,uds_List] :=
  Sort[WeaklyConnectedComponents[AsynchronousUdsBAField[{regra,2,1},uds]]];

AgruparUDSPorBaciasIdenticas[regra_Integer,udsLista_List] :=
  Module[{baciasComUDS,gruposEquivalentes},
    baciasComUDS=Map[{#,getSortedBoA[regra,#]}&, udsLista];
    gruposEquivalentes=GroupBy[baciasComUDS, Last->First];
    Values[gruposEquivalentes]];

ObterPadroesParticaoPorGrau[regra_Integer,listaUDSPorGrau_List] :=
  Module[{resultadosPorGrau={}},
    Scan[Module[{grau=#[[1]], udsDoGrau=#[[2]], grupos, tamanhos},
      grupos=AgruparUDSPorBaciasIdenticas[regra,udsDoGrau];
      tamanhos=Sort[Length /@ grupos, Greater];
      AppendTo[resultadosPorGrau,{grau,tamanhos}]]&,
      listaUDSPorGrau];
    SortBy[resultadosPorGrau,First]];

ProcessarTodasRegrasComParticoes[regrasLista_List,listaUDSPorGrau_List] :=
  Module[{resultados={}},
    Scan[Module[{regra=#, padroesGrau, assinatura},
      padroesGrau=ObterPadroesParticaoPorGrau[regra,listaUDSPorGrau];
      assinatura=Map[{#[[1]],#[[2]]}&, padroesGrau];
      AppendTo[resultados, regra->assinatura]]&,
      regrasLista];
    resultados];

AgruparRegrasPorAssinaturaParticao[resultadosRegras_List] :=
  Module[{agrupados},
    agrupados=GroupBy[resultadosRegras, Last->First];
    Normal[agrupados]];

(* ============================================================ *)
(* MAIN EXECUTION *)
(* ============================================================ *)

Print["Lattice size: ", latticeSize];
Print["Computing interaction graph..."];

intG = CAInteractionGraph[radius, latticeSize];

Print["Computing rotation-equivalent UDS representatives..."];
minReps = RotationEquivUDSReps[intG];
Print["Total UDS representatives: ", Length[minReps]];

Print["Computing ECA equivalence classes..."];
equivalenceECAClasses = First /@ GatherBy[Range[0,255], DynamicalEquivalenceClass[#]&];
Print["Total ECA classes: ", Length[equivalenceECAClasses]];

Print["Grouping UDSs by asynchronism degree..."];
totalLista = Normal[GroupBy[minReps, TotalAsynchronismDegree[radius,#]&]];
listaGruposAssincronos = Rest[totalLista];

Print["Processing all rules with partition signatures..."];
resultadosParticoes = ProcessarTodasRegrasComParticoes[equivalenceECAClasses, listaGruposAssincronos];
regrasAgrupadasPorParticao = AgruparRegrasPorAssinaturaParticao[resultadosParticoes];
regrasAgrupadasOrdenadas = SortBy[regrasAgrupadasPorParticao, {-Length[#[[2]]]&, First[Sort[#[[2]]]]&}];

Print["\n(Format: Degree -> (sizes of groups with identical basins))\n"];

gridDataParticoes = Prepend[
  Map[{Column[Map[Row[{#[[1]], " \[RightArrow] ", ToString[#[[2]]]}]&, #[[1]]],
    Spacings->0.2, Alignment->Left], ToString[Sort[#[[2]]]]}&,
    regrasAgrupadasOrdenadas],
  {Style["Partition Signature", Bold], Style["ECA Rules", Bold]}];

Grid[gridDataParticoes, Frame->All, Alignment->{Left,Top},
  Spacings->{1,0.8}, Background->{None,{LightGray,White}}]
