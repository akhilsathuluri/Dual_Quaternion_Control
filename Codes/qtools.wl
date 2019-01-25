(* ::Package:: *)

(*This package contains the modules required for dealing with quaternion*)


BeginPackage["qtools`","rotationtools`"]


quat::usage="Given angle and an axis returns the corresponding unit quaternion"
qmult::usage="Multiplies two input quaternions"
cquat::usage="Returns the conjugate of a given quaternion"
unitdq::usage="Returns a unit dual quaternion (of translation followed by rotation) given translation, angle, and axis"
cdualq::usage="Returns the conjugate of a dual quaternion"
dqmult::usage="Multiplies two input dual quaternion"
Hp::usage="Gives the corresponding H+ Hamilton operator of the quaternion"
Hn::usage="Gives the corresponding H- Hamilton operator of the quaternion"
dHp::usage="Gives the corresponding H+ Hamilton operator of the quaternion"
dHn::usage="Gives the corresponding H- Hamilton operator of the dual quaternion"
decq::usage="Decomposes a given quaternion into angle,  axis form"
logdq::usage="Gives the logarithm (se(3)) element of a dual quaternion (SE(3))"
explogdq::usage="Finds the exponent of a se(3) element and gives an element from SE(3)"
R82dual::usage="Takes an element from \!\(\*SuperscriptBox[\(R\), \(\(8\)\(\\\ \)\)]\)and returns a dual quaternion"
R62dual::usage="Takes an element from \!\(\*SuperscriptBox[\(R\), \(\(6\)\(\\\ \)\)]\)and returns the corresponfing se(3) element of a dual quaternion"
frameplotter::usage="Takes in an se(3) element corresponding to a dual quaternion and plots its frame with respecto the origin"


Begin["`Private`"]; 


(*Quaternion Modules*)
quat[\[Theta]_,k_ ]:=Module[{},
Return[{Cos[\[Theta]/2], Sin[\[Theta]/2]{k[[1]], k[[2]], k[[3]]}}];
]
qmult[q1_, q2_]:=Module[{},
Return[{q1[[1]]*q2[[1]]-q1[[2]].q2[[2]], q1[[1]]*q2[[2]]+q1[[2]]*q2[[1]]+Cross[q1[[2]], q2[[2]]]}];
]
cquat[q_]:=Module[{},
Return[{q[[1]], -q[[2]]}];
]


(*Dual Quaternion Modules*)
(*Module for unit dual quaternions, given a unit rotation \[Theta], axis, k and translation t. This module is to produce a dual quaternion representing a translation followed by rotation*)
unitdq[ t_, \[Theta]_, k_]:=Module[{qrot},
qrot = quat[\[Theta], k];
Return[{qrot, qmult[{0,t}, qrot]/2}];
]
cdualq[dq1_]:=Module[{},
Return[{cquat[dq1[[1]]], cquat[dq1[[2]]]}];
];
dqmult[dq1_, dq2_]:=Module[{},
Return[ {qmult[dq1[[1]], dq2[[1]]], qmult[dq1[[2]], dq2[[1]]]+qmult[dq1[[1]], dq2[[2]]]}];
]


(*Testing the Hamiltonian Operator for given quaternion or a dual quaternion*)
Hp[q_] := Module[{a0, a1, a2, a3},
a0 = q[[1]];a1 = q[[2]][[1]];a2 = q[[2]][[2]];a3 = q[[2]][[3]];
Return[{{a0, -a1, -a2, -a3},{a1, a0, -a3, a2},{a2, a3, a0, -a1},{a3, -a2, a1, a0}}];
]
Hn[q_] := Module[{a0, a1, a2, a3},
a0 = q[[1]];a1 = q[[2]][[1]];a2 = q[[2]][[2]];a3 = q[[2]][[3]];
Return[{{a0, -a1, -a2, -a3},{a1, a0, a3, -a2},{a2, -a3, a0, a1},{a3, a2, -a1, a0}}];
]


(*Hailtonian Operators extended to dual quaternion*)
dHp[dq1_]:=Module[{},
Return[N[ArrayFlatten[{{Hp[dq1[[1]]], ConstantArray[0,{4,4}]},{Hp[dq1[[2]]], Hp[dq1[[1]]]}}]]];
]
dHn[dq1_]:=Module[{},
Return[N[ArrayFlatten[{{Hn[dq1[[1]]], ConstantArray[0,{4,4}]},{Hn[dq1[[2]]], Hn[dq1[[1]]]}}]]];
]


(*Module to extract the orientation from a quaternion*)
decq[q_]:=Module[{sine, \[Phi],k},
sine = Sqrt[q[[2]].q[[2]]];
\[Phi] = 2*ArcTan[q[[1]], sine];
k = q[[2]]/Sin[\[Phi]/2];
Return[{\[Phi], k}];
]


(*To implement this, we need to extract the position and angle at each instant from the FK quaternion*)
logdq[dq1_]:=Module[{tx, ty, tz, solt, decqv, logqv},
decqv = decq[dq1[[1]]];
solt =Solve[(N[(qmult[{0,{tx, ty, tz}}, dq1[[1]]]/2)-dq1[[2]]][[2]])==0, {tx, ty, tz}][[1]];
logqv = {decqv[[1]]*decqv[[2]],{tx, ty, tz}/.solt}/2;
Return[logqv];
]
(*Hence the normal logarithmic identities do not work for dual quaternions as the multiplication here has different rules*)
explogdq[logdq_]:=Module[{temp\[Theta], tempk},
temp\[Theta] = Norm[logdq[[1]]];
tempk = logdq[[1]]/temp\[Theta];
Return[unitdq[2*logdq[[2]], 2*temp\[Theta], tempk]];
];


(*Since the integrated error will be in R^8 we need a way to map it back to dual quaternions*)
R82dual[rq_]:=Module[{},
Return[{{rq[[1]],rq[[2;;4]]},{rq[[5]], rq[[6;;8]]}}];
]

R62dual[z_]:=Module[{},
Return[{z[[1;;3]],z[[4;;6]]}];
]


(*Plots the frame given the log of a dual quaternion or the se(3) element*)
frameplotter[logdq_]:=Module[{\[Theta]vec, \[Theta]val, kvec, Amat, Rmat, tvec, origin, l, point, r},
(*l is the length of the axis that needs to be plotted and the radius of the point*)
l = 0.2;
r = 0.01;
(*Show[Graphics3D[{Blue, Arrow[{{0,0,0},{0,0,1}}]}], Graphics3D[{Green, Arrow[{{0,0,0},{0,1,0}}]}], Graphics3D[{Red, Arrow[{{0,0,0},{1,0,0}}]}]]*)
\[Theta]vec = 2*logdq[[1;;3]];
tvec = 2*logdq[[4;;6]];
\[Theta]val = Norm[\[Theta]vec];
kvec = \[Theta]vec/\[Theta]val;
(*Using the Rodrigues formulation for finding the rotation matrix. So formulating the Rodrigues parameters*)
Amat = vec2SkewMat[kvec*Tan[\[Theta]val/2]];
Rmat = (Inverse[IdentityMatrix[3]-Amat].(IdentityMatrix[3]+Amat));
(*Plots z, y, x, axis respectively*)
origin = Show[Graphics3D[{Black, Sphere[{0,0,0}, r]}],Graphics3D[{Blue, Arrow[{{0,0,0},{0,0,l}}]}], Graphics3D[{Green, Arrow[{{0,0,0},{0,l,0}}]}], Graphics3D[{Red, Arrow[{{0,0,0},{l,0,0}}]}]];
(*We also need to mark the coordinate to which the frame is attached*)
point = Graphics3D[{Black, Sphere[tvec,r]}];
(*Returning the transformed frame and the origin*)
Return[Show[origin, Graphics3D[{Blue, Arrow[{tvec,tvec+Rmat.{0,0,l}}]}], Graphics3D[{Green, Arrow[{tvec,tvec+Rmat.{0,l,0}}]}], Graphics3D[{Red, Arrow[{tvec,tvec+Rmat.{l,0,0}}]}], point]];
]


End[]
EndPackage[]
