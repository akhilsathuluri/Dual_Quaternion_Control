(* ::Package:: *)

BeginPackage["rotationtools`"]


rotationtools::usage = " 'rotationtools' package povides the following functions"
rotX::usage = "Gives the rotation matrix with respect to X"
rotY::usage = "Gives the rotation matrix with respect to Y"
rotZ::usage = "Gives the rotation matrix with respect to Z"
axisToRotMat::usage = "Gives the rotation matrix about the axis"
rotMatToAxis::usage = "Gives the axis and rotation angle given a rotation matrix"
solveEulerZYX::usage= "Solves for individual rotations \[Theta]1, \[Theta]2, \[Theta]3 given ZYX rotation"
solveEulerZYZ::usage = "Solves for individual angles \[Theta]1, \[Theta]2, \[Theta]3 given a ZYZ rotation"
SkewMat2vec::usage = "Given a Skew symmetric matrix, this function returns the corresponding vector"
vec2SkewMat::usage = "Given a vec it gives the skew mat"
omegaJacobianBodyFixed::usage = "Given a rotation matrix R and the corresponding angles \[Theta], it returns the angular velocity Jacobian of the manipulator in the body fixed frame"
omegaJacobianSpaceFixed::usage = "Given a rotation matrix R and the corresponding angles \[Theta], it returns the angular velocity Jacobian of the manipulator in the space fixed frame"
Hom::usage = "Returns a Homogenous transformation matrix given rotation and translation"
DH2HomTrans::usage = "Returns a homogenous transformation matrix between the i-1th and ith link"
multiplyHomTrans::usage = "Multiplies two homogenous transformation matrices to reduce computational time"
DHTable2HomTrans::usage = "Returns the end effector orientation given a table of DH parameters"
omegaJacobianEndEffector::usage = "Returns the angular velovity Jaconian of the end effector given a homogeneous transformation matrix"
linearJacobianEndEffector::usage = "Returns the velocity Jacobian of the endeffector given a homogeneous transformation matrix"
MMat2CMat::usage = "Module to calculate the C matrix given a mass matrix M"
IKPUMA::usage = "Completes the inverse kinematics of the PUMA 560 robot and returns the all 8 possible orientations"


Begin["`Private`"]; 


rotX[\[Theta]_]:=
	Module[{},
		Return[{{1, 0, 0}, {0, Cos[\[Theta]], -Sin[\[Theta]]},{0, Sin[\[Theta]], Cos[\[Theta]]}}];
	];


rotY[\[Theta]_]:=
	Module[{roty},
		roty = {{Cos[\[Theta]], 0, Sin[\[Theta]]}, {0, 1, 0}, {-Sin[\[Theta]], 0, Cos[\[Theta]]}};
		Return[roty];
	];


rotZ[\[Theta]_]:=
	Module[{rotz},
		rotz = {{Cos[\[Theta]], -Sin[\[Theta]], 0}, {Sin[\[Theta]], Cos[\[Theta]], 0}, {0, 0, 1}};
		Return[rotz];
	];


axisToRotMat[{k1_, k2_, k3_}, \[Phi]_]:=
	Module[{K, Rk\[Phi], Iden, k1s, k2s, k3s},
	k1s = k1;
	k2s = k2;
	k3s = k3;
	K = {{0, -k3s, k2s}, {k3s, 0, -k1s},{-k2s, k1s, 0}};
	(*K = {{0, k3, -k2}, {-k3, 0, k1},{k2, -k1, 0}};*)
	Iden = {{1,0,0},{0,1,0},{0,0,1}};	
	Rk\[Phi] = Iden+(K)*Sin[\[Phi]]+(K.K)*(1-Cos[\[Phi]]);
	Return[Rk\[Phi]];
	];


rotMatToAxis[Rk\[Phi]_]:=
Module[{r11, r12, r13, r21, r22, r23, r31, r32, r33, Cos\[Phi], Sin\[Phi], \[Phi], k1, k2, k3},
		r11 = Rk\[Phi][[1]][[1]];
		r12 = Rk\[Phi][[1]][[2]];
		r13 = Rk\[Phi][[1]][[3]];
		r21 = Rk\[Phi][[2]][[1]];
		r22 = Rk\[Phi][[2]][[2]];
		r23 = Rk\[Phi][[2]][[3]];
		r31 = Rk\[Phi][[3]][[1]];
		r32 = Rk\[Phi][[3]][[2]];
		r33 = Rk\[Phi][[3]][[3]];
		Cos\[Phi] = (r11+r22+r33-1)/2;
		(*Sin\[Phi] = {Sqrt[((r12-r21)/2)^2+((r31-r13)/2)^2+((r23-r32)/2)^2], -Sqrt[((r12-r21)/2)^2+((r31-r13)/2)^2+((r23-r32)/2)^2]};*)
		\[Phi] = ArcCos[Cos\[Phi]];

		k3 = (r21-r12)/(2*Sin[\[Phi]]);
		k2 = (r13-r31)/(2*Sin[\[Phi]]);
		k1 = (r32-r23)/(2*Sin[\[Phi]]);
		Return[{k1, k2, k3, \[Phi]}];
		];
(*	k3 = (r12-r21)/(2*Sin\[Phi]);
		k2 = (r31-r13)/(2*Sin\[Phi]);
		k1 = (r23-r32)/(2*Sin\[Phi]);
		Return[{k1, k2, k3, \[Phi]}];*)
		(*If[\[Phi]\[Equal]0,
		Return[Print[StringForm["Error: There is no rotation, hence cannot find the rotation axis."]]]]
		];*)


solveEulerZYX[R_]:=
	Module[{r11, r12, r13, r21, r22, r23, r31, r32, r33, \[Theta]1, \[Theta]2, \[Theta]3},
		r11 = R[[1]][[1]];
		r12 = R[[1]][[2]];
		r13 = R[[1]][[3]];
		r21 = R[[2]][[1]];
		r22 = R[[2]][[2]];
		r23 = R[[2]][[3]];
		r31 = R[[3]][[1]];
		r32 = R[[3]][[2]];
		r33 = R[[3]][[3]];

		If[r13 !=  +1  && r13 != -1 , 
			\[Theta]2 = {ArcTan[+Sqrt[r11^2+r12^2], r13], ArcTan[-Sqrt[r11^2+r12^2], r13]}; 
			\[Theta]1 = ArcTan[r33/Cos[\[Theta]2], -r23/Cos[\[Theta]2]];
			\[Theta]3 = ArcTan[r11/Cos[\[Theta]2], -r12/Cos[\[Theta]2]];
			Return[{\[Theta]1, \[Theta]2, \[Theta]3}];
			];
		If[r13 == 1, 
			\[Theta]1 = ArcTan[r22, r21];
			\[Theta]2 = \[Pi]/2;
			\[Theta]3 = 0;
			Return[{\[Theta]1, \[Theta]2, \[Theta]3}];
			];
		If[r13 == -1, 
			\[Theta]1 = -ArcTan[r22, r21];
			\[Theta]2 = -\[Pi]/2;
			\[Theta]3 = 0;
			Return[{\[Theta]1, \[Theta]2, \[Theta]3}];
			];
	];


solveEulerZYZ[R_]:=
	Module[{r11, r12, r13, r21, r22, r23, r31, r32, r33, \[Theta]1, \[Theta]2, \[Theta]3},
		r11 = R[[1]][[1]];
		r12 = R[[1]][[2]];
		r13 = R[[1]][[3]];
		r21 = R[[2]][[1]];
		r22 = R[[2]][[2]];
		r23 = R[[2]][[3]];
		r31 = R[[3]][[1]];
		r32 = R[[3]][[2]];
		r33 = R[[3]][[3]];
		If[r33 != 1||r33!=-1, 
			\[Theta]2 = {ArcTan[r33, Sqrt[r31^2+r32^2]],ArcTan[r33, -Sqrt[r31^2+r32^2]]};
			\[Theta]1 = ArcTan[r13/Sin[\[Theta]2], r23/Sin[\[Theta]2]];
			\[Theta]3 = ArcTan[-r31/Sin[\[Theta]2], r32/Sin[\[Theta]2]];
			Return[{\[Theta]1, \[Theta]2, \[Theta]3}];
			];
		If[r33 ==  1, 
			\[Theta]1 = 0;
			\[Theta]2 = 0;
			\[Theta]3 = ArcTan[r11, -r12];
			Return[{\[Theta]1, \[Theta]2, \[Theta]3}];
			];
		If[r13 ==  -1, 
			\[Theta]1 = 0;
			\[Theta]2 = \[Pi];
			\[Theta]3 = ArcTan[-r11, r12];
			Return[{\[Theta]1, \[Theta]2, \[Theta]3}];
			];
	];


(*The following functions deal with the homogenous transformation matrices and Jacobian matrices*)


SkewMat2vec[R_]:=
Module[{},
Return[{R[[3, 2]], R[[1, 3]], R[[2, 1]]}];
];


vec2SkewMat[v_]:=
Module[{},
Return[{{0, -v[[3]], v[[2]]}, {v[[3]], 0, -v[[1]]},{-v[[2]], v[[1]], 0}}];
] ;


(*SkewMat2vec[R_]:=
Module[{r11, r12, r13, r21, r22, r23, r31, r32, r33},
		r11 = R[[1]][[1]];
		r12 = R[[1]][[2]];
		r13 = R[[1]][[3]];
		r21 = R[[2]][[1]];
		r22 = R[[2]][[2]];
		r23 = R[[2]][[3]];
		r31 = R[[3]][[1]];
		r32 = R[[3]][[2]];
		r33 = R[[3]][[3]];
(*If[r11==0&&r22==0&&r33==0&&r21==-r12&&r13==-r31&&r23==-r32,
Return[{r32, r13, r21}],
Print["The entered matrix is not a Skew Symmetric Matrix"]*)
Return[{r32, r13, r21}];
] ;*)


omegaJacobianBodyFixed[R_, \[Theta]_]:=
Module[{i, n, Jw},
n = Dimensions[\[Theta]][[1]];
Jw = ConstantArray[0, {n, 3}];
For[i = 1, i<=n, i = i+1, 
Jw[[i]] = SkewMat2vec[Simplify[Transpose[R].D[R, \[Theta][[i]]]]];
];
Return[Transpose[Jw]];
];


omegaJacobianSpaceFixed[R_, \[Theta]_]:=
			Module[{n, i, Rtemp, Jwb},
			n = Dimensions[R][[1]];
			Jwb = omegaJacobianBodyFixed[R, \[Theta]];
			Rtemp = R[[1]];
			For[i = 2, i<=n, i=i+1, 
			Rtemp = Rtemp.R[[i]];
				];
			Return[Rtemp.Jwb.Transpose[Rtemp]];
			];


Hom[R_, a_]:=
Module[{H, Zer},
Zer = {0, 0, 0, 1};
H = {{R[[1]][[1]], R[[1]][[2]], R[[1]][[3]], a[[1]]}, {R[[2]][[1]], R[[2]][[2]], R[[2]][[3]], a[[2]]}, {R[[3]][[1]], R[[3]][[2]], R[[3]][[3]], a[[3]]}, Zer};
Return[H];
];


DH2HomTrans[aim1_, \[Alpha]im1_, di_, \[Theta]i_]:=
Module[{Trans, Iden, Zer, R, Twist, a, Offset, LinkR, HM},
(*First we move Zi-1 to Zi. So only translation*)
Iden = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
Zer = {0, 0, 0, 1};
R = Iden;
a = {aim1, 0, 0};
Trans = {{R[[1]][[1]], R[[1]][[2]], R[[1]][[3]], a[[1]]}, {R[[2]][[1]], R[[2]][[2]], R[[2]][[3]], a[[2]]}, {R[[3]][[1]], R[[3]][[2]], R[[3]][[3]], a[[3]]}, Zer};
(*Now the translated Zi-1 axis will twist an angle of twist, which is give by rotation matrix Rot*)
R = {{1, 0, 0}, {0, Cos[\[Alpha]im1], -Sin[\[Alpha]im1]}, {0, Sin[\[Alpha]im1], Cos[\[Alpha]im1]}};
a = {0, 0, 0};
Twist = {{R[[1]][[1]], R[[1]][[2]], R[[1]][[3]], a[[1]]}, {R[[2]][[1]], R[[2]][[2]], R[[2]][[3]], a[[2]]}, {R[[3]][[1]], R[[3]][[2]], R[[3]][[3]], a[[3]]}, Zer};
(*Now it oriented same as Zi, but has an offset of di*)
R = Iden;
a = {0, 0, di};
(*Since we will need to describe di in our coordinate frame*)
Offset = {{R[[1]][[1]], R[[1]][[2]], R[[1]][[3]], a[[1]]}, {R[[2]][[1]], R[[2]][[2]], R[[2]][[3]], a[[2]]}, {R[[3]][[1]], R[[3]][[2]], R[[3]][[3]], a[[3]]}, Zer};
(*Now its the link rotation*)
R = {{Cos[\[Theta]i], -Sin[\[Theta]i], 0}, {Sin[\[Theta]i], Cos[\[Theta]i], 0}, {0, 0, 1}};
a = {0, 0, 0};
LinkR = {{R[[1]][[1]], R[[1]][[2]], R[[1]][[3]], a[[1]]}, {R[[2]][[1]], R[[2]][[2]], R[[2]][[3]], a[[2]]}, {R[[3]][[1]], R[[3]][[2]], R[[3]][[3]], a[[3]]}, Zer};
HM = Trans.Twist.Offset.LinkR;
Return[HM];
];


(*We write another code to multiply the homogenous transformation matrices as using direct multiplication might cause increase in time taken for computation*)
multiplyHomTrans[T_]:=
Module[{R, R1, R2, p1, p2, n, MulT, p, Linear, Rot, i},
n = Dimensions[T][[1]];
R = ConstantArray[0, {n, 3, 3}];
MulT = ConstantArray[0, {4, 4}];
p = ConstantArray[0, {2, 3}];
For[i = 1, i<=n, i = i+1,
R[[i]][[1]][[1]] = T[[i]][[1]][[1]];
R[[i]][[1]][[2]] = T[[i]][[1]][[2]];
R[[i]][[1]][[3]] = T[[i]][[1]][[3]];
R[[i]][[2]][[1]] = T[[i]][[2]][[1]];
R[[i]][[2]][[2]] = T[[i]][[2]][[2]];
R[[i]][[2]][[3]] = T[[i]][[2]][[3]];
R[[i]][[3]][[1]] = T[[i]][[3]][[1]];
R[[i]][[3]][[2]] = T[[i]][[3]][[2]];
R[[i]][[3]][[3]] = T[[i]][[3]][[3]];
p[[i]][[1]] = T[[i]][[1]][[4]];
p[[i]][[2]] = T[[i]][[2]][[4]];
p[[i]][[3]] = T[[i]][[3]][[4]];
];
Linear = R[[1]].p[[2]]+p[[1]];
Rot = R[[1]].R[[2]];
MulT[[1]][[1]] = Rot[[1]][[1]];
MulT[[1]][[2]] = Rot[[1]][[2]];
MulT[[1]][[3]] = Rot[[1]][[3]];
MulT[[1]][[4]] = Linear[[1]];
MulT[[2]][[4]] = Linear[[2]];
MulT[[3]][[4]] = Linear[[3]];
MulT[[2]][[1]] = Rot[[2]][[1]];
MulT[[2]][[2]] = Rot[[2]][[2]];
MulT[[2]][[3]] = Rot[[2]][[3]];
MulT[[3]][[1]] = Rot[[3]][[1]];
MulT[[3]][[2]] = Rot[[3]][[2]];
MulT[[3]][[3]] = Rot[[3]][[3]];
MulT[[4]][[1]] = 0;
MulT[[4]][[2]] = 0;
MulT[[4]][[3]] = 0;
MulT[[4]][[4]] = 1;
Return[MulT];
];


(*The table is expected to have rows corresponding to each link and coulmns with parameters ordered as ai-1, \[Alpha]i-1, di, \[Theta]i*)
DHTable2HomTrans[DHtab_]:=
Module[{Dim, Hn, n, H1, R1, \[Theta]1, a1, an, Rn, Iden, i, Htemp},
Dim = Dimensions[DHtab];
Iden = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
n = Dim[[1]];
\[Theta]1 = DHtab[[1]][[4]];
R1 = rotZ[\[Theta]1];
a1 = {0, 0, 0};
H1 = Hom[R1, a1];
Htemp = H1;
For[i = 2, i <= n, i = i+1, 
Htemp = Htemp.DH2HomTrans[DHtab[[i]][[1]], DHtab[[i]][[2]], DHtab[[i]][[3]], DHtab[[i]][[4]]];
];
Return[Htemp];
];


omegaJacobianEndEffector[T_, \[Theta]_]:=
Module[{R, i},
R = ConstantArray[0, {3, 3}];
R[[1]][[1]] = T[[1]][[1]];
R[[1]][[2]] = T[[1]][[2]];
R[[1]][[3]] = T[[1]][[3]];
R[[2]][[1]] = T[[2]][[1]];
R[[2]][[2]] = T[[2]][[2]];
R[[2]][[3]] = T[[2]][[3]];
R[[3]][[1]] = T[[3]][[1]];
R[[3]][[2]] = T[[3]][[2]];
R[[3]][[3]] = T[[3]][[3]];
Return[omegaJacobianBodyFixed[R, \[Theta]]];
];


linearJacobianEndEffector[T_, \[Theta]_]:=
Module[{ax, ay, az, p, n, i, Jv},
n = Dimensions[\[Theta]][[1]];
Jv = ConstantArray[0,{n,3}];
ax = Simplify[T[[1]][[4]]];
ay = Simplify[T[[2]][[4]]];
az = Simplify[T[[3]][[4]]];
p = {ax, ay, az};
For[i = 1, i<=n, i = i+1,
Jv[[i]] = D[p, \[Theta][[i]]];
];
Return[Transpose[Jv]];
];


(*Backup*)

(*DHTable2HomTrans[DHtab_]:=
Module[{Dim, Hn, n, H1, R1, \[Theta]1, a1, an, Rn, Iden, i, Htemp},
Dim = Dimensions[DHtab];
Iden = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
n = Dim[[1]];
\[Theta]1 = DHtab[[1]][[4]];
R1 = rotZ[\[Theta]1];
a1 = {0, 0, 0};
H1 = Hom[R1, a1];
Htemp = H1;
For[i = 2, i <= n, i = i+1, 
Htemp = Htemp.DH2HomTrans[DHtab[[i]][[1]], DHtab[[i]][[2]], DHtab[[i]][[3]], DHtab[[i]][[4]]];
];
Return[Htemp.Hom[Iden,{DHtab[[n]][[1]], 0, 0}]];
];*)


(*Module to calculate the C matrix given a mass matrix M*)
MMat2CMat[M_, q_] :=
Module[{C, n, i, j, k, tempM, ii, jj, kk, tempC},
(*Compute the derivatives of the upper triangular elements of the mass matrix*)
(*Computational advantage is achived by computing the dericatives of only the upper triangle elements of the mass matrix and storing them and 
using them later for the computation of the C matrix unlike in the three loop case where the same derivative is computed again and again here they are already stored and are just reused any number of 
times reducing time for recomputation*)
n = Length[q];
tempM = ConstantArray[0, {Dimensions[M][[1]],Dimensions[M][[2]],n}];
tempC = ConstantArray[0, {Dimensions[M][[1]],Dimensions[M][[2]],n}];
For[i=1,i<n,
For[j=i, j<n, 
For[k=1, k<n,
tempM[[i,j,k]] = D[M[[i,j]],q[[k]]];
tempM[[j,i,k]] = tempM[[j,i,k]];
k++];
j++];
i++];
For [ii = 1, ii <= n, 
For[jj = 1, jj <= n,
For[kk = 1, kk<= n, 
tempC[[ii,jj]] +=1/2*( tempM[[i, j , k]] + tempM[[i, k, j]] - tempM[[j, k , i]])*D[q[[kk]],t] ;
kk++];
jj++];
ii++];
Return[tempC];
];


IKPUMA[R_,p_]:=
Module[{x, y, z, \[Theta]1, \[Theta]3, \[Theta]2, \[Theta]4, \[Theta]5, \[Theta]6, a2, a3, d3, d4, K, Rc1, Rc2, Rc3, Rc4, dh1, dh2, dh3, Hc1, Hc2, Hc3, Hc4, \[Theta]51, \[Theta]52, \[Theta]53, \[Theta]54, \[Theta]41, \[Theta]42, \[Theta]43, \[Theta]44, \[Theta]61, \[Theta]62, \[Theta]63, \[Theta]64},
a2 = 0.4318;
a3 = 0.019;
d3 = 0.125;
d4 = 0.432;
x = p[[1]];
y = p[[2]];
z = p[[3]];
K = (1/(2*a2))*(x^2+y^2+z^2-d3^2-a2^2-a3^2-d4^2);


\[Theta]1 = {2*ArcTan[(-x+Sqrt[x^2+y^2-d3^2])/(y+d3)],2*ArcTan[(-x-Sqrt[x^2+y^2-d3^2])/(y+d3)]};
\[Theta]3 = {2*ArcTan[(-d4+Sqrt[d4^2+a3^2-K^2])/(K+a3)],2*ArcTan[(-d4-Sqrt[d4^2+a3^2-K^2])/(K+a3)]};
\[Theta]2 ={2*ArcTan[(z-a3*Sin[\[Theta]3[[1]]]-d4*Cos[\[Theta]3[[1]]]),(-a2-a3*Cos[\[Theta]3[[1]]]+d4*Sin[\[Theta]3[[1]]]+Sqrt[a2^2+a3^2+d4^2+2*a2*(a3*Cos[\[Theta]3[[1]]]-d4*Sin[\[Theta]3[[1]]])-z^2])],2*ArcTan[(z-a3*Sin[\[Theta]3[[2]]]-d4*Cos[\[Theta]3[[2]]]),(-a2-a3*Cos[\[Theta]3[[2]]]+d4*Sin[\[Theta]3[[2]]]+Sqrt[a2^2+a3^2+d4^2+2*a2*(a3*Cos[\[Theta]3[[2]]]-d4*Sin[\[Theta]3[[2]]])-z^2])],2*ArcTan[(z-a3*Sin[\[Theta]3[[1]]]-d4*Cos[\[Theta]3[[1]]]),(-a2-a3*Cos[\[Theta]3[[1]]]+d4*Sin[\[Theta]3[[1]]]-Sqrt[a2^2+a3^2+d4^2+2*a2*(a3*Cos[\[Theta]3[[1]]]-d4*Sin[\[Theta]3[[1]]])-z^2])],2*ArcTan[(z-a3*Sin[\[Theta]3[[2]]]-d4*Cos[\[Theta]3[[2]]]),(-a2-a3*Cos[\[Theta]3[[2]]]+d4*Sin[\[Theta]3[[2]]]-Sqrt[a2^2+a3^2+d4^2+2*a2*(a3*Cos[\[Theta]3[[2]]]-d4*Sin[\[Theta]3[[2]]])-z^2])]};

(*For case 1*)
dh1 = DH2HomTrans[0,0,0,\[Theta]1[[2]]];
dh2 = DH2HomTrans[0,-\[Pi]/2,0,\[Theta]2[[4]]];
dh3 = DH2HomTrans[a2,0,d3,\[Theta]3[[2]]];

Hc1 = dh1.dh2.dh3;
Rc1 = Transpose[Hc1[[1;;3,1;;3]]].R;

If[Rc1[[2,3]]!= 1&&Rc1[[2,3]]!= 1,
\[Theta]51 = {ArcTan[Rc1[[2,3]],Sqrt[Rc1[[2,1]]^2+Rc1[[2,2]]^2]],ArcTan[Rc1[[2,3]],-Sqrt[Rc1[[2,1]]^2+Rc1[[2,2]]^2]]};
\[Theta]41 = ArcTan[-Rc1[[1,3]]/Sin[\[Theta]51],Rc1[[3,3]]/Sin[\[Theta]51]];
\[Theta]61 = ArcTan[Rc1[[2,1]]/Sin[\[Theta]51],-Rc1[[2, 2]]/Sin[\[Theta]51]];
];
If[Rc1[[2,3]] == 1,
\[Theta]41 = {0,0};
\[Theta]51 = {0,0};
\[Theta]61 = ArcTan[Rc1[[1,1]],-Rc1[[1,2]]];
];
If[Rc1[[2,3]] == -1,
\[Theta]41 = {0,0};
\[Theta]51 = {\[Pi], \[Pi]};
\[Theta]61 = {-ArcTan[-Rc1[[1,1]],Rc1[[1,2]]],-ArcTan[-Rc1[[1,1]],Rc1[[1,2]]]};
];

(*For case 2*)
dh1 = DH2HomTrans[0,0,0,\[Theta]1[[1]]];
dh2 = DH2HomTrans[0,-\[Pi]/2,0,\[Theta]2[[2]]];
dh3 = DH2HomTrans[a2,0,d3,\[Theta]3[[2]]];

Hc2 = dh1.dh2.dh3;
Rc2 = Transpose[Hc2[[1;;3,1;;3]]].R;

If[Rc2[[2,3]]!= 1&&Rc2[[2,3]]!= 1,
\[Theta]52 = {ArcTan[Rc2[[2,3]],Sqrt[Rc2[[2,1]]^2+Rc2[[2,2]]^2]],ArcTan[Rc2[[2,3]],-Sqrt[Rc2[[2,1]]^2+Rc2[[2,2]]^2]]};
\[Theta]42 = ArcTan[-Rc2[[1,3]]/Sin[\[Theta]52],Rc2[[3,3]]/Sin[\[Theta]52]];
\[Theta]62 = ArcTan[Rc2[[2,1]]/Sin[\[Theta]52],-Rc2[[2, 2]]/Sin[\[Theta]52]];
];
If[Rc2[[2,3]] == 1,
\[Theta]42 = {0,0};
\[Theta]52 = {0,0};
\[Theta]62 = {ArcTan[Rc2[[1,1]],-Rc2[[1,2]]],ArcTan[Rc2[[1,1]],-Rc2[[1,2]]]};
];
If[Rc2[[2,3]] == -1,
\[Theta]42 = {0,0};
\[Theta]52 = {\[Pi],\[Pi]};
\[Theta]62 = {-ArcTan[-Rc2[[1,1]],Rc2[[1,2]]],-ArcTan[-Rc2[[1,1]],Rc2[[1,2]]]};
];

(*For case 3*)
dh1 = DH2HomTrans[0,0,0,\[Theta]1[[2]]];
dh2 = DH2HomTrans[0,-\[Pi]/2,0,\[Theta]2[[3]]];
dh3 = DH2HomTrans[a2,0,d3,\[Theta]3[[1]]];

Hc3 = dh1.dh2.dh3;
Rc3 = Transpose[Hc3[[1;;3,1;;3]]].R;

If[Rc3[[2,3]]!= 1&&Rc3[[2,3]]!= 1,
\[Theta]53 = {ArcTan[Rc3[[2,3]],Sqrt[Rc3[[2,1]]^2+Rc3[[2,2]]^2]],ArcTan[Rc3[[2,3]],-Sqrt[Rc3[[2,1]]^2+Rc3[[2,2]]^2]]};
\[Theta]43 = ArcTan[-Rc3[[1,3]]/Sin[\[Theta]53],Rc3[[3,3]]/Sin[\[Theta]53]];
\[Theta]63 = ArcTan[Rc3[[2,1]]/Sin[\[Theta]53],-Rc3[[2, 2]]/Sin[\[Theta]53]];
];
If[Rc3[[2,3]] ==  1,
\[Theta]43 = {0,0};
\[Theta]53 = {0,0};
\[Theta]63 = {ArcTan[Rc3[[1,1]],-Rc3[[1,2]]],ArcTan[Rc3[[1,1]],-Rc3[[1,2]]]};
];
If[Rc3[[2,3]] == -1,
\[Theta]43 = {0,0};
\[Theta]53 = {\[Pi], \[Pi]};
\[Theta]63 = {-ArcTan[-Rc3[[1,1]],Rc3[[1,2]]],-ArcTan[-Rc3[[1,1]],Rc3[[1,2]]]};
];

(*For case 4*)
dh1 = DH2HomTrans[0,0,0,\[Theta]1[[1]]];
dh2 = DH2HomTrans[0,-\[Pi]/2,0,\[Theta]2[[1]]];
dh3 = DH2HomTrans[a2,0,d3,\[Theta]3[[1]]];

Hc4 = dh1.dh2.dh3;
Rc4 = Transpose[Hc4[[1;;3,1;;3]]].R;

If[Rc4[[2,3]]!=  1&&Rc4[[2,3]]!= 1,
\[Theta]54 = {ArcTan[Rc4[[2,3]],Sqrt[Rc4[[2,1]]^2+Rc4[[2,2]]^2]],ArcTan[Rc4[[2,3]],-Sqrt[Rc4[[2,1]]^2+Rc4[[2,2]]^2]]};
\[Theta]44 = ArcTan[-Rc4[[1,3]]/Sin[\[Theta]54],Rc4[[3,3]]/Sin[\[Theta]54]];
\[Theta]64 = ArcTan[Rc4[[2,1]]/Sin[\[Theta]54],-Rc4[[2, 2]]/Sin[\[Theta]54]];
];
If[Rc4[[2,3]] ==  1,
\[Theta]44 = {0,0};
\[Theta]54 = {0,0};
\[Theta]64 = {ArcTan[Rc4[[1,1]],-Rc4[[1,2]]],ArcTan[Rc4[[1,1]],-Rc4[[1,2]]]};
];
If[Rc4[[2,3]] ==  -1,
\[Theta]44 = {0,0};
\[Theta]54 = {\[Pi], \[Pi]};
\[Theta]64 = {-ArcTan[-Rc4[[1,1]],Rc4[[1,2]]],-ArcTan[-Rc4[[1,1]],Rc4[[1,2]]]};
];

Return[(180/(\[Pi]))*{{\[Theta]1[[2]], \[Theta]2[[4]], \[Theta]3[[2]], \[Theta]41[[1]], \[Theta]51[[1]], \[Theta]61[[1]]}, {\[Theta]1[[2]], \[Theta]2[[4]], \[Theta]3[[2]], \[Theta]41[[2]], \[Theta]51[[2]], \[Theta]61[[2]]},{\[Theta]1[[1]], \[Theta]2[[2]], \[Theta]3[[2]], \[Theta]42[[1]], \[Theta]52[[1]], \[Theta]62[[1]]},{\[Theta]1[[1]], \[Theta]2[[2]], \[Theta]3[[2]], \[Theta]42[[2]], \[Theta]52[[2]], \[Theta]62[[2]]}, {\[Theta]1[[2]], \[Theta]2[[3]], \[Theta]3[[1]], \[Theta]43[[1]], \[Theta]53[[1]], \[Theta]63[[1]]},{\[Theta]1[[2]], \[Theta]2[[3]], \[Theta]3[[1]], \[Theta]43[[2]], \[Theta]53[[2]], \[Theta]63[[2]]},{\[Theta]1[[1]], \[Theta]2[[1]], \[Theta]3[[1]], \[Theta]44[[1]], \[Theta]54[[1]], \[Theta]64[[1]]},{\[Theta]1[[1]], \[Theta]2[[1]], \[Theta]3[[1]], \[Theta]44[[2]], \[Theta]54[[2]], \[Theta]64[[2]]}}];
];


End[]
EndPackage[]
