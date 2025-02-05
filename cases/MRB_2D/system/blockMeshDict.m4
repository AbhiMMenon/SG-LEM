/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])

convertToMeters 1e-3;
// sizes go here
m4_define(cellSize, 1.0)  // mm

m4_define(D1, 2.6)
m4_define(D2, 7)
m4_define(D3, 40)
m4_define(D4, 60)
m4_define(D5, 80)
m4_define(D6, 400)
m4_define(L, 250)
m4_define(alpha, 26)
m4_define(zmax, 5)
m4_define(k, 10)

m4_define(r1, calc(D1/2))
m4_define(r2, calc(D2/2))
m4_define(r3, calc(D3/2))
m4_define(r4, calc(D4/2))
m4_define(r5, calc(D5/2))
m4_define(r6, calc(D6/2))

m4_define(p3, calc(0-(r3-r2)*tan(deg2rad(26))))
m4_define(p4, calc(p3-(r4-r3)*tan(deg2rad(26))))
m4_define(p5, calc(p4-(r5-r4)*tan(deg2rad(26))))
m4_define(p6, calc(p4-k))


//m4_define(LCell,calc(int(L/cellSize)))
//m4_define(r1Cell,calc(int(r1/cellSize*3)))
//m4_define(r2Cell,calc(int(r2/cellSize*2)))
//m4_define(r3Cell,calc(int(r3/cellSize)))
//m4_define(r4Cell,calc(int(r4/cellSize*0.6)))

m4_define(LCell,calc(int(L/cellSize)))
m4_define(r1Cell,20)
m4_define(r2Cell,20)
m4_define(r3Cell,80)
m4_define(r4Cell,40)
m4_define(r5Cell,40)
m4_define(r6Cell,40)

//m4_define(zCell,calc(int(zmax/cellSize)))
m4_define(zCell,10)

vertices
(

    (0 0  0)    //0
    (L 0  0)    //1
    (0 r1 0)    //2
    (L r1 0)    //3
    (0 r2 0)    //4
    (L r2 0)    //5
    (p3 r3 0)    //6
    (L r3 0)    //7

    (0 0  zmax)    //8
    (L 0  zmax)    //9
    (0 r1 zmax)    //10
    (L r1 zmax)    //11
    (0 r2 zmax)    //12
    (L r2 zmax)    //13
    (p3 r3 zmax)    //14
    (L r3 zmax)    //15
    
    (p4 r4    0)    //16
    (L r4     0)    //17
    (p4 r4 zmax)    //18
    (L r4  zmax)    //19
    
    (p5 r5    0)    //20
    (L  r5     0)   //21
    (p5 r5 zmax)    //22
    (L  r5  zmax)   //23
    
    (p5 r6    0)    //24
    (L  r6     0)   //25
    (p5 r6 zmax)    //26
    (L  r6  zmax)   //27
    
    (p6 r5    0)    //28
    (p6 r6    0)    //29
    (p6 r5    zmax)    //30
    (p6 r6    zmax)    //31
    
    (p6 r3    0)    //32
    (p6 r4    0)    //23
    (p6 r3    zmax)    //34
    (p6 r4    zmax)    //35


);

blocks
(
    hex (0 1 3 2 8 9 11 10)(LCell r1Cell zCell) simpleGrading (1 1 1) // 0
    hex (2 3 5 4 10 11 13 12)(LCell r2Cell zCell) simpleGrading (1 1 1) // 1
    hex (4 5 7 6 12 13 15 14)(LCell r3Cell zCell) simpleGrading (1 1 1) // 2
    hex (6 7 17 16 14 15 19 18)(LCell r4Cell zCell) simpleGrading (1 1 1) // 3
    hex (16 17 21 20 18 19 23 22)(LCell r5Cell zCell) simpleGrading (1 1 1) // 4
    hex (20 21 25 24 22 23 27 26)(LCell r6Cell zCell) simpleGrading (1 40 1) // 5
    hex (28 20 24 29 30 22 26 31)(10 r6Cell zCell) simpleGrading (1 40 1) // 
    hex (32 6 16 33 34 14 18 35)(40 r4Cell zCell) simpleGrading (1 1 1) // 


);

boundary
(
  INLET_JET
  {
  	type patch;
  	faces
  	(
  	    (0 8 10 2)
  	    
  	);
  }
  
  INLET_SLOT_1
  {
  	type patch;
  	faces
  	(
  	    (2 10 12 4)
  	    
  	);
  }
  
 INLET_SLOT_2
  {
  	type patch;
  	faces
  	(
  	    (32 34 35 33)
  	    
  	);
  }
  
 INLET_COFLOW
  {
  	type patch;
  	faces
  	(
  	    (28 30 31 29)
  	    
  	);
  }
 
WALL_BB_1
  {
  	type wall;
  	faces
  	(
  	    (4 12 14 6)
  	    (6 14 34 32)
  	    
  	);
  }

WALL_BB_2
  {
  	type wall;
  	faces
  	(
  	    (33 35 16 18)
  	    (16 18 22 20)
  	    (28 20 22 30)
  	    
  	);
  }

SIDE
  {
  	type wall;
  	faces
  	(
  	    (29 24 26 31)
  	    ( 24 25 27 26)
  	    
  	);
  }

OUTLET
  {
  	type patch;
  	faces
  	(
  	    (1 9 11 3 )
  	    ( 3 11 13 5)
  	    ( 5 13 15 7)
  	    (7 15 19 17)
  	    (17 19 23 21)
  	    (21 23 27 25)
  	    
  	);
  }

BACK
  {
  	type cyclic;
    neighbourPatch FRONT;
  	faces
  	(
  	    (0 2 3 1)
  	    (2 4  5 3)
  	    ( 4 6 7 5)
  	    (6 16 17 7)
  	    (32 33 16 6)
  	    (16 20 21 17)
  	    (20 24 25 21)
  	    (28 29 24 20)
  	    
  	);
  }

FRONT
  {
  	type cyclic;
    neighbourPatch BACK;
  	faces
  	(
  	    (8 9 11 10)
  	    ( 10 11 13 12)
  	    ( 12 13 15 14)
  	    (14 15 19 18 )
  	    (34 14 18 35)
  	    ( 18 19 23 22)
  	    ( 22 23 27 26)
  	    (30 22 26 31)
  	    
  	);
  }

AXIS
  {
  	type symmetryPlane;
  	faces
  	(
  	    (8 9 1 0)
  	    
  	);
  }



);

mergePatchPairs
(
);

// ************************************************************************* //
