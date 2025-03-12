import sys

n = int(sys.argv[1])

file = open('blockMeshDict', 'w')

file.write(r"""FoamFile
{
	version		2.0;
	format		ascii;
	class		dictionary;
	object		blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 1;

x_min -16.0;
x_max  16.0;

y_min -16.0;
y_max  16.0;

z_min  0;
z_max  1;

n_x    """ + str(n) + """;
n_y    """ + str(n) + """;
n_z     1;

vertices
(
	($x_min $y_min $z_min)
	($x_max $y_min $z_min)
	($x_max $y_max $z_min)
	($x_min $y_max $z_min)
	($x_min $y_min $z_max)
	($x_max $y_min $z_max)
	($x_max $y_max $z_max)
	($x_min $y_max $z_max)
);

blocks
(
	hex	(0 1 2 3 4 5 6 7)	($n_x $n_y $n_z)	simpleGrading	(1 1 1)
);

edges
(
);

patches
(
	patch		front		(  (0 4 7 3)  )
	patch		rear		(  (1 2 6 5)  )
	empty		ground		(  (0 3 2 1)  )
	patch		left		(  (0 1 5 4)  )
	patch		right		(  (3 7 6 2)  )
	empty		top			(  (4 5 6 7)  )
);""")

file.close()

