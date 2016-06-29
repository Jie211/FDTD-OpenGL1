all:
	gcc FDTD2d_tm.c -framework GLUT -framework OpenGL -I/usr/X11/include -L/usr/X11/lib -lglut -lGLU -lGL -I/usr/local/include -L/usr/local/lib -Wall
