main: main.cpp griditerator.cpp griditerator.h partitioninfo.h
	g++ -lstdc++ -Wall -Werror -g -O3  -MMD -MP griditerator.cpp main.cpp ../FreeImage/Dist/libfreeimage.a -I../FreeImage/Dist -o main