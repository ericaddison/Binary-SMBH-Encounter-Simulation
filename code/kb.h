//kbhit function for linux
//taken from http://linux-sxs.org/programming/kbhit.html

#ifndef KB_H
#define KB_H



void init_keyboard();
void close_keyboard();
int kbhit();
int getch();
void anykey();

#endif

