#ifndef _UTIL_H
#define _UTIL_H

extern int check_pg_exist(string pg_name);
bool fexists(string filename);
extern void my_system(string cmd);
extern bool plsaMatrixExist(string data_path);
extern void show_version();

#endif /* _UTIL_H */
