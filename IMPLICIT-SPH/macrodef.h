#ifndef __MACRODEF_H__
#define __MACRODEF_H__

using namespace std;

#define ISPH_EPSILON 1.0e-24

#define ISPH_MLS_MAX_ORDER 10
#define ISPH_BDF_MAX_ORDER 4
#define ISPH_MAX_CONCENTRATION 4

#define ISPH_COMPUTE_STATUS_MAX_SIZE 20

#undef  FALSE
#define FALSE 0
  
#undef  TRUE
#define TRUE 1

#undef  LAMMPS_FAILURE
#define LAMMPS_FAILURE -1
  
#undef  LAMMPS_SUCCESS
#define LAMMPS_SUCCESS 0
  
#undef  __FUNCT__
#define __FUNCT__ "NONAME"

#ifdef DEBUG_VERBOSE
#define FUNCT_ENTER(me) if (me == 0)                                    \
    cout << (unsigned long)this                                         \
         << setw(10) << right << string(++g_funct_counter, '+')  << ", " \
         << __FUNCT__ << endl;
#define FUNCT_EXIT(me)  if (me == 0)                                    \
    cout << (unsigned long)this                                         \
         << setw(10) << right << string(g_funct_counter--, '-')  << ", " \
         << __FUNCT__ << endl;
#else
#define FUNCT_ENTER
#define FUNCT_EXIT
#endif 

#define QUICK_RETURN(cond, r_val) if ((cond)) { return r_val; } 

#ifdef ERROR_CHECK
#define IS_VALID(cond, pid, what)               \
  if (!(cond)) {                                \
    if (pid == 0) cout << ERROR(what) << endl;  \
    throw -1;                                   \
  }
#else
#define IS_VALID(cond, pid, what) 
#endif
  
#undef  ERROR
#define ERROR(what) what << ", "<< __FILE__ << ", " << __LINE__
  
#undef  LOOKUP
#define LOOKUP(A,i,j) A[i][j]
  
#undef  VIEW2
#define VIEW2(A, dim, i, j) A[j*dim+i]

#endif
