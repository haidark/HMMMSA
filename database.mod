GFORTRAN module version '0' created from database.f90 on Tue Dec 16 11:57:39 2014
MD5:8254035c3611acb579821e37e8536bf5 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () () ()
() () ())

()

()

()

()

(2 'database' 'database' 'database' 1 ((MODULE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ()
() () 0 0)
3 'database_build_fraglib' 'database' 'database_build_fraglib' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 4 0 (5 6 7 8) () 0 () () ()
0 0)
9 'database_dump_xyz' 'database' 'database_dump_xyz' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 UNKNOWN ()) 10 0 (11 12 13) () 0 () () () 0 0)
14 'database_free_bsf' 'database' 'database_free_bsf' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 UNKNOWN ()) 15 0 (16) () 0 () () () 0 0)
17 'database_read_bsf' 'database' 'database_read_bsf' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 UNKNOWN ()) 18 0 (19 20 21) () 0 () () () 0 0)
22 'frag' 'database' 'frag' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((
23 'istart' (DERIVED 24 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN POINTER) UNKNOWN-ACCESS ()) (25 'dme' (
REAL 4 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN) UNKNOWN-ACCESS ()) (26 'score' (REAL 4 0 0 REAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ())
(27 'next' (DERIVED 22 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN POINTER) UNKNOWN-ACCESS ())) PUBLIC (() ())
() 0 0)
28 'gamma_node' 'database' 'gamma_node' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0
0 () () 0 ((29 'i' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ()) (30 'g'
(REAL 4 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN) UNKNOWN-ACCESS ()) (31 'next' (DERIVED 28 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER)
UNKNOWN-ACCESS ())) PUBLIC (() ()) () 0 0)
24 'struct_node' 'database' 'struct_node' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0
0 () () 0 ((32 'gam' (DERIVED 28 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER) UNKNOWN-ACCESS ())
(33 'n_gamma' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ()) (34 'restype' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN)
UNKNOWN-ACCESS ()) (35 'dssp' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ()) (36 'delta' (REAL 4 0 0
REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN)
UNKNOWN-ACCESS ()) (37 'theta' (REAL 4 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ()) (38 'phi'
(REAL 4 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN) UNKNOWN-ACCESS ()) (39 'psi' (REAL 4 0 0 REAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ())
(40 'omg' (REAL 4 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) UNKNOWN-ACCESS ()) (41 'calpha' (REAL 4 0
0 REAL ()) (1 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '3')) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN DIMENSION) UNKNOWN-ACCESS ()) (42 'code' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '5')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN)
UNKNOWN-ACCESS ()) (43 'next' (DERIVED 24 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER)
UNKNOWN-ACCESS ())) PUBLIC (() ()) () 0 0)
19 'bsf' '' 'bsf' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN DUMMY) (
CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
20 'bsf_root' '' 'bsf_root' 18 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN POINTER DUMMY) (DERIVED 24 0 0 DERIVED ()) 0 0 () () 0 () () ()
0 0)
21 'nnodes' '' 'nnodes' 18 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
5 'bsf_root' '' 'bsf_root' 4 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
POINTER DUMMY) (DERIVED 24 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
6 'frags' '' 'frags' 4 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
ALLOCATABLE DIMENSION DUMMY) (DERIVED 22 0 0 DERIVED ()) 0 0 () (1
DEFERRED () ()) 0 () () () 0 0)
7 'nfrag' '' 'nfrag' 4 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
8 'frag_len' '' 'frag_len' 4 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
16 'root' '' 'root' 15 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
POINTER DUMMY) (DERIVED 24 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
11 'ifrag' '' 'ifrag' 10 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
POINTER DUMMY) (DERIVED 24 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
12 'xyz' '' 'xyz' 10 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
DIMENSION DUMMY) (REAL 4 0 0 REAL ()) 0 0 () (2 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0
0 INTEGER ()) 0 13 ())) 0 () () () 0 0)
13 'frag_len' '' 'frag_len' 10 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
)

('database' 0 2 'database_build_fraglib' 0 3 'database_dump_xyz' 0 9
'database_free_bsf' 0 14 'database_read_bsf' 0 17 'frag' 0 22 'gamma_node'
0 28 'struct_node' 0 24)
