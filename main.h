/*
 *  main.h
 *  Assignment One
 *
 *  Created by Derek Williams on 10-01-30.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#define FORWARD_SPACE		0
#define BACKWARD_SPACE		1
#define CENTRAL_SPACE		2
#define	LAX_FREDRICHS		3
#define LEAPFROG			4
#define EQUILIBRIUM			5
#define LAX_WENDROFF		6
#define CRANK_NICOLSON		7
#define NUM_SCHEMES			8

#define FORWARD_SPACE_NAME	"Forward-Time Forward-Space"
#define BACKWARD_SPACE_NAME	"Forward-Time Backward-Space"
#define CENTRAL_SPACE_NAME	"Forward-Time Central-Space"
#define	LAX_FREDRICHS_NAME	"Lax Fredrichs"
#define LEAPFROG_NAME		"Leapfrog"
#define EQUILIBRIUM_NAME	"Equilibrium"
#define LAX_WENDROFF_NAME	"Lax Wendroff"
#define CRANK_NICOLSON_NAME	"Crank Nicolson"

#define A_CONST				1000
#define A_SPACE				1001

#define A_CONST_NAME		"constant"
#define A_SPACE_NAME		"1 + 1/4(3 - x)(1 + x)"

#define INIT_ONE			0
#define INIT_TWO			1
#define INIT_THREE			2
#define INIT_FOUR			3
#define INIT_FIVE			4
#define INIT_SIX			5
#define INIT_SEVEN			6
#define INIT_EIGHT			7
#define NUM_INIT			8

#define INIT_ONE_NAME		"cos^2(π x), |x| <= 1/2"
#define INIT_TWO_NAME		"1 - |x|, |x| <= 1"
#define INIT_THREE_NAME		"sin(2π x)"
#define INIT_FOUR_NAME		"cos(2π x)"
#define INIT_FIVE_NAME		"cos(ξπ x) cos^2(π/2 x)"
#define INIT_SIX_NAME		"1 - |x| if |x| < 1/2, 1/4 if |x| = 1/2"
#define INIT_SEVEN_NAME		"sin(1.2(x - y))cosh(x + 2y)"
#define INIT_EIGHT_NAME		"exp(-(x^2 + y^2)"

#define MAX_VM_N			3000
#define SAVE_ID				3001

#define F_ZERO				4000
#define F_ONE				4001

#define F_ZERO_NAME			"F(t, x) = 0"
#define F_ONE_NAME			"F(t, x) = sin(π (x-t))"

#define SETTINGS_ID			5000

#define BOUND_ONE			0
#define BOUND_TWO			1
#define BOUND_THREE			2
#define BOUND_FOUR			3
#define BOUND_FIVE			4
#define NUM_BOUND			5

#define BOUND_ONE_NAME		"v(n+1,M) = v(n,M-1)"
#define BOUND_TWO_NAME		"v(n+1,M) - v(n,M) + λ(v(n+1,M) - v(n+1,M-1)) = kv(n+1,M)"
#define BOUND_THREE_NAME	"v(n+1,M) = 2v(n+1,M-1) - v(n+1,M-2)"
#define BOUND_FOUR_NAME		"Neumann"
#define BOUND_FIVE_NAME		"Exact"

#define DERIV_VANILLA_CALL	0
#define DERIV_VANILLA_PUT	1
#define NUM_DERIVS			2

#define DERIV_VANILLA_CALL_NAME	"Vanilla Call Option"
#define DERIV_VANILLA_PUT_NAME	"Vanilla Put Option"

#define TIME_START			0

#define BLOWUP_VAL			5.0

#define ALIGN_CENTER		0
#define ALIGN_RIGHT			1
#define ALIGN_LEFT			-1

#define TINY				-1e100

#define BUFLEN				1024

#define CONVERG				1e-10
#define	SUM_MAX				1e6

#define SCHEME_CTRL			1
#define H_CTRL				2
#define TIME_CTRL			3
#define UTERM_CTRL			4
#define LAMBDA_CTRL			5
#define XI_CTRL				6
#define IC_CTRL				7
#define TOP_CTRL			8
#define BOTTOM_CTRL			9
#define RIGHT_CTRL			10
#define LEFT_CTRL			11
#define REPEAT_CTRL			12
#define BC_CTRL				13
#define SCALE_CTRL			14
#define MU_CTRL				15
#define A_CTRL				16
#define B_CTRL				17
#define DIV_CTRL			18
#define SAVE_CTRL			19
#define LOAD_CTRL			20
#define S_CTRL				21
#define R_CTRL				22
#define K_CTRL				23
#define DERIV_CTRL			24

#define DEFAULT_A			&aConst
#define DEFAULT_B			&bConst

#define DEFAULT_LAMBDA		0.8

#define BOUNDS_BOTTOM		0
#define BOUNDS_TOP			100

#define	DEFAULT_THREEDIM	0		

#define TIME_END			1

#define BOUNDS_LEFT			0
#define BOUNDS_RIGHT		100

#define DEFAULT_SCHEME	BACKWARD_SPACE
#define DEFAULT_BOUNDARY &BZero
#define DEFAULT_INIT		&initOne
#define DEFAULT_MU		0.0

#define DEFAULT_INIT_NUM INIT_ONE
#define DEFAULT_BOUND	8


#define CF_SCHEME		0
#define CF_LAMBDA		1
#define CF_MU			2
#define CF_TOP			3
#define CF_BOTTOM		4
#define CF_LEFT			5
#define CF_RIGHT		6
#define CF_BOUNDARY		7
#define CF_TIME			8
#define CF_THREEDIM		9
#define CF_WIRE			10
#define CF_ACONST		11
#define CF_BCONST		12
#define CF_A_FUNC		13
#define CF_B_FUNC		14
#define CF_INIT_FUNC	15
#define CF_BOUNDS_FUNC	16
#define CF_H			17
#define CF_PBC			18
#define CF_RHO			19
#define CF_THETA		20
#define CF_PHI			21
#define CF_R			22
#define CF_SIGMA		23
#define CF_K			24
#define CF_DERIV		25

/* handy vector operations */

#define VSUB3(A, B, C)				\
	A[0] = B[0] - C[0];				\
	A[1] = B[1] - C[1];				\
	A[2] = B[2] - C[2];

#define VCROSS3(A, B, C)			\
	A[0] =   B[1]*C[2] - B[2]*C[1];	\
	A[1] = - B[0]*C[2] + B[2]*C[0];	\
	A[2] =   B[0]*C[1] - B[1]*C[0];

#define VEQ3(A, B)					\
	A[0] = B[0];					\
	A[1] = B[1];					\
	A[2] = B[2];

#define VSET2(A, X, Y)				\
	A[0] = X;						\
	A[1] = Y;

#define VSET3(A, X, Y, Z)			\
	A[0] = X;						\
	A[1] = Y;						\
	A[2] = Z;

